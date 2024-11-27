#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>
#include <unistd.h>

namespace nb = nanobind;

#define DasBufferKey 2032
#define NCHANNELS 4096
#define NParallel 5
#define NSerial 2
#define BITREDUCTION 2
#define NBeams (NParallel * NSerial)
#define FFT_Samps_per_Block 800
#define TEL_SHM_BlockSize (FFT_Samps_per_Block * NCHANNELS)
#define TEL_to_FRB_Block_factor 32
#define DataSize (TEL_to_FRB_Block_factor * TEL_SHM_BlockSize)
#define MaxDataBlocks 12
#define total_bin_in_FRBblock 25600
#define bin_size 0.00131072 // Size of each bin in bytes

typedef struct {
    unsigned int active, status, is_buf_empty;
    double pc_time, ref_time, rec_time;
    struct timeval timestamp_gps[MaxDataBlocks];
    double blk_nano[MaxDataBlocks];
    unsigned int flag, curBlock, curRecord, blockSize, nBeams;
    int overFlow;
    unsigned char data[(long)NBeams * (long)(DataSize) * (long)(MaxDataBlocks)];
} Buffer;

nb::ndarray<uint8_t> get_data_as_numpy_array(int count, int offset) {
    // Attach to shared memory
    int idbuf = shmget(DasBufferKey, sizeof(Buffer), SHM_RDONLY);
    if (idbuf < 0) {
        throw std::runtime_error("Shared memory does not exist.");
    }

    Buffer* BufRead = static_cast<Buffer*>(shmat(idbuf, nullptr, SHM_RDONLY));
    if (BufRead == reinterpret_cast<Buffer*>(-1)) {
        throw std::runtime_error("Could not attach to shared memory.");
    }

    // Calculate block and bin locations
    int binbeg = offset;
    int binend = count + offset;
    int binbeg_block_loc = binbeg / total_bin_in_FRBblock;
    int binbeg_block_loc_cycle = binbeg_block_loc % MaxDataBlocks; // Cycle within available blocks
    int binbeg_bin_loc = binbeg % total_bin_in_FRBblock;
    int binend_block_loc = binend / total_bin_in_FRBblock;
    int binend_block_loc_cycle = binend_block_loc % MaxDataBlocks; // Cycle within available blocks
    int binend_bin_loc = binend % total_bin_in_FRBblock;

    int nBeams = BufRead->nBeams;
    size_t total_data_size = 0;

    // Calculate total size for all beams
    for (int beam = 0; beam < nBeams; ++beam) {
        if (binbeg_block_loc_cycle == binend_block_loc_cycle) {
            total_data_size += bin_size * (binend_bin_loc - binbeg_bin_loc);
        } else {
            for (int block = binbeg_block_loc_cycle; block <= binend_block_loc_cycle; ++block) {
                if (block == binbeg_block_loc_cycle) {
                    total_data_size += DataSize - bin_size * binbeg_bin_loc;
                } else if (block == binend_block_loc_cycle) {
                    total_data_size += bin_size * binend_bin_loc;
                } else {
                    total_data_size += DataSize;
                }
            }
        }
    }
        
    // Allocate temporary buffer
    uint8_t* temp_buffer = static_cast<uint8_t*>(malloc(total_data_size));
    if (!temp_buffer) {
        throw std::bad_alloc();
    }

    size_t global_offset = 0;

    // Extract data for each beam
    for (int beam = 0; beam < nBeams; ++beam) {
        if (binbeg_block_loc_cycle == binend_block_loc_cycle) {
            size_t segment_size = bin_size * (binend_bin_loc - binbeg_bin_loc);
            memcpy(temp_buffer + global_offset,
                   BufRead->data + (long)DataSize * binbeg_block_loc_cycle * nBeams + (long)DataSize * beam + (long)bin_size * binbeg_bin_loc,
                   segment_size);
            global_offset += segment_size;
        } else {
            for (int block = binbeg_block_loc_cycle; block <= binend_block_loc_cycle; ++block) {
                if (block == binbeg_block_loc_cycle) {
                    size_t segment_size = (long)DataSize - bin_size * binbeg_bin_loc;
                    memcpy(temp_buffer + global_offset,
                           BufRead->data + (long)DataSize * block * nBeams + (long)DataSize * beam + (long)bin_size * binbeg_bin_loc,
                           segment_size);
                    global_offset += segment_size;
                } else if (block == binend_block_loc_cycle) {
                    size_t segment_size = bin_size * binend_bin_loc;
                    memcpy(temp_buffer + global_offset,
                           BufRead->data + (long)DataSize * block * nBeams + (long)DataSize * beam,
                           segment_size);
                    global_offset += segment_size;
                } else {
                    size_t segment_size = (long)DataSize;
                    memcpy(temp_buffer + global_offset,
                           BufRead->data + (long)DataSize * block * nBeams + (long)DataSize * beam,
                           segment_size);
                    global_offset += segment_size;
                }
            }
        }
    }

    // Create ndarray
    size_t num_samples = total_data_size / (NCHANNELS * nBeams);
    nb::capsule free_when_done(temp_buffer, [](void* p) noexcept { free(p); });
    return nb::ndarray<uint8_t>(temp_buffer, {nBeams, NCHANNELS, num_samples}, free_when_done);
}
