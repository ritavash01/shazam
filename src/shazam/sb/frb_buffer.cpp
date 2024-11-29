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
#define NSerial 10
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

nb::ndarray<nb::numpy, uint8_t, nb::shape<-1, -1, -1>> get_data_as_numpy_array(int count, int offset) {
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
   

    // Calculate total size for all beams (This will give me the size of the temporary buffer to whole the total data)
 
    size_t total_data_size_per_beam = 0;
    if (binbeg_block_loc_cycle == binend_block_loc_cycle) {
        total_data_size_per_beam += bin_size * (binend_bin_loc - binbeg_bin_loc);
    } else {
        for (int block = binbeg_block_loc_cycle; block <= binend_block_loc_cycle; ++block) {
            if (block == binbeg_block_loc_cycle) {
                total_data_size_per_beam += DataSize - bin_size * binbeg_bin_loc;
            } else if (block == binend_block_loc_cycle) {
                total_data_size_per_beam += bin_size * binend_bin_loc;
            } else {
                total_data_size_per_beam += DataSize;
            }
        }
    }

   //Reshaping the data
   size_t nt = total_data_size_per_beam / NCHANNELS; // Number of time samples
    if (nt * NCHANNELS != total_data_size_per_beam) {
        throw std::runtime_error("Data size is not evenly divisible by nf");
    }
   // Total data size for all beams
    size_t total_data_size = nBeams * total_data_size_per_beam;

// Allocate buffer for all beams
uint8_t* temp_buffer = (uint8_t*)malloc(total_data_size);
if (!temp_buffer) {
    throw std::runtime_error("Memory allocation failed");
}

// Extract data for each beam and store in its section of the buffer
size_t global_offset = 0;
for (int beam = 0; beam < nBeams; ++beam) {
    size_t beam_offset = global_offset;

    if (binbeg_block_loc_cycle == binend_block_loc_cycle) {
        size_t segment_size = bin_size * (binend_bin_loc - binbeg_bin_loc);
        memcpy(temp_buffer + beam_offset,
               BufRead->data + (long)DataSize * binbeg_block_loc_cycle * nBeams + (long)DataSize * beam + (long)bin_size * binbeg_bin_loc,
               segment_size);
        global_offset += segment_size;
    } else {
        for (int block = binbeg_block_loc_cycle; block <= binend_block_loc_cycle; ++block) {
            if (block == binbeg_block_loc_cycle) {
                size_t segment_size = (long)DataSize - bin_size * binbeg_bin_loc;
                memcpy(temp_buffer + beam_offset,
                       BufRead->data + (long)DataSize * block * nBeams + (long)DataSize * beam + (long)bin_size * binbeg_bin_loc,
                       segment_size);
                beam_offset += segment_size;
            } else if (block == binend_block_loc_cycle) {
                size_t segment_size = bin_size * binend_bin_loc;
                memcpy(temp_buffer + beam_offset,
                       BufRead->data + (long)DataSize * block * nBeams + (long)DataSize * beam,
                       segment_size);
                beam_offset += segment_size;
            } else {
                size_t segment_size = (long)DataSize;
                memcpy(temp_buffer + beam_offset,
                       BufRead->data + (long)DataSize * block * nBeams + (long)DataSize * beam,
                       segment_size);
                beam_offset += segment_size;
            }
        }
    }
}


 // Define shape for 3D ndarray
    size_t shape[3] = {nBeams, nt, NCHANNELS};

    // Create the ndarray
    nb::capsule free_when_done(temp_buffer, [](void* p) noexcept { free(p); });
    return nb::ndarray<nb::numpy, std::uint8_t, nb::shape<-1, -1, -1>>(
        temp_buffer, 3, shape, free_when_done);
}
