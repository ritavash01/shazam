#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
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

/* Structure that describes each shared memory block */
typedef struct {
    unsigned int active, status, is_buf_empty;
    double pc_time, ref_time, rec_time;
    struct timeval timestamp_gps[MaxDataBlocks];
    double blk_nano[MaxDataBlocks];
    unsigned int flag, curBlock, curRecord, blockSize, nBeams;
    int overFlow;
    unsigned char data[(long)NBeams * (long)(DataSize) * (long)(MaxDataBlocks)];
} Buffer;

// Function to retrieve data as a NumPy array
nb::ndarray<uint8_t> get_data_as_numpy_array(int count, int offset) {
    int idbuf;
    Buffer *BufRead;

    // Attach to the shared memory
    idbuf = shmget(DasBufferKey, sizeof(Buffer), SHM_RDONLY);
    if (idbuf < 0) {
        throw std::runtime_error("Shared memory does not exist.");
    }
    
    BufRead = (Buffer *)shmat(idbuf, nullptr, SHM_RDONLY);
    if (BufRead == (Buffer *)-1) {
        throw std::runtime_error("Could not attach to shared memory.");
    }

    // Calculate bin and block positions
    int binbeg = offset;
    int binend = count + offset; 
    int binbeg_block_loc = binbeg / total_bin_in_FRBblock;
    int binbeg_bin_loc = binbeg % total_bin_in_FRBblock;
    int binend_block_loc = binend / total_bin_in_FRBblock;
    int binend_bin_loc = binend % total_bin_in_FRBblock;

    if (BufRead->curRecord - binbeg_block_loc < 0 || BufRead->curRecord - binbeg_block_loc >= MaxDataBlocks) {
        shmdt(BufRead);
        throw std::runtime_error("Data has been overwritten in shared memory.");
    }

    // Calculate total size of data to retrieve
    size_t total_size = (binend - binbeg) * bin_size;
    uint8_t *buffer = (uint8_t *)malloc(total_size);
    if (!buffer) {
        shmdt(BufRead);
        throw std::runtime_error("Memory allocation failed.");
    }

    // Copy data from shared memory
    size_t offset_ptr = 0;
    for (int block = binbeg_block_loc; block <= binend_block_loc; ++block) {
        size_t start = (block == binbeg_block_loc) ? binbeg_bin_loc * bin_size : 0;
        size_t end = (block == binend_block_loc) ? binend_bin_loc * bin_size : DataSize;
        size_t segment_size = end - start;

        memcpy(buffer + offset_ptr, BufRead->data + (long)block * DataSize + start, segment_size);
        offset_ptr += segment_size;
    }


    

    // Wrap buffer in an ndarray
    size_t nf = NCHANNELS;
    size_t num_samples = total_size / nf;
    nb::capsule free_when_done(buffer, [](void *p) { free(p); });

    return nb::ndarray<uint8_t>(buffer, {nf, num_samples}, free_when_done);
}

NB_MODULE(internals, m) {
    m.def("get_data_as_numpy_array", &get_data_as_numpy_array, "Retrieve data from shared memory as a NumPy array");
}
