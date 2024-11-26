#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

namespace nb = nanobind;

#ifndef SOFTCORR_H
#define SOFTCORR_H

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
    idbuf = shmget(DasBufferKey, sizeof(Buffer), SHM_RDONLY);
    
    if (idbuf < 0) {
        printf("Shared memory does not exist.\n");
        return nb::ndarray<uint8_t>();
    }
    
    BufRead = (Buffer *)shmat(idbuf, nullptr, SHM_RDONLY);
    if (BufRead == (Buffer *)-1) {
        printf("Could not attach to shared memory.\n");
        return nb::ndarray<uint8_t>();
    }
    
    // Define block and bin locations for the beginning and end of data
    int binbeg = offset;
    int binend = count + offset; 
    int binbeg_block_loc = binbeg / total_bin_in_FRBblock; //This will increment forever during the observation period and it is important for checking if SHM is overwritten
    int binbeg_block_loc_cycle = binbeg_block_loc / 12 ; // This will cycle between 0 and 11
    int binbeg_bin_loc = binbeg % total_bin_in_FRBblock;
    int binend_block_loc = binend / total_bin_in_FRBblock;
    int binend_block_loc_cycle = binend_block_loc / 12 ; 
    int binend_bin_loc = binend % total_bin_in_FRBblock;
    nBeams = BufRead->nBeams;
    RecNum = dataBuffer_FRB->curRecord;


if(0 =< RecNum - binbeg_block_loc < 12 ){
    // Calculate total size of data to retrieve
    size_t total_size = 0;
    for (int block = binbeg_block_loc_cycle; block <= binend_block_loc_cycle; ++block) {
        if (block == binbeg_block_loc_cycle) {
            total_size += DataSize - bin_size * binbeg_bin_loc_cycle; // Start block
        } else if (block == binend_block_loc_cycle) {
            total_size += bin_size * binend_bin_loc_cycle; // End block
        } else {
            total_size += DataSize; // Full intermediate blocks
        }
    }

    // Allocate a temporary buffer to hold the data
    uint8_t *buffer = (uint8_t *)malloc(total_size);
    if (!buffer) {
        printf("Memory allocation failed.\n");
        shmdt(BufRead);
        return nb::ndarray<uint8_t>();
    }

    // Copy data into buffer
    size_t offset = 0;
    for (int block = binbeg_block_loc_cycle; block <= binend_block_loc_cycle; ++block) {
        size_t segment_size;
        if (block == binbeg_block_loc_cycle) {
            segment_size = DataSize - bin_size * binbeg_bin_loc_cycle;
            memcpy(buffer + offset, BufRead->data + (long)DataSize * block * NBeams + (long)binbeg_bin_loc_cycle * bin_size, segment_size);
        } else if (block == binend_block_loc_cycle) {
            segment_size = bin_size * binend_bin_loc_cycle;
            memcpy(buffer + offset, BufRead->data + (long)DataSize * block * NBeams, segment_size);
        } else {
            segment_size = DataSize;
            memcpy(buffer + offset, BufRead->data + (long)DataSize * block * NBeams, segment_size);
        }
        offset += segment_size;
    }



    // Determine dimensions for reshaping: here, `nf` (frequency channels) is defined as NCHANNELS
    size_t nf = NCHANNELS;
    size_t num_samples = total_size / nf;

    // Wrap the buffer as a NumPy array with specified shape
    nb::ndarray<uint8_t> result = nb::ndarray<uint8_t>(buffer, {nf, num_samples}, nb::capsule(buffer, free));

    return result;
}
else{printf("Your data has been overwritten already")}  //Replace by better error messege

NB_MODULE(your_module_name, m) {
    m.def("get_data_as_numpy_array", &get_data_as_numpy_array, "Retrieve data from shared memory as a NumPy array");
}
}
