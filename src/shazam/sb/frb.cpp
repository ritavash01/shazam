#include <cstdint>
#include <cstring>
#include <sys/shm.h>
#include <sys/types.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

#define MAXBLKS 16
#define HDRKEY 2031
#define BUFKEY 2032
#define BLKSIZE (32 * 512 * 4096)
#define TOTALSIZE (long)(BLKSIZE) * (long)(MAXBLKS)

/* Struct for storing the ring buffer's header. */
typedef struct {
  unsigned int active;
  unsigned int status;
  double comptime;
  double datatime;
  double reftime;
  struct timeval timestamp[MAXBLKS];
  struct timeval timestamp_gps[MAXBLKS];
  double blk_nano[MAXBLKS];
} Header;

/* Struct for storing data from the ring buffer. */
typedef struct {
  unsigned int flag;
  unsigned int curr_blk;
  unsigned int curr_rec;
  unsigned int blk_size;
  int overflow;
  double comptime[MAXBLKS];
  double datatime[MAXBLKS];
  unsigned char data[TOTALSIZE];
} Buffer;

NB_MODULE(internals, m) {
  m.def(
      "_readfrbshm",
      [](int offset, int nrecords) {
        int idHdr = shmget(HDRKEY, sizeof(Header), SHM_RDONLY);
        int idBuf = shmget(BUFKEY, sizeof(Buffer), SHM_RDONLY);
        if (idHdr < 0 || idBuf < 0) exit(1);

        Header *Hdr = (Header *)shmat(idHdr, 0, 0);
        Buffer *Buf = (Buffer *)shmat(idBuf, 0, 0);
        if ((Buf) == (Buffer *)-1) exit(1);

        size_t nf = 4096;
        size_t nt = static_cast<size_t>(BLKSIZE * nrecords / nf);
        size_t shape[2] = {nt, nf};
        return nb::ndarray<nb::numpy, std::uint8_t, nb::shape<-1, -1>>(
            Buf->data + (long)BLKSIZE * (long)offset, 2, shape, nb::handle());
      },
      "Read from the pipeline's ring buffer, FRB_SHM.\n"
      "\n"
      "SHARED MEMORY SHENANIGANS!\n"
      "==========================\n"
      "\n"
      "For a sampling time of 1.31072ms, the shared memory at the\n"
      "telescope is structured as 32 blocks, with each block being\n"
      "512 samples, or 0.67108864s, long. The entire shared memory\n"
      "at the telescope is 21.4748364s long, with a size of 64 MB.\n"
      "We then form another shared memory when we wish to search for\n"
      "FRBs, where each block is 21.47483648s long, and there are\n"
      "16 blocks. This makes this shared memory 343.59738368s long,\n"
      "with a size of 1 GB.");
  m.def("_streamfrbshm", []() {});
}
