#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

namespace nb = nanobind; 
nb::ndarray<uint8_t> get_data_as_numpy_array(int count, int offset); 

