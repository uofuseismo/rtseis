# Installation Instructions

Below are instructions for building RTSeis from source on Linux.

## Dependencies

Prior to building RTSeis the following dependencies must be cleared:

   1. A C11 or higher C compiler.
   2. A C++11 or higher C++ compiler.
   3. [CMake](https://cmake.org/) V2.8 or higher.
   4. [Intel's Performance Primitives](https://software.intel.com/en-us/intel-ipp) performs most of the heavy lifting.
   5. [Intel's MKL](https://software.intel.com/en-us/mkl) backfills some of the matrix algebra and any massiv correlations.
   6. [cJSON](https://github.com/DaveGamble/cJSON) parses JSON messages.
   7. Python3 with NumPy and ctypes.

## Configuring

After clearing the dependencies configure CMake

