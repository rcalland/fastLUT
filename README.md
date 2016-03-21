# fastLUT
Fast lookup table C++ code with linear interpolation and z-curve caching optimizations.

Currently, the script go.sh will compile a shared library for you to link with your code. The executable "test" will run some benchmarks to show how it can be used. The code requires the cern library ROOT to compile, however it is not fundamental to the code, and will probably be removed later. 

The 3D data will be arranged according to a z-curve algorithm (morton numbers) to preserve as much locality as possible, and thus optimize the CPU cache. This only really becomes worthwhile with large arrays.
