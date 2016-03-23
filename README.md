# fastLUT
Fast lookup table C++ code with linear interpolation, GPU acceleration and z-curve caching optimizations.

Currently, the script go.sh will compile a shared library for you to link with your code. The executable "test" will run some benchmarks to show how it can be used. 

# CPU mode
The 3D data can be arranged according to a z-curve algorithm (morton numbers) to preserve as much locality as possible, and thus optimize the CPU cache. This only really becomes worthwhile with large arrays.

# GPU mode
The 3D data is stored on the GPU VRAM and accessed via the texture memory hardware. This will perform fast linear interpolation but with limited precision (enough for most use cases).