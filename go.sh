FLAGS="-O2" # -mtune=native" # -DHEMI_CUDA_DISABLE"

nvcc ${FLAGS} -Xcompiler -fPIC -std=c++11 fastLUT3D.cc -c
g++ ${FLAGS} -shared -o libfastLUT3D.so fastLUT3D.o

#g++ ${FLAGS} -L. -std=c++11 -Wall -o test test.cc -lfastLUT3D `root-config --glibs --cflags`
nvcc -L. -ccbin g++ -std=c++11 -o test_gpu test_gpu.cc -L/usr/local/cuda-7.5/lib64 -lcudart -lfastLUT3D --expt-extended-lambda
