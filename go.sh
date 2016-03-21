FLAGS="-O2 -fno-align-functions -fno-align-loops -mtune=native"

g++ ${FLAGS} `root-config --glibs --cflags` -std=c++11 -fpic -Wall fastLUT3D.cc -c
g++ ${FLAGS} `root-config --glibs --cflags` -shared -o libfastLUT3D.so fastLUT3D.o

g++ ${FLAGS} -L/home/rcalland/work/fastLUT -std=c++11 -Wall -o test test.cc -lfastLUT3D `root-config --glibs --cflags`
