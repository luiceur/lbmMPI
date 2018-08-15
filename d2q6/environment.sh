# Shell script to contain environment settings for various systems
# Usage: 'source environment.sh <system>'
SYSTEM=$1
case $SYSTEM in
idun)
    module load intel
    mkdir -p obj
    export CC=mpiicc
    export ARCH=intel
    export FFMPEG=${HOME}/tools/bin/ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -b:v 16384k"
    ;;
epic)
    module load GCC/5.4.0-2.26
    module load OpenMPI/1.10.3
    mkdir -p obj
    export PATH+=:/usr/local/cuda/bin
    export NVCC=nvcc
    export CC=nvcc
    export ARCH=cuda
    export FFMPEG=${HOME}/tools/bin/ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -b:v 16384k"
    ;;
local)
    mkdir -p obj
    export CC=mpicc
    export ARCH=generic
    export FFMPEG=ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -tune grain -b:v 16384k"
    ;;
*)
    echo "Environment not predefined for system '${SYSTEM}'"
    ;;
esac
