git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build
module laod profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
cmake ..
make -j
