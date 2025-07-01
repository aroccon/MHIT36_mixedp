#git clone https://github.com/NVIDIA/cuDecomp
echo "Upload cuDecomp manually, MN5 has no connection"
cd cuDecomp
mkdir build
cd build
module load nano
module load nvidia-hpc-sdk/25.3
module load nccl
cmake ..
make -j
