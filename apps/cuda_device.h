#pragma once

#ifdef USE_TCAPI_CUDA
#include <cuda_runtime_api.h>
#include <iostream>
#include <mpi.h>

// Launch each process with one visible GPU (e.g. Slurm --gpus-per-task=1).
// Do this before creating any TCAPI context, including secondary real contexts.
inline void initialize_tensor_device(MPI_Comm comm) {
    int rank=0;
    MPI_Comm_rank(comm,&rank);
    int count=0;
    auto status=cudaGetDeviceCount(&count);
    if(status!=cudaSuccess || count!=1) {
        std::cerr << "rank " << rank
                  << ": tcapi-cuda requires one visible GPU per process; "
                  << "use launcher GPU binding (e.g. --gpus-per-task=1). "
                  << "Visible GPUs=" << count << ", CUDA=" << cudaGetErrorString(status) << '\n';
        MPI_Abort(comm,1);
        return;
    }
    status=cudaSetDevice(0);
    if(status!=cudaSuccess) {
        std::cerr << "rank " << rank << ": cudaSetDevice: " << cudaGetErrorString(status) << '\n';
        MPI_Abort(comm,1);
    }
}
#else
#include <mpi.h>
inline void initialize_tensor_device(MPI_Comm) {}
#endif
