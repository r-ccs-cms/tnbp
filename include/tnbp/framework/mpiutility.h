/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/mpiutility.h
@brief Utilities for MPI
*/
#ifndef TNBP_FRAMEWORK_MPIUTILITY_H
#define TNBP_FRAMEWORK_MPIUTILITY_H

namespace tnbp {

  void get_range(int size, int rank, int & i_start, int & i_end)
  {
    int i_all = i_end - i_start;

    int i_div = i_all / size;
    int i_res = i_all % size;

    int j = i_start;
    int j_start;
    int j_end;
    if( rank < i_res )
      {
        j_start = i_start + ( i_div + 1 ) * rank;
        j_end = i_start + ( i_div + 1 ) * ( rank + 1 );
      }
    else
      {
        j_start = i_start + ( i_div + 1 ) * i_res + i_div * ( rank - i_res );
        j_end = i_start + ( i_div + 1 ) * i_res + i_div * ( rank + 1 - i_res );
      }
    i_start = j_start;
    i_end = j_end;
  }

  void MpiBcast(std::string &str, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // share the size
    int size = 0;
    if (rank == root) {
      size = static_cast<int>(str.size());
    }
    MPI_Bcast(&size, 1, MPI_INT, root, comm);
    
    // reserve buffer for non-root processs
    if (rank != root) {
      str.resize(size);
    }
    
    // share the data
    if (size > 0) {
      MPI_Bcast(str.data(), size, MPI_CHAR, root, comm);
    }
  }
  
}

#endif
