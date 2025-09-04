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
  
}

#endif
