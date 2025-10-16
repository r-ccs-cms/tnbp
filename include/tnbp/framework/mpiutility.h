/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/mpiutility.h
@brief Utilities for MPI
*/
#ifndef TNBP_FRAMEWORK_MPIUTILITY_H
#define TNBP_FRAMEWORK_MPIUTILITY_H

#include "mpi.h"

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

  template <typename T> MPI_Datatype GetMpiType();

  template<> inline MPI_Datatype GetMpiType<int>()    { return MPI_INT32_T; }
  template<> inline MPI_Datatype GetMpiType<float>()   { return MPI_FLOAT; }
  template<> inline MPI_Datatype GetMpiType<double>()   { return MPI_DOUBLE; }
  template<> inline MPI_Datatype GetMpiType<std::complex<float>>()   { return MPI_CXX_FLOAT_COMPLEX; }
  template<> inline MPI_Datatype GetMpiType<std::complex<double>>()  { return MPI_CXX_DOUBLE_COMPLEX; }
  
  // coor -> address (row-major)
  template <typename ShapeT, typename CoorT>
  inline std::ptrdiff_t address_from_coor(const ShapeT& shape,
					  const CoorT& coor) {
    if (shape.empty()) return 0;
    size_t addr = coor[shape.size()-1];
    for (size_t k = shape.size()-1; k > 0; --k) {
      addr = addr * shape[k-1] + coor[k-1];
    }
    return static_cast<std::ptrdiff_t>(addr);
  }

  template <typename TenT>
  void MpiBcast(context_handle_t<TenT> &ctx,
		TenT &A, int root, MPI_Comm comm) {
    using ShapeT   = typename tci::tensor_traits<TenT>::shape_t;
    using SizeT    = typename tci::tensor_traits<TenT>::ten_size_t;
    using ElemT    = typename tci::tensor_traits<TenT>::elem_t;
    using RankT    = typename tci::tensor_traits<TenT>::rank_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    
    int mpi_rank = -1;
    MPI_Comm_rank(comm, &mpi_rank);
    
    uint32_t size = 0, rank = 0;
    std::vector<uint32_t> shape;
    std::vector<ElemT>    data;
    
    if (mpi_rank == root) {
      const SizeT  size_A  = tci::size(ctx, A);
      const RankT  rank_A  = tci::rank(ctx, A);
      const ShapeT shape_A = tci::shape(ctx, A);
      
      size = static_cast<uint32_t>(size_A);
      rank = static_cast<uint32_t>(rank_A);
      
      shape.resize(rank);
      {
	auto it = shape.begin();
	for (auto d : shape_A) *it++ = static_cast<uint32_t>(d);
      }
      
      data.resize(static_cast<size_t>(size_A));
      auto it_data = data.begin();
      tci::to_range(ctx, A, it_data,
	   [shape_A](const auto& coor) -> std::ptrdiff_t {
	     return address_from_coor(shape_A, coor);
	   });
    }
    
    // bcast meta-data and shape
    MPI_Bcast(&size, 1, MPI_UINT32_T, root, comm);
    MPI_Bcast(&rank, 1, MPI_UINT32_T, root, comm);
    
    if (mpi_rank != root) shape.resize(rank);
    MPI_Bcast(shape.data(), rank, MPI_UINT32_T, root, comm);
    
    // data
    if (mpi_rank != root) data.resize(size);
    MPI_Datatype DataT = GetMpiType<ElemT>();
    MPI_Bcast(data.data(), size, DataT, root, comm);
    
    // Reconstruct A in reciever 
    if (mpi_rank != root) {
      ShapeT shape_A(rank);
      {
	auto it = shape_A.begin();
	for (auto d : shape) *it++ = static_cast<BondDimT>(d);
      }
      auto it_data = data.begin();
      tci::assign_from_range(
	   ctx, shape_A, it_data,
	   [shape_A](const auto& coor) -> std::ptrdiff_t {
	     return address_from_coor(shape_A, coor);
	   },A);
    }
  }


  template <typename TenT>
  void MpiSend(context_handle_t<TenT> &ctx,
	       const TenT &A, int dst, MPI_Comm comm) {
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using SizeT  = typename tci::tensor_traits<TenT>::ten_size_t;
    using ElemT  = typename tci::tensor_traits<TenT>::elem_t;
    using RankT  = typename tci::tensor_traits<TenT>::rank_t;
    
    // shape / rank / size
    const SizeT  size_A  = tci::size(ctx, A);
    const RankT  rank_A  = tci::rank(ctx, A);
    const ShapeT shape_A = tci::shape(ctx, A);
    
    std::vector<uint32_t> shape;
    shape.resize(static_cast<size_t>(rank_A));
    {
      auto it = shape.begin();
      for (auto dim : shape_A) *it++ = static_cast<uint32_t>(dim);
    }
    
    std::vector<ElemT> data;
    data.resize(static_cast<size_t>(size_A));
    {
      auto it_data = data.begin();
      tci::to_range(
	   ctx, A, it_data,
	   [shape_A](const auto& coor) -> std::ptrdiff_t {
	     return address_from_coor(shape_A, coor);
	   });
    }

    const uint32_t size = static_cast<uint32_t>(size_A);
    const uint32_t rank = static_cast<uint32_t>(rank_A);
  
    // Send data
    MPI_Send(&size,  1, MPI_UINT32_T, dst, 0, comm);
    MPI_Send(&rank,  1, MPI_UINT32_T, dst, 1, comm);
    MPI_Send(shape.data(), rank, MPI_UINT32_T, dst, 2, comm);
    
    MPI_Datatype DataT = GetMpiType<ElemT>();
    MPI_Send(data.data(), size, DataT, dst, 3, comm);
  }
  
  template <typename TenT>
  void MpiRecv(context_handle_t<TenT> &ctx,
	       TenT &A, int src, MPI_Comm comm) {
    using ShapeT   = typename tci::tensor_traits<TenT>::shape_t;
    using ElemT    = typename tci::tensor_traits<TenT>::elem_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    
    uint32_t size = 0;
    uint32_t rank = 0;
    std::vector<uint32_t> shape;
    std::vector<ElemT>    data;
    MPI_Status status{};
    
    // Recieve data
    MPI_Recv(&size, 1, MPI_UINT32_T, src, 0, comm, &status);
    MPI_Recv(&rank, 1, MPI_UINT32_T, src, 1, comm, &status);
    
    shape.resize(rank);
    MPI_Recv(shape.data(), rank, MPI_UINT32_T, src, 2, comm, &status);
    
    data.resize(size);
    MPI_Datatype DataT = GetMpiType<ElemT>();
    MPI_Recv(data.data(), size, DataT, src, 3, comm, &status);
    
    // Reconstruction of TenT
    ShapeT shape_A(rank);
    {
      auto it = shape_A.begin();
      for (auto d : shape) *it++ = static_cast<BondDimT>(d);
    }

    auto it_data = data.begin();
    tci::assign_from_range(
       ctx, shape_A, it_data,
       [shape_A](const auto& coor) -> std::ptrdiff_t {
	 return address_from_coor(shape_A, coor);
       },
       A);
  }
  
}

#endif
