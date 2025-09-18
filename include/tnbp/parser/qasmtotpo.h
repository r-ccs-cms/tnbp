/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/parser/qasmtotpo.h
@brief Function to transform qasm to tpo
*/
#ifndef TNBP_PARSER_QASMTOTPO_H
#define TNBP_PARSER_QASMTOTPO_H

#include <cmath>
#include "qasm/ir.hpp"

namespace tnbp {

  /**
     Function to construct matrix from qasm::Instruction
     @param[in] ins: qasm::Instruction corresponding to the gate
   */
  template<typename ElemT>
  std::vector<ElemT> InstructionMatrix(const qasm::Instruction & ins) {
    const double pi = 3.14159265358979323846;

    switch (ins.op) {
      // ===== 1 qubit family =====
    case qasm::Op::U3: {
      double th=ins.params[0], ph=ins.params[1], la=ins.params[2];
      return {
	ElemT(std::cos(th/2),0), -std::exp(ElemT(0,la))*std::sin(th/2),
	std::exp(ElemT(0,ph))*std::sin(th/2), std::exp(ElemT(0,ph+la))*std::cos(th/2)
      };
    }
    case qasm::Op::U2: {
      double ph=ins.params[0], la=ins.params[1];
      double th=pi/2;
      double s=1/std::sqrt(2.0);
      return {
	ElemT(s,0), -std::exp(ElemT(0,la))*s,
	std::exp(ElemT(0,ph))*s, std::exp(ElemT(0,ph+la))*s
      };
    }
    case qasm::Op::U1: {
      double la=ins.params[0];
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), std::exp(ElemT(0,la)) };
    }
    case qasm::Op::RX: {
      double th=ins.params[0];
      return {
	ElemT(std::cos(th/2),0), ElemT(0,-std::sin(th/2)),
	ElemT(0,-std::sin(th/2)), ElemT(std::cos(th/2),0)
      };
    }
    case qasm::Op::RY: {
      double th=ins.params[0];
      return {
	ElemT(std::cos(th/2),0), ElemT(-std::sin(th/2),0),
	ElemT(std::sin(th/2),0), ElemT(std::cos(th/2),0)
      };
    }
    case qasm::Op::RZ: {
      double th=ins.params[0];
      return {
	std::exp(ElemT(0,-th/2)), ElemT(0,0),
	ElemT(0,0), std::exp(ElemT(0,th/2))
      };
    }
    case qasm::Op::H: {
      double s=1/std::sqrt(2.0);
      return { ElemT(s,0), ElemT(s,0), ElemT(s,0), ElemT(-s,0) };
    }
    case qasm::Op::X:
      return { ElemT(0,0), ElemT(1,0), ElemT(1,0), ElemT(0,0) };
    case qasm::Op::Y:
      return { ElemT(0,0), ElemT(0,-1), ElemT(0,1), ElemT(0,0) };
    case qasm::Op::Z:
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), ElemT(-1,0) };
    case qasm::Op::S:
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), ElemT(0,1) };
    case qasm::Op::SDG:
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), ElemT(0,-1) };
    case qasm::Op::T:
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), std::exp(ElemT(0,pi/4)) };
    case qasm::Op::TDG:
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), std::exp(ElemT(0,-pi/4)) };
    case qasm::Op::ID:
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), ElemT(1,0) };
      
      // ===== 2 qubit gates (4x4) =====
    case qasm::Op::CX: {
      // |00>->|00>, |01>->|01>, |10>->|11>, |11>->|10>
      std::vector<ElemT> m(16,ElemT(0,0));
      m[0*4+0]=1; m[1*4+1]=1; m[2*4+3]=1; m[3*4+2]=1;
      return m;
    }
    case qasm::Op::CZ: {
      // |11>に -1 の位相
      std::vector<ElemT> m(16,ElemT(0,0));
      m[0*4+0]=1; m[1*4+1]=1; m[2*4+2]=1; m[3*4+3]=-1;
      return m;
    }
    case qasm::Op::SWAP: {
      // |01> <-> |10>
      std::vector<ElemT> m(16,ElemT(0,0));
      m[0*4+0]=1; m[1*4+2]=1; m[2*4+1]=1; m[3*4+3]=1;
      return m;
    }
      
      // ===== 3 qubit gates (8x8) =====
    case qasm::Op::CCX: {
      // Toffoli: |110>->|111>, |111>->|110>
      std::vector<ElemT> m(64,ElemT(0,0));
      for(int i=0;i<6;i++) m[i*8+i]=1; // |000>..|101> 固定
      m[6*8+7]=1; // |110>->|111>
      m[7*8+6]=1; // |111>->|110>
      return m;
    }
    case qasm::Op::CSWAP: {
      // Fredkin: 制御=1のとき |101><->|110>
      std::vector<ElemT> m(64,ElemT(0,0));
      // |000>..|100>,|111> はそのまま
      m[0*8+0]=1; m[1*8+1]=1; m[2*8+2]=1; m[3*8+3]=1;
      m[4*8+4]=1; m[7*8+7]=1;
      // swap |101> <-> |110>
      m[5*8+6]=1;
      m[6*8+5]=1;
      return m;
    }
      
    default:
      throw std::runtime_error("matrix for this op not implemented");
    }
  }

  /**
     Function to get the size of qubits of the gate defined by qasm::Instruction
     @param[in] ins: qasm::Instruction corresponding to the gate
   */
  inline int OpQubitCount(const qasm::Instruction & ins) {
    switch (ins.op) {
    case Op::U3: case Op::U2: case Op::U1:
    case Op::RX: case Op::RY: case Op::RZ:
    case Op::H:  case Op::X:  case Op::Y:  case Op::Z:
    case Op::S:  case Op::SDG: case Op::T: case Op::TDG:
    case Op::ID:
      return 1;
    case Op::CX: case Op::CZ: case Op::SWAP:
      return 2;
    case Op::CCX: case Op::CSWAP:
      return 3;
    case Op::MEASURE:
      return 1;
    case Op::RESET:
      return 1;
    case Op::BARRIER:
    case Op::CUSTOM:
      return ins.qubits.size(); // 可変長は実データ依存
    }
    return ins.qubits.size();
  }  
  
  /**
     Function to get the tensor form of the gate defined by qasm::Instruction
     @param[in] ins: qasm::Instruction corresponding to the gate
   */
  template <typename TenT>
  TenT InstructionTensor(context_handle_t<TenT> & ctx,
			 const qasm::Instruction & ins) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    std::vector<BondDimT> shape(2*OpQubitCount(ins),2);
    std::vector<ElemT> data = InstructionMatrix(ins);
    auto it_data = data.begin();
    TenT res;
    tci::asign_from_container(
	 ctx,shape,it_data,
	 [&shape](const CoorsT & coor) {
	   return address_from_coor(shape,coor);
	 },res);
    return res;
  }

  /**
     Function to construct tensor product operator from qasm::Program
     @param[in] ctx: context_handle for tensor operations
     @param[in] program: circuit composed of qasm2/qasm3 instruction
     @param[in] edges: edges to define the tensor network
     @param[in] num_gates: number of gates stacked on the one-slice of the tensor-product-operator
  */
  template <typename TenT>
  std::vector<std::vector<TenT>> QasmToTPO(
		  context_handle_t<TenT> & ctx,
		  const qasm::Program & program,
		  const std::vector<std::pair<int,int>> & edges,
		  const std::vector<int> & num_gates) {
    
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tci::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tci::create_context(ctx_r);
    
    auto site = GetSiteIndexFromBond(edges);
    std::vector<std::vector<TenT>> T(num_gates,
		std::vector<TenT>(site.size()));
    int gate_count = 0;
    int m = 0;
    for(auto const & ins : program.instructions) {

      // initialization of tensor product operator
      if( gate_count == 0 ) {
	for(auto const & i : site) {
	  auto virtualbond = GetSurroundingBondIndex(i,edges);
	  auto num_virtualbond = virtualbond.size();
	  ShapeT bond(num_virtualbond+2,1);
	  bond[num_virtualbond+0] = 2;
	  bond[num_virtualbond+1] = 2;
	  std::vector<ElemT> data(4,ElemT(0.0));
	  data[0] = ElemT(1.0);
	  data[3] = ElemT(1.0);
	  tci::assign_from_container(ctx,bond,it_data,
	      [&bond](const CoorsT & c) {
		return address_from_coor(bond,c);
	      }, T[m][i]);
	}
      }
      
      TenT gate = InstructionTensor(ctx,ins);
      auto rank_tensor = tci::rank(ctx,gate);
      int num_qubits = rank_tensor/2;
      if( num_qubits == 1 ) {
	auto site_a = static_cast<int>(ins.qubits[0].index);
	auto virtualbond = GetSurroundingBondIndex(site_a,edges);
	auto num_virtualbond = virtualbond.size();
	List<BondLabelT> IdxT(num_virtualbond+2);
	List<BondLabelT> IdxG(2);
	List<BondLabelT> IdxR(num_virtualbond+2);
	std::iota(IdxT.begin(),IdxT.end(),0);
	std::iota(IdxR.begin(),IdxR.end(),0);
	IdxT[num_virtualbond+0] = static_cast<BondLabelT>(-1);
	IdxT[num_virtualbond+1] = static_cast<BondLabelT>(num_virtualbond)+1;
	IdxG[0] = static_cast<BondLabelT>(num_virtualbond);
	IdxG[1] = static_cast<BondLabelT>(-1);
	tci::contract(ctx,T[m][site_a],IdxT,gate,IdxG,T[m][site_a],IdxR);
      } else if ( num_qubits == 2 ) {

	// we first perform svd
	List<BondIdxT> new_order(4);
	new_order[0] = static_cast<BondLabelT>(0);
	new_order[1] = static_cast<BondLabelT>(2);
	new_order[2] = static_cast<BondLabelT>(1);
	new_order[3] = static_cast<BondLabelT>(3);
	tci::transpose(ctx,gate,new_order);
	TenT A;
	TenT B;
	RealTenT S;
	RankT num_row_bonds = 2;
	tci::svd(ctx,gate,num_row_bonds,A,S,B);
	tci::for_each(ctx_r,S,[](ElemT & elem) {
	  if( elem > 0.0 ) { elem = std::sqrt(elem); }
	  else { elem = 0.0; }});
	TenT D;
	tci::convert(ctx_r,S,ctx,D);
	tci::diag(ctx,D);
	List<BondLabelT> IdxG(3);
	List<BondLabelT> IdxD(2);
	List<BondLabelT> IdxGnew(3);
	IdxG[0] = static_cast<BondLabelT>(0);
	IdxG[1] = static_cast<BondLabelT>(1);
	IdxG[2] = static_cast<BondLabelT>(-1);
	IdxD[0] = static_cast<BondLabelT>(-1);
	IdxD[1] = static_cast<BondLabelT>(2);
	std::iota(IdxGnew.begin(),IdxGnew.end(),0);
	tci::contract(ctx,A,IdxG,D,IdxD,A,IdxGnew);
	IdxD[0] = static_cast<BondLabelT>(0);
	IdxD[1] = static_cast<BondLabelT>(-1);
	IdxG[0] = static_cast<BondLabelT>(-1);
	IdxG[1] = static_cast<BondLabelT>(1);
	IdxG[2] = static_cast<BondLabelT>(2);
	tci::contract(ctx,D,IdxD,B,IdxG,B,IdxGnew);
	ShapeT shape_D = tci::shape(ctx,D);
	auto vdim = shape_D[0];

	// Assign  A and B on the tensor-product-operators
	auto site_a = static_cast<int>(ins.qubits[0].index);
	auto site_b = static_cast<int>(ins.qubits[1].index);
	auto path = FindShortestPath(edges,site_a,site_b);

	for(size_t i=0; i < path.size(); i++) {

	  TenT X;
	  
	  if( 0 < i && i < path.size()-1 ) {
	    
	    ShapeT dimX(1,vdim*2);
	    std::vector<ElemT> dataX(2*vdim,ElemT(1.0));
	    auto it_dataX = dataX.begin();
	    tci::assign_from_container(ctx,dimX,it_dataX,
		      [](const Coors & c) {
			return c[0]; },X);
	    tci::diag(ctx,X);
	    ShapeT shapeX(4);
	    shapeX[0] = static_cast<BondDimT>(vdim);
	    shapeX[1] = static_cast<BondDimT>(2);
	    shapeX[2] = static_cast<BondDimT>(vdim);
	    shapeX[3] = static_cast<BondDimT>(2);
	    tci::reshape(ctx,X,shapeX);
	    List<BondIdxT> new_order_X(4);
	    new_order_X[0] = static_cast<BondIdxT>(0);
	    new_order_X[1] = static_cast<BondIdxT>(2);
	    new_order_X[2] = static_cast<BondIdxT>(1);
	    new_order_X[3] = static_cast<BondIdxT>(3);
	    tci::transpose(ctx,X,new_order_X);
	  } else if ( i == 0 ) {
	    tci::copy(ctx,A,X);
	  } else if ( i == path.size()-1 ) {
	    tci::copy(ctx,B,X);
	  }
	  auto vb_a = GetSurroundingBondIndex(path[i],edges);
	  auto rank_X = tci::rank(ctx,X);
	  List<BondLabelT> IdxA(vb_a.size()+2,1);
	  List<BondLabelT> IdxX(rank_X);
	  int target_bond_m = static_cast<int>(vb_a.size())+2;
	  if( i > 0 ) {
	    target_bond_m = 0;
	    for(auto const & k : vb_a) {
	      if( edges[k].first == path[i]
		  && edges[k].second == path[i-1] ) {
		break;
	      }
	      if( edges[k].second == path[i]
		  && edges[k].first == path[i-1] ) {
		break;
	      }
	      target_bond_m++;
	    }
	  }
	  int target_bond_p = vb_a.size()+2;
	  if( i < path.size()-1 ) {
	    target_bond_p = 0;
	    for(auto const & k : vb_a) {
	      if( edges[k].first == path[i]
		  && edges[k].second == path[i+1] ) {
		break;
	      }
	      if( edges[k].second == path[i]
		  && edges[k].first == path[i+1] ) {
		break;
	      }
	      target_bond_p++;
	    }
	  }
	  int k = 0;
	  int kx = 0;
	  ShapeT shapeA = tci::shape(ctx,T[m][path[i]]);
	  ShapeT new_shapeA(vb_a.size()+2);
	  for(int b=0; b < vb_a.size(); b++) {
	    if( b == target_bond_m ) {
	      IdxA[b] = static_cast<BondLabelT>(k);
	      k++;
	      IdxX[kx] = static_cast<BondLabelT>(k);
	      k++;
	      kx++;
	      new_shapeA[b] = static_cast<BondDimT>(vdim*shapeA[b]);
	    } else if ( b == target_bond_p ) {
	      IdxA[b] = static_cast<BondLabelT>(k);
	      k++;
	      IdxX[kx] = static_cast<BondLabelT>(k);
	      k++;
	      kx++;
	      new_shapeA[b] = static_cast<BondDimT>(vdim*shapeA[b]);
	    } else {
	      IdxA[b] = static_cast<BondLabelT>(k);
	      k++;
	      new_shapeA[b] = shapeA[b];
	    }
	  }
	  IdxX[kx] = static_cast<BondLabelT>(vb_a.size());
	  IdxX[kx+1] = static_cast<BondLabelT>(-1);
	  IdxA[vb_a.size()+0] = static_cast<BondLabelT>(-1);
	  IdxA[vb_a.size()+1] = static_cast<BondLabelT>(vb_a.size()+1);
	  new_shapeA[vb_a.size()+0] = shapeA[vb_a.size()+0];
	  new_shapeA[vb_a.size()+1] = shapeA[vb_a.size()+1];
	  List<BondLabelT> IdxP(vb_a.size()+2);
	  std::iota(IdxP.begin(),IdxP.end(),0);
	  tci::contract(ctx,X,IdxX,T[m][path[i]],IdxA,T[m][path[i]],IdxP);
	  tci::reshape(ctx,T[m][path[i]],new_shapeA);
	}
	
      } else if ( num_qubits == 3 ) {

	// choose shortest path
	auto site_a = static_cast<int>(ins.qubits[0].index);
	auto site_b = static_cast<int>(ins.qubits[1].index);
	auto site_c = static_cast<int>(ins.qubits[2].index);

	auto path_ab = FindShortestPath(edges,site_a,site_b);
	auto path_bc = FindShortestPath(edges,site_b,site_c);
	auto path_ca = FindShortestPath(edges,site_c,site_a);

	std::vector<size_t> path_length(3);
	path_length[0] = path_ab.size()+path_bc.size()-1;
	path_length[1] = path_bc.size()+path_ca.size()-1;
	path_length[2] = path_ca.size()+path_ab.size()-1;

	size_t min_length = path_length[0];
	int which_path = 0;
	
	if( path_length[1] < path_length[0] ) {
	  which_path = 1;
	  min_length = path_length[1];
	}
	if( path_length[2] < path_length[which_path] ) {
	  which_path = 2;
	  min_length = path_length[2];
	}

	std::vector<int> path(min_length);
	std::vector<BondIdxT> new_order(6);
	if( which_path == 0 ) {
	  for(int i=0; i < path_length[0]-1; i++) {
	    path[i] = path_ab[i];
	  }
	  for(int i=0; i < path_length[1]; i++) {
	    path[i+path_length[0]] = path_bc[i];
	  }
	  new_order[0] = 0;
	  new_order[1] = 3;
	  new_order[2] = 1;
	  new_order[3] = 4;
	  new_order[4] = 2;
	  new_order[5] = 5;
	} else if ( which_path == 1 ) {
	  for(int i=0; i < path_length[1]-1; i++) {
	    path[i] = path_bc[i];
	  }
	  for(int i=0; i < path_length[2]; i++) {
	    path[i+path_length[1]] = path_ca[i];
	  }
	  new_order[0] = 1;
	  new_order[1] = 4;
	  new_order[2] = 2;
	  new_order[3] = 5;
	  new_order[4] = 0;
	  new_order[5] = 3;
	  site_a = static_cast<int>(ins.qubits[1].index);
	  site_b = static_cast<int>(ins.qubits[2].index);
	  site_c = static_cast<int>(ins.qubits[0].index);
	} else {
	  for(int i=0; i < path_length[2]-1; i++) {
	    path[i] = path_ca[i];
	  }
	  for(int i=0; i < path_length[0]; i++) {
	    path[i+path_length[2]] = path_ab[i];
	  }
	  new_order[0] = 2;
	  new_order[1] = 5;
	  new_order[2] = 0;
	  new_order[3] = 3;
	  new_order[4] = 1;
	  new_order[5] = 4;
	  site_a = static_cast<int>(ins.qubits[2].index);
	  site_b = static_cast<int>(ins.qubits[0].index);
	  site_c = static_cast<int>(ins.qubits[1].index);
	}
	tci::transpose(ctx,gate,new_order);

	// decompose a 3-qubit unitary into network of three 1-qubit tensors
	TenT Ga;
	TenT Gb;
	TenT Gc;
	RealTenT S;
	TenT Gv; // temporary
	RankT num_row_bonds = 2;
	tci::svd(ctx,gate,num_row_bonds,Ga,S,Gv);
	tci::for_each(ctx_r,S,[](ElemT & elem) {
	  if( elem > 0.0 ) { elem = std::sqrt(elem); }
	  else             { elem = 0.0; }
	});
	TenT D;
	tci::convert(ctx_r,S,ctx,D);
	tci::diag(ctx,D);
	List<BondLabelT> IdxGa(3);
	List<BondLabelT> IdxD(2);
	List<BondLabelT> IdxGa_new(3);
	IdxGa[0] = static_cast<BondLabelT>(0);
	IdxGa[1] = static_cast<BondLabelT>(1);
	IdxGa[2] = static_cast<BondLabelT>(-1);
	IdxD[0] = static_cast<BondLabelT>(-1);
	IdxD[1] = static_cast<BondLabelT>(2);
	std::iota(IdxGa_new.begin(),IdxGa_new.end(),0);
	tci::contract(ctx,Ga,IdxGa,D,IdxD,Ga,IdxGa_new);
	List<BondLabelT> IdxGv(5);
	List<BondLabelT> IdxGv_new(5);
	IdxD[0] = static_cast<BondLabelT>(0);
	IdxD[1] = static_cast<BondLabelT>(-1);
	IdxGv[0] = static_cast<BondLabelT>(-1);
	IdxGv[1] = static_cast<BondLabelT>(1);
	IdxGv[2] = static_cast<BondLabelT>(2);
	IdxGv[3] = static_cast<BondLabelT>(3);
	IdxGv[4] = static_cast<BondLabelT>(4);
	std::iota(IdxGv_new.begin(),IdxGv_new.end(),0);
	tci::contract(ctx,D,IdxD,Gv,IdxGv,Gv,IdxGv_new);
	ShapeT shape_Dab = tci::shape(ctx,D);
	num_row_bonds = 3;
	tci::svd(ctx,Gv,num_row_bonds,Gb,S,Gc);
	tci::for_each(ctx,S,[](ElemT & elem) {
	  if( elem > 0.0 ) { elem = std::sqrt(elem); }
	  else             { elem = 0.0; }
	});
	tci::convert(ctx_r,S,ctx,D);
	tci::diag(ctx,D);
	List<BondLabelT> IdxGb(4);
	List<BondLabelT> IdxGb_new(4);
	IdxGb[0] = static_cast<BondLabelT>(0);
	IdxGb[1] = static_cast<BondLabelT>(2);
	IdxGb[2] = static_cast<BondLabelT>(3);
	IdxGb[3] = static_cast<BondLabelT>(-1);
	IdxD[0] = static_cast<BondLabelT>(-1);
	IdxD[1] = static_cast<BondLabelT>(1);
	std::iota(IdxGb_new.begin(),IdxGb_new.end(),0);
	tci::contract(ctx,Gb,IdxGb,D,IdxD,Gb,IdxGb_new);
	List<BondLabelT> IdxGc(3);
	List<BondLabelT> IdxGc_new(3);
	IdxD[0] = static_cast<BondLabelT>(0);
	IdxD[1] = static_cast<BondLabelT>(-1);
	IdxGc[0] = static_cast<BondLabelT>(-1);
	IdxGc[1] = static_cast<BondLabelT>(1);
	IdxGc[2] = static_cast<BondLabelT>(2);
	std::iota(IdxGc_new.begin(),IdxGc_new.end(),0);
	tci::contract(ctx,D,IdxD,Gc,IdxGc,Gc,IdxGc_new);
	ShapeT shape_Dbc = tci::shape(ctx,D);

	auto vdim_ab = shape_Dab[0];
	auto vdim_bc = shape_Dbc[0];

	auto vdim = vdim_ab;
	size_t ib=path.size();
	// attach unitary onto the tpo
	for(size_t i=0; i < path.size(); i++) {

	  TenT X;
	  if( path[i] == site_a ) {
	    tci::copy(ctx,Ga,X);
	  } else if ( path[i] == site_b ) {
	    tci::copy(ctx,Gb,X);
	    ib = i;
	  } else if ( path[i] == site_c ) {
	    tci::copy(ctx,Gc,X);
	  } else {
	    if( i > ib ) {
	      vdim = vdim_bc;
	    }
	    ShapeT dimX(1,vdim*2);
	    std::vector<ElemT> dataX(2*vdim,ElemT(1.0));
	    auto it_dataX = dataX.begin();
	    tci::assign_from_container(ctx,dimX,it_dataX,
		  [](const CoorsT & c) {
		    return c[0];
		  }, X);
	    tci::diag(ctx,X);
	    ShapeT shapeX(4);
	    shapeX[0] = static_cast<BondDimT>(vdim);
	    shapeX[1] = static_cast<BondDimT>(2);
	    shapeX[2] = static_cast<BondDimT>(vdim);
	    shapeX[3] = static_cast<BondDimT>(2);
	    tci::reshape(ctx,X,shapeX);
	    List<BondIdxT> new_order_X(4);
	    new_order_X[0] = static_cast<BondIdxT>(0);
	    new_order_X[1] = static_cast<BondIdxT>(2);
	    new_order_X[2] = static_cast<BondIdxT>(1);
	    new_order_X[3] = static_cast<BondIdxT>(3);
	    tci::transpose(ctx,X,new_order_X);
	  }
	  auto vb_A = GetSurroundingBondIndex(path[i],edges);
	  auto rank_X = tci::rank(ctx,X);
	  List<BondLabelT> IdxA(vb_A.size()+2,1);
	  List<BondLabelT> IdxX(rank_X);
	  List<BondLabelT> IdxN(vb_A.size()+rank_X,1);
	  int target_bond_m = static_cast<int>(vb_A.size())+2;
	  if( i > 0 ) {
	    target_bond_m = 0;
	    for(auto const & k : vb_A) {
	      if( edges[k].first == path[i]
		  && edges[k].second == path[i-1] ) {
		break;
	      }
	      if( edges[k].second == path[i]
		  && edges[k].first == path[i-1] ) {
		break;
	      }
	      target_bond_m++;
	    }
	  }
	  int target_bond_p = vb_A.size()+2;
	  if( i < path.size()-1 ) {
	    target_bond_p = 0;
	    for(austo const & k : vb_a) {
	      if( edges[k].first == path[i]
		  && edges[k].second == path[i+1] ) {
		break;
	      }
	      if( edges[k].second == path[i]
		  && edges[k].first == path[i+1] ) {
		break;
	      }
	      target_bond_p++;
	    }
	  }
	  int k=0;
	  int kx = 0;
	  ShapeT shapeA = tci::shape(ctx,T[m][path[i]]);
	  ShapeT new_shapeA(vb_A.size()+2);
	  for(int b=0; b < vb_A.size(); b++) {
	    if( b == target_bond_m ) {
	      IdxA[b] = static_cast<BondLabelT>(k);
	      k++;
	      IdxX[kx] = static_cast<BondLabelT>(k);
	      k++;
	      kx++;
	      new_shapeA[b] = static_cast<BondDimT>(vdim*shapeA[b]);
	    } else if ( b == target_bond_p ) {
	      IdxA[b] = static_cast<BondLabelT>(k);
	      k++;
	      IdxX[kx] = static_cast<BondLabelT>(k);
	      k++;
	      kx++;
	      new_shapeA[b] = static_cast<BondDimT>(vdim*shapeA[b]);
	    } else {
	      IdxA[b] = static_cast<BondLabelT>(k);
	      k++;
	      new_shapeA[b] = shapeA[b];
	    }
	  }
	  IdxX[kx] = static_cast<BondLabelT>(vb_A.size());
	  IdxX[kx+1] = static_cast<BondLabelT>(-1);
	  IdxA[vb_A.size()+0] = static_cast<BondLabelT>(-1);
	  IdxA[vb_A.size()+1] = static_cast<BondLabelT>(vb_A.size()+1);
	  std::iota(IdxN.begin(),IdxN.end(),0);
	  tci::contract(ctx,X,IdxX,T[m][path[i]],IdxA,T[m][path[i]],IdxN);
	  new_shapeA[vb_A.size()+0] = shapeA[vb_A.size()+0];
	  new_shapeA[vb_A.size()+1] = shapeA[vb_A.size()+1];
	  tci::reshape(ctx,T[m][path[i]],new_shapeA);
	}
      } // else if ( num_qubits == 3 )

      gate_count++;
      if( gate_count == num_gates[m] ) {
	m++;
	gate_count=0;
      }
    } // end for(auto const & ins : p.instructions)
    return T;
  } // end function

}

#endif
