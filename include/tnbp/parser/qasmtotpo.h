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
  
  template <typename TenT>
  TenT InstructionTensor(context_handle_t<TenT> & ctx,
			 const qasm::Instruction & ins) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    std::vector<int> bond(2*OpQubitCount(ins),2);
    std::vector<ElemT> data = InstructionMatrix(ins);
    return tci::initialize<TenT>(ctx,bond,data);
  }

  /**
     Function to construct tensor product operator
  */
  template <typename TenT>
  std::vector<std::vector<TenT>> QasmToTPO(
		  context_handle_t<TenT> & ctx,
		  const qasm::Program & program,
		  const std::vector<std::pair<int,int>> & edges,
		  std::vector<int> & num_gates) {

    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using IntT = typename tci::tensor_traits<TenT>::bond_label_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    
    auto site = GetSiteIndexFromBond(edges);
    std::vector<std::vector<TenT>> res(num_gates,std::vector<TenT>(site.size()));
    int gate_count = 0;
    int m = 0;
    for(auto const & ins : p.instructions) {

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
	  T[m][i] = tci::initialize<TenT>(ctx,bond,data);
	}
      }
      
      TenT gate = InstructionTensor(ctx,ins);
      auto rank_tensor = tci::rank(ctx,gate);
      int num_qubits = rank_tensor/2;
      if( num_qubits == 1 ) {
	auto site_a = static_cast<int>(ins.qubits[0].index);
	auto virtualbond = GetSurroundingBondIndex(site_a,edges);
	auto num_virtualbond = virtualbond.size();
	std::vector<IntT> IdxT(num_virtualbond+2);
	std::vector<IntT> IdxG(2);
	std::iota(IdxT.begin(),IdxT.end(),0);
	IdxT[num_virtualbond+0] = -1;
	IdxT[num_virtualbond+1] = static_cast<IntT>(num_virtualbond)+1;
	IdxG[0] = static_cast<IntT>(num_virtualbond);
	IdxG[1] = -1;
	tci::contract(ctx,T[m][site_a],IdxT,gate,IdxG,T[m][site_a]);
      } else if ( num_qubits == 2 ) {

	// we first perform svd
	std::vector<BondIdxT> new_order(4);
	new_order[0] = 0;
	new_order[1] = 2;
	new_order[2] = 1;
	new_order[3] = 3;
	tci::transpose(ctx,gate,new_order);
	TenT A;
	TenT B;
	RealTenT D;
	RankT num_row_bonds = 2;
	tci::svd(ctx,gate,num_row_bonds,A,D,B);
	tci::for_each(ctx,D,[](ElemT & elem) {
	  if( elem > 0.0 ) { elem = std::sqrt(elem); }
	  else { elem = 0.0; }});
	tci::diag(ctx,D);
	std::vector<IntT> IdxG(3);
	std::vector<IntT> IdxD(2);
	IdxG[0] = 0;
	IdxG[1] = 1;
	IdxG[2] = -1;
	IdxD[0] = -1;
	IdxD[1] = 2;
	tci::contract(ctx,A,IdxG,D,IdxD,A);
	IdxD[0] = 0;
	IdxD[1] = -1;
	IdxG[0] = -1;
	IdxG[1] = 1;
	IdxG[2] = 2;
	tci::contract(ctx,D,IdxD,B,IdxG,B);
	ShapeT shape_D = tci::shape(ctx,D);
	auto vdim = shape_D[0];

	// Assign  A and B on the tensor-product-operators
	auto site_a = static_cast<int>(ins.qubits[0].index);
	auto site_b = static_cast<int>(ins.qubits[1].index);
	auto path = FindShortestPath(edges,site_a,site_b);

	for(size_t i=0; i < path.size(); i++) {

	  if( 0 < i && i < path.size()-1 ) {
	    ShapeT dimX(1,vdim*2);
	    std::vector<RealT> dataX(2*vdim,RealT(1.0));
	    RealTenT X = tci::initialize<RealTenT>(ctx,dimX,dataX);
	    tci::diag(ctx,X);
	    ShapeT shapeX(4);
	    shapeX[0] = vdim;
	    shapeX[1] = 2;
	    shapeX[2] = vdim;
	    shapeX[3] = 2;
	    tci::reshape(ctx,X,shapeX);
	    std::vector<BondIdxT> new_order_X(4);
	    new_order_X[0] = 0;
	    new_order_X[1] = 2;
	    new_order_X[2] = 1;
	    new_order_X[3] = 3;
	    tci::transpose(ctx,X,new_order_X);
	    std::vector<IntT> IdxX(4);
	    
	    auto vb_a = GetSurroundingBondIndex(path[i],edges);
	    std::vector<IntT> IdxA(vb_a.size()+2,1);
	    int target_bond_m = 0;
	    for(auto const & m : vb_a) {
	      if( edges[m].first == path[i]
		  && edges[m].second == path[i-1] ) {
		break;
	      }
	      if( edges[m].second == path[i]
		  && edges[m].first == path[i-1] ) {
		break;
	      }
	      target_bond_m++;
	    }
	    int target_bond_p = 0;
	    for(auto const & m : vb_a) {
	      if( edges[m].first == path[i]
		  && edges[m].second == path[i+1] ) {
		break;
	      }
	      if( edges[m].second == path[i]
		  && edges[m].first == path[i+1] ) {
		break;
	      }
	      target_bond_p++;
	    }
	    int k=0;
	    ShapeT shapeA = tci::shape(ctx,T[m][path[i]]);
	    ShapeT new_shapeA(vb_a.size()+2);
	    for(int b=0; b < vb_a.size(); b++) {
	      if( b == target_bond_m ) {
		IdxA[b] = k;
		k++;
		IdxX[0] = k;
		k++;
		new_shapeA[b] = vdim*shapeA[b];
	      } else if ( b == target_bond_p ) {
		IdxA[b] = k;
		k++;
		IdxX[1] = k;
		k++;
		new_shapeA[b] = vdim*shapeA[b];
	      } else {
		IdxA[b] = k;
		k++;
		new_shapeA[b] = shapeA[b];
	      }
	    }
	    IdxX[2] = static_cast<IntT>(vb_a.size());
	    IdxX[3] = -1;
	    IdxA[vb_a.size()+0] = -1;
	    IdxA[vb_a.size()+1] = static_cast<IntT>(vb_a.size()+1);
	    new_shapeA[vb_a.size()+0] = shapeA[vb_a.size()+0];
	    new_shapeA[vb_a.size()+1] = shapeA[vb_a.size()+1];
	    tci::contract(ctx,T[m][path[i]],IdxA,X,IdxX,T[m][path[i]]);
	    tci::reshape(ctx,T[m][path[i]],new_shapeA);
	  } else if ( i == 0 ) {
	    // now writing
	  } else if ( i == path.size()-1 ) {
	    // now writing
	  }
	}

	

	
      } else if ( num_qubits == 3 ) {
	
      }

      gate_count++;
      if( gate_count == num_gates[m] ) {
	m++;
	gate_count=0;
      }
    }
  }

  /**
   */
  
}

#endif
