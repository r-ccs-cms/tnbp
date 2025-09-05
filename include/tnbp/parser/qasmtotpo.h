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
	  std::vector<int> bond(num_virtualbond+2,1);
	  bond[num_virtualbond+0] = 2;
	  bond[num_virtualbond+1] = 2;
	  std::vector<ElemT> data(4,ElemT(0.0));
	  data[0] = ElemT(1.0);
	  data[3] = ElemT(1.0);
	  T[m][i] = tci::initialize<TenT>(ctx,bond,data);
	}
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
