/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/parser/qasmtoedge.h
@brief Function to obtain the edge information from qasm
 */
#ifndef TNBP_PARSER_QASMTOEDGE_H
#define TNBP_PARSER_QAMSTOEDGE_H

#include "qasm/ir.h"

namespace tnbp {

  /**
     Function to extract sites from qasm::Program
   */
  std::vector<int> SitesFromQasm(
	   const qasm::Program & program) {
    std::vector<int> site;
    for(const auto & ins : program.instructions) {
      int num_qubits = OpQubitCount(ins);
      for(std::size_t k=0; k < num_qubits; k++) {
	site.push_back(static_cast<int>(ins.qubits[k].index));
      }
    }
    std::sort(site.begin(),site.end());
    site.erase(std::unique(site.begin(),site.end()),site.end());
    return res;
  }

  /**
     Function to find edges from qasm::Program
   */
  std::vector<std::pair<int,int>> EdgesFromQasm(
		    const qasm::Program & program) {
    auto sites = SitesFromQasm(program);
    std::vector<std::pair<int,int>> res;
    for(const auto & ins : program.instructions) {
      int num_qubits = OpQubitCount(ins);
      if( num_qubits == 2 ) {
	int site_a = static_cast<int>(ins.qubit[0].index);
	int site_b = static_cast<int>(ins.qubit[1].index);
	res.push_back(std::make_pair(site_a,site_b));
      }
    }
    std::sort(res.begin(),res.end());
    res.erase(std::unique(res.begin(),res.end()),res.end());
    return res;
  }

  
}
#endif

