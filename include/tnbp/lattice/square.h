/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/lattice/square.h
@brief define the bonds and pairs for two-dimensional square lattice
*/
#ifndef TNBP_LATTICE_SQUARE_H
#define TNBP_LATTICE_SQUARE_H

namespace tnbp {

  /**
     Function to define the bonds for sqaure lattice
   */
  std::vector<std::pair<int,int>> bond_square_lattice(int Lx, int Ly) {
    size_t N = 2*Lx*Ly-Lx-Ly;
    std::vector<std::pair<int,int>> res(N);
    size_t m=0;
    for(int iy=0; iy < Ly; iy++) {
      for(int ix=0; ix < Lx; ix++) {
	int site_i = ix + Lx * iy;
	int jx = ix+1;
	int jy = iy+1;
	if( jx < Lx ) {
	  int site_j = jx + Lx * iy;
	  res[m++] = std::make_pair(site_i,site_j);
	}
	if( jy < Ly ) {
	  int site_j = ix + Lx * jy;
	  res[m++] = std::make_pair(site_i,site_j);
	}
      }
    }
    return res;
  }

}

#endif
