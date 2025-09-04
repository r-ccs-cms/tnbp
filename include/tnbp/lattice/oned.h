//// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/lattice/oned.h
@brief define the bonds and pairs for one-dimensional lattice
*/
#ifndef TNBP_LATTICE_ONED_H
#define TNBP_LATTICE_ONED_H

namespace {

  /**
     Function to define the bonds for one-dimensional lattice
   */
  std::vector<std::pair<int,int>> Bond_OneDLattice(int Lx, int Px) {
    size_t L = Lx;
    size_t N = Lx-1;
    if( Px != 0 ) {
      N += 1;
    }
    std::vector<std::pair<int,int>> res(N);
    size_t m = 0;
    for(int i=0; i < L; i++) {
      int site_i = i;
      int site_j = i+1;
      if( Px != 0 ) {
        site_j = (i+1) % Lx;
      }
      if( site_j < L ) {
        res[m] = std::make_pair(site_i,site_j);
        m++;
      }
    }
    return res;
  }

  /**
     Function to define the parallel bonds for one-dimensional lattice
   */
  std::vector<std::vector<std::pair<int,int>>> ParallelBond_OneDLattice(int Lx, int Px) {
    std::vector<std::vector<std::pair<int,int>>> res;
    size_t size_A;
    size_t size_B;
    size_t size_C = 0;
    if( Lx % 2 == 0 ) {
      // Lx = 4 (OBC):   * -- * -- * -- *
      // A-bond:         * -- *    * -- *
      // B-bond:              * -- *    *
      // Lx = 4 (PBC): - * -- * -- * -- * -
      // A-bond:         * -- *    * -- *
      // B-bond:       - *    * -- *    * -
      size_A = Lx/2;
      size_B = Lx/2-1;
      if( Px != 0 ) {
	size_B += 1;
      }
      res.resize(2);
    } else {
      // Lx = 5 (OBC):   * -- * -- * -- * -- *
      // A-bond:         * -- *    * -- *
      // B-bond:              * -- *    * -- *

      // Lx = 5 (PBC): - * -- * -- * -- * -- * -
      // A-bond:         * -- *    * -- *
      // B-bond:              * -- *    * -- *
      // C-bond:       - *                   * -
      
      // Lx = 7:   * -- * -- * -- * -- * -- * -- *
      // A-bond:   * -- *    * -- *    * -- *
      // B-bond:        * -- *    * -- *    * -- *
      size_A = Lx/2;
      size_B = Lx/2;
      if( Px != 0 ) {
	size_C = 1;
	res.resize(3);
      } else {
	res.resize(2);
      }
    }

    res[0].resize(size_A);
    res[1].resize(size_B);
    if( size_C != 0 ) {
      res[2].resize(size_C);
    }

    size_t mA = 0;
    size_t mB = 0;
    size_t mC = 0;

    for(int i=0; i < Lx; i++) {
      int site_i = i;
      int site_j = i+1;
      if( Px != 0 ) {
	site_j = (i+1) % Lx;
      }
      if( site_j < Lx ) {
	if( site_i % 2 == 0 ) {
	  if ( site_j == 0 && size_C != 0 ) {
	    res[2][mC] = std::make_pair(site_i,site_j);
	    mC++;
	  } else {
	    res[0][mA] = std::make_pair(site_i,site_j);
	    mA++;
	  }
	} else {
	  res[1][mB] = std::make_pair(site_i,site_j);
	  mB++;
	}
      }
    }
    return I;
  }
  
  
}

#endif
