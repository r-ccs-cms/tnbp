//// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/lattice/honeycomb.h
@brief define the bonds and pairs for honeycomb lattice
*/
#ifndef TNBP_LATTICE_HONEYCOMB_H
#define TNBP_LATTICE_HONEYCOMB_H

namespace {

  /**
     Function to define the bonds for honeycomb lattice
   */
  std::vector<std::pair<int,int>> Bond_HoneycombLattice(int Lx, int Ly) {
    int L = 2 * Lx * Ly;
    int N = Lx*Ly+Lx*(Ly-1)+(Lx-1)*Ly;
    std::vector<std::pair<int,int>> res(N);
    size_t m=0;
    for(int y=0; y < Ly; y++) {
      for(int x=0; x < Lx; x++) {
        int site_a = 0 + 2*x + 2*Lx*y;
        int site_b = 1 + 2*x + 2*Lx*y;
        int site_x = 0 + 2*(x+1) + 2*Lx*y;
        int site_y = 0 + 2*x + 2*Lx*(y+1);
        res[m] = std::make_pair(site_a,site_b);
        m++;
        if( x < Lx-1 ) {
          res[m] = std::make_pair(site_b,site_x);
          m++;
        }
        if( y < Ly-1 ) {
          res[m] = std::make_pair(site_b,site_y);
          m++;
        }
      }
    }
    return res;
  }

  /**
     Function to define the parallel bonds for honeycomb lattice
   */
  std::vector<std::vector<std::pair<int,int>>> ParallelBond_HoneycombLattice(int Lx, int Ly) {

    std::vector<std::vector<std::pair<int,int>>> res(3);
    size_t size_A = Lx*Ly;
    size_t size_B = (Lx-1)*Ly;
    size_t size_C = Lx*(Ly-1);
    res[0].resize(size_A);
    res[1].resize(size_B);
    res[2].resize(size_C);

    size_t mA = 0;
    size_t mB = 0;
    size_t mC = 0;

    for(int y=0; y < Ly; y++) {
      for(int x=0; x < Lx; x++) {
        int site_a = 0 + 2*x + 2*Lx*y;
        int site_b = 1 + 2*x + 2*Lx*y;
        int site_x = 0 + 2*(x+1) + 2*Lx*y;
        int site_y = 0 + 2*x + 2*Lx*(y+1);
        res[0][mA] = std::make_pair(site_a,site_b);
        mA++;
        if( x < Lx-1 ) {
          res[1][mB] = std::make_pair(site_b,site_x);
          mB++;
        }
        if( y < Ly-1 ) {
          res[2][mC] = std::make_pair(site_b,site_y);
          mC++;
        }
      }
    }
    assert( size_A == mA );
    assert( size_B == mB );
    assert( size_C == mC );
    return res;
  }
  
  
}

#endif
