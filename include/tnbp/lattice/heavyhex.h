//// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/lattice/heavyhex.h
@brief define the bonds and pairs for heavy-hexagonal lattice
*/
#ifndef TNBP_LATTICE_HEAVYHEX_H
#define TNBP_LATTICE_HEAVYHEX_H

namespace tnbp {

  std::vector<std::pair<int,int>> Bond_HeavyHexLattice(int Lx, int Ly) {

    size_t NumV = 2 * Lx + 1 + (2*Lx+2)*(Ly-1) + 2 * Lx + 1;
    size_t NumE = 2 * Lx + (2*Lx+1)*(Ly-1) + 2 * Lx + (Lx+1)*Ly;    
    size_t L = NumV + NumE;
    size_t N = 4 * Lx + ( 4 * Lx +2 ) * (Ly-1) + 4 * Lx + 2*(Lx+1) * Ly;
    std::vector<std::pair<int,int>> I(N);
    
    std::vector<int> x_begin(Ly+1);
    std::vector<int> x_end(Ly+1);
    for(int y=0; y < Ly; y++) {
      if( y == 0 ) {
	x_begin[y] = 0;
	x_end[y] = 4*Lx+1;
      } else {
	x_begin[y] = x_end[y-1] + Lx+1;
	x_end[y] = x_end[y-1] + 5 * Lx + 4;
      }
    }
    x_begin[Ly] = x_end[Ly-1] + Lx+1;
    x_end[Ly] = x_begin[Ly] + 4*Lx+1;

    size_t m=0;
    for(int y=0; y < Ly+1; y++) {
      for(int x=x_begin[y]; x < x_end[y]-1; x++) {
	I[m] = std::make_pair(x,x+1);
	m++;
      }
      if( y < Ly ) {
	int x_margin = 2;
	if( y == 0 ) {
	  x_margin = 0;
	}
	for(int x=0; x < Lx+1; x++) {
	  I[m] = std::make_pair(x_begin[y]+x_margin+4*x,x_end[y]+x);
	  m++;
	  I[m] = std::make_pair(x_end[y]+x,x_begin[y+1]+4*x);
	  m++;
	}
      }
    }
    return I;
  }

  std::vector<std::vector<std::pair<int,int>>> ParallelBond_HeavyHexLattice(int Lx, int Ly) {
    std::vector<int> x_begin(Ly+1);
    std::vector<int> x_end(Ly+1);
    for(int y=0; y < Ly; y++) {
      if( y == 0 ) {
	x_begin[y] = 0;
	x_end[y] = 4*Lx+1;
      } else {
	x_begin[y] = x_end[y-1] + Lx+1;
	x_end[y] = x_end[y-1] + 5 * Lx + 4;
      }
    }
    x_begin[Ly] = x_end[Ly-1] + Lx+1;
    x_end[Ly] = x_begin[Ly] + 4*Lx+1;

    std::vector<std::vector<std::pair<int,int>>> I(3);

    int x_c = 0;
    int x_margin = 0;
    for(int y=0; y < Ly+1; y++) {
      for(int x=x_begin[y]; x < x_end[y]-1; x++) {
	int color = ( x + x_c ) % 3;
	I[color].push_back(std::make_pair(x,x+1));
      }
      if( y < Ly ) {
	if( y == 0 ) {
	  x_margin = 0;
	} else {
	  x_margin = 2;
	}
	for(int x=0; x < Lx+1; x++) {
	  int color = ( x_begin[y]+x_margin+4*x + x_c + 1 ) % 3;
	  I[color].push_back(std::make_pair(x_begin[y]+x_margin+4*x,x_end[y]+x));
	  color = (color+1) % 3;
	  I[color].push_back(std::make_pair(x_end[y]+x,x_begin[y+1]+4*x));
	  if( x == 0 ) {
	    if( (x_begin[y+1])%3 == 0 ) {
	      x_c = (color+1)%3;
	    } else if ( (x_begin[y+1])%3 == 1 ) {
	      x_c = color%3;
	    } else {
	      x_c = (color+2)%3;
	    }
	  }
	}
      }
    }
    return I;
  }
  
  
}

#endif
  
