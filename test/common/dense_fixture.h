#pragma once
#include <complex>
#include <vector>
using Scalar = std::complex<double>;
inline Scalar entry(int site, int physical, const std::vector<int>& bonds) {
    double x = 0.7 + 0.13 * site + 0.23 * physical;
    double y = 0.03 * (site + 1) * (physical ? -1 : 1);
    for (std::size_t k=0; k<bonds.size(); ++k) {
        x += (physical ? -0.17 : 0.19) * (k+1) * bonds[k];
        y += 0.07 * (k+1) * bonds[k];
    }
    return {x,y};
}
