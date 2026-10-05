#include "link_check.h"
#include <iostream>

int main() {
    // External inline functions denote one entity across translation units.
    // This also catches replacing inline with static to hide link failures.
    if (addresses_a() != addresses_b()) {
        std::cerr << "FAIL: different function identities across translation units\n";
        return 1;
    }
    int first=0,last=5;
    tnbp::get_range(2,1,first,last);
    if(first!=3 || last!=5) return 1;
    std::cout << "PASS: header-only link and 11 shared function identities\n";
}
