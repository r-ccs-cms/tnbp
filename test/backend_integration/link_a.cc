#include "link_check.h"

FunctionAddresses addresses_a() {
    return {
        &tnbp::get_range,
        &tnbp::bond_oned_lattice,
        &tnbp::parallel_bond_oned_lattice,
        &tnbp::bond_honeycomb_lattice,
        &tnbp::parallel_bond_honeycomb_lattice,
        &tnbp::Bond_HeavyHexLattice,
        &tnbp::ParallelBond_HeavyHexLattice,
        &tnbp::bond_square_lattice,
        &tnbp::SitesFromQasm,
        &tnbp::EdgesFromQasm,
        static_cast<void (*)(std::string&, int, MPI_Comm)>(&tnbp::MpiBcast)};
}
