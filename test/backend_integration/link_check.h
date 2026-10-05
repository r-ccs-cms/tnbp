#pragma once
#include "qasm/utility.h"
#include "qasm/any.h"
#include "pauli/pauli_string.h"
#include "tnbp/tnbp.h"
#include <tuple>

using FunctionAddresses = std::tuple<
    decltype(&tnbp::get_range),
    decltype(&tnbp::bond_oned_lattice),
    decltype(&tnbp::parallel_bond_oned_lattice),
    decltype(&tnbp::bond_honeycomb_lattice),
    decltype(&tnbp::parallel_bond_honeycomb_lattice),
    decltype(&tnbp::Bond_HeavyHexLattice),
    decltype(&tnbp::ParallelBond_HeavyHexLattice),
    decltype(&tnbp::bond_square_lattice),
    decltype(&tnbp::SitesFromQasm),
    decltype(&tnbp::EdgesFromQasm),
    decltype(static_cast<void (*)(std::string&, int, MPI_Comm)>(&tnbp::MpiBcast))>;
FunctionAddresses addresses_a();
FunctionAddresses addresses_b();
