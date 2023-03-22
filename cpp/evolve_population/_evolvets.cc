// Wright-Fisher simulation for a fwdpy11::DiploidPopulation with
// tree sequences.
//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <core/evolve_discrete_demes/evolvets.hpp>
//#include <fwdpy11/discrete_demography/simulation/demographic_model_state.hpp>

namespace py = pybind11;

void
init_evolve_with_tree_sequences(py::module &m)
{
    py::class_<evolve_with_tree_sequences_options>(m,
                                                   "_evolve_with_tree_sequences_options")
        .def(py::init<>())
        .def_readwrite("preserve_selected_fixations",
                       &evolve_with_tree_sequences_options::preserve_selected_fixations)
        .def_readwrite("suppress_edge_table_indexing",
                       &evolve_with_tree_sequences_options::suppress_edge_table_indexing)
        .def_readwrite("record_gvalue_matrix",
                       &evolve_with_tree_sequences_options::record_gvalue_matrix)
        .def_readwrite(
            "track_mutation_counts_during_sim",
            &evolve_with_tree_sequences_options::track_mutation_counts_during_sim)
        .def_readwrite(
            "remove_extinct_mutations_at_finish",
            &evolve_with_tree_sequences_options::remove_extinct_mutations_at_finish)
        .def_readwrite("reset_treeseqs_to_alive_nodes_after_simplification",
                       &evolve_with_tree_sequences_options::
                           reset_treeseqs_to_alive_nodes_after_simplification)
        .def_readwrite("preserve_first_generation",
                       &evolve_with_tree_sequences_options::preserve_first_generation)
        .def_readwrite("allow_residual_selfing",
                       &evolve_with_tree_sequences_options::allow_residual_selfing);

    m.def("evolve_with_tree_sequences", &evolve_with_tree_sequences);
}
