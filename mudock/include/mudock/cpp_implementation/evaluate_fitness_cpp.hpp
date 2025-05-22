#pragma once

#include <mudock/chem/ligand_maps.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#include <stdint.h>

namespace mudock {
  void evaluate_fitness( const fp_type* __restrict__ ligand_x,
                    const fp_type* __restrict__ ligand_y,
                    const fp_type* __restrict__ ligand_z,
                    const fp_type* __restrict__ ligand_vol,
                    const fp_type* __restrict__ ligand_solpar,
                    const fp_type* __restrict__ ligand_charge,

                    /// Vina ligand data
                    const int* __restrict__ ligand_is_hbond_acceptor,
                    const int* __restrict__ ligand_is_hbond_donor,
                    const int* __restrict__ ligand_is_hydrophobic,
                    const fp_type* __restrict__ ligand_vdw_radius,
                    const int* __restrict__ interacting_pairs_first,
                    const int* __restrict__ interacting_pairs_second,
                    const size_t num_interacting_pairs,

                    const int* __restrict__ map_ligand_offsets,
                    const int num_atoms,
                    const int num_rotamers,
                    const int* __restrict__ frag_masks,
                    const int* __restrict__ frag_start_indexes,
                    const int* __restrict__ frag_stop_indexes,
                    const int num_nonbond,
                    const int* __restrict__ non_bond_list_a1,
                    const int* __restrict__ non_bond_list_a2,
                    const fp_type* __restrict__ cA_list,
                    const fp_type* __restrict__ cB_list,
                    const int* __restrict__ xB_list,
                    const fp_type* __restrict__ grid_maps,
                    const fp_type* __restrict__ electro_map,
                    const fp_type* __restrict__ desolv_map,

                    /// Vina protein data
                    const fp_type* __restrict__ protein_x,
                    const fp_type* __restrict__ protein_y,
                    const fp_type* __restrict__ protein_z,
                    const int* __restrict__ p_is_hbond_acceptor,
                    const int* __restrict__ p_is_hbond_donor,
                    const int* __restrict__ p_is_hydrophobic,
                    const fp_type* __restrict__ p_vdw_radius,
                    const int protein_num_atoms,

                    const int num_generations,
                    const int population_size,
                    const int tournament_length,
                    const fp_type mutation_prob,
                    const fp_type* __restrict__ minimum,
                    const fp_type* __restrict__ maximum,
                    const fp_type* __restrict__ center,
                    const int map_index_x,
                    const int map_index_xy,
                    individual* __restrict__ population_buffer1,
                    individual* __restrict__ population_buffer2,
                    const int seed);
} // namespace mudock
