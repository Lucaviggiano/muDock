#include <memory>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/calc_energy_cpp.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/evaluate_fitness_cpp.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/virtual_screen.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>

#include <mudock/cpp_implementation/vina.hpp>

#define TRANSFORM 1

namespace mudock {

 
  virtual_screen_cpp::virtual_screen_cpp(std::shared_ptr<dynamic_molecule>& _protein, const knobs& knobs){

    protein = _protein;
    auto prt = *_protein;
    #if VINA
      auto tp  = mudock::point3D{};
      auto itp = mudock::index3D{1, 1, 1};
      auto sv  = std::vector<mudock::grid_atom_map>{};
      for (int i = 0; i < mudock::num_ligand_map_types(); ++i)
        sv.emplace_back(mudock::autodock_type_from_map(static_cast<mudock::ligand_map_types>(i)), prt);
      grid_atom_maps    = std::make_shared<const mudock::grid_atom_mapper>(sv);
      electro_map = std::make_shared<const mudock::grid_map>(prt);
      desolv_map   = std::make_shared<const mudock::grid_map>(prt);
    #else
      grid_atom_maps = std::make_shared<const grid_atom_mapper>(generate_atom_grid_maps(ptr));
      electro_map = std::make_shared<const grid_map>(generate_electrostatic_grid_map(ptr));
      desolv_map = std::make_shared<const grid_map>(generate_desolvation_grid_map(ptr));
    #endif
    population =  std::vector<individual>(knobs.population_number);;
    next_population =  std::vector<individual>(knobs.population_number);;
    configuration = knobs;
  }

  void virtual_screen_cpp::operator()(static_molecule& ligand) {
    const auto seed =
    static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    
    // Place the molecule to the center of the target protein
    const int num_atoms = ligand.num_atoms();
    const auto x = ligand.get_x(), y = ligand.get_y(), z = ligand.get_z();

    #if TRANSFORM
    const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
    translate_molecule(x.data(),
                       y.data(),
                       z.data(),
                       num_atoms,
                       electro_map->center.x - ligand_center_of_mass.x,
                       electro_map->center.y - ligand_center_of_mass.y,
                       electro_map->center.z - ligand_center_of_mass.z);
    #endif

    const auto lig_is_hbond_ac = ligand.get_is_hbond_acceptor();
    const auto lig_is_hbond_dn = ligand.get_is_hbond_donor();
    const auto lig_is_hydro = ligand.get_is_hydrophobic();
    const auto lig_vdw_radius = ligand.get_vdw_radius();

    /// Parse the protein data
    const auto protein_x = protein.get()->get_x();
    const auto protein_y = protein.get()->get_y();
    const auto protein_z = protein.get()->get_z();
    const auto protein_is_hbond_ac = protein.get()->get_is_hbond_acceptor();
    const auto protein_is_hbond_dn = protein.get()->get_is_hbond_donor();
    const auto protein_is_hydro = protein.get()->get_is_hydrophobic();
    const auto protein_vdw_radius = protein.get()->get_vdw_radius();
    const auto protein_num_atoms = protein.get()->num_atoms();

    /// Parse interacting pairs of the ligand
    auto [ip_first, ip_second] = get_interactive_pairs(ligand);
    size_t num_interacting_pairs = ip_first.size();

    // Find out the rotatable bonds in the ligand
    auto graph = make_graph(ligand.get_bonds(), ligand.num_atoms());
    const auto ligand_fragments =
        std::make_unique<fragments<static_containers>>(graph, ligand.get_bonds(), ligand.num_atoms());

    const auto num_rotamers = ligand_fragments.get()->get_num_rotatable_bonds();

    // Get weed bonds and non bonds lists
    grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
    nonbonds(nbmatrix, ligand.get_bonds(), num_atoms);
    std::vector<int> non_bond_list_a1, non_bond_list_a2;
    weed_bonds(nbmatrix, non_bond_list_a1, non_bond_list_a2, num_atoms, *ligand_fragments.get());
    const auto non_bond_size = non_bond_list_a1.size();
    std::vector<fp_type> cA_v, cB_v;
    std::vector<int> xB_v;
    precompute_lennard_jones(non_bond_size, cA_v, cB_v, xB_v, ligand, non_bond_list_a1, non_bond_list_a2);

    const fp_type minimum[3]        = {electro_map.get()->minimum_coord.x,
                                       electro_map.get()->minimum_coord.y,
                                       electro_map.get()->minimum_coord.z};
    const fp_type maximum[3]        = {electro_map.get()->maximum_coord.x,
                                       electro_map.get()->maximum_coord.y,
                                       electro_map.get()->maximum_coord.z};
    const fp_type center[3]         = {electro_map.get()->center.x,
                                       electro_map.get()->center.y,
                                       electro_map.get()->center.z};
    const int atom_map_size         = grid_atom_maps.get()->get_single_map_size();
    const fp_type* atom_map_pointer = grid_atom_maps.get()->get_fused_maps().data();

    std::vector<int> frag_masks;
    std::vector<int> frag_start_indexes;
    std::vector<int> frag_stop_indexes;
    frag_masks.resize(num_atoms * num_rotamers);
    frag_start_indexes.resize(num_rotamers);
    frag_stop_indexes.resize(num_rotamers);
    for (int rot = 0; rot < num_rotamers; ++rot) {
      std::memcpy((frag_masks.data() + num_atoms * rot),
                  ligand_fragments.get()->get_mask(rot).data(),
                  num_atoms * sizeof(int));
      const auto [start_index, stop_index] = ligand_fragments.get()->get_rotatable_atoms(rot);
      frag_start_indexes.data()[rot]       = start_index;
      frag_stop_indexes.data()[rot]        = stop_index;
    }

    std::vector<int> map_ligand_offsets;
    map_ligand_offsets.resize(num_atoms);
    for (int i = 0; i < num_atoms; i++)
      map_ligand_offsets[i] =
          static_cast<int>(map_from_autodock_type(ligand.autodock_type(i))) * atom_map_size;
    // Simulate the population evolution for the given amount of time
    evaluate_fitness(x.data(),
                     y.data(),
                     z.data(),
                     ligand.get_vol().data(),
                     ligand.get_solpar().data(),
                     ligand.get_charge().data(),

                     /// Vina ligand data
                     lig_is_hbond_ac.data(),
                     lig_is_hbond_dn.data(),
                     lig_is_hydro.data(),
                     lig_vdw_radius.data(),
                     ip_first.data(),
                     ip_second.data(),
                     num_interacting_pairs,

                     map_ligand_offsets.data(),
                     num_atoms,
                     ligand_fragments.get()->get_num_rotatable_bonds(),
                     frag_masks.data(),
                     frag_start_indexes.data(),
                     frag_stop_indexes.data(),
                     non_bond_size,
                     non_bond_list_a1.data(),
                     non_bond_list_a2.data(),
                     cA_v.data(),
                     cB_v.data(),
                     xB_v.data(),
                     atom_map_pointer,
                     electro_map.get()->data(),
                     desolv_map.get()->data(),

                     /// Vina protein data
                     protein_x.data(),
                     protein_y.data(),
                     protein_z.data(),
                     protein_is_hbond_ac.data(),
                     protein_is_hbond_dn.data(),
                     protein_is_hydro.data(),
                     protein_vdw_radius.data(),
                     protein_num_atoms,

                     configuration.num_generations,
                     configuration.population_number,
                     configuration.tournament_length,
                     configuration.mutation_prob,
                     minimum,
                     maximum,
                     center,
                     electro_map.get()->index.size_x(),
                     electro_map.get()->index.size_xy(),
                     population.data(),
                     next_population.data(),
                     seed);

    // update the ligand position with the best one that we found
    const std::vector<individual>& last_population =
        (configuration.num_generations % 2 == 0) ? next_population : population;

    const auto best_individual_it =
        std::min_element(std::begin(last_population),
                         std::end(last_population),
                         [](const auto a, const auto b) { return a.score < b.score; });

    #if TRANSFORM
    apply(x.data(),
          y.data(),
          z.data(),
          best_individual_it->genes,
          num_atoms,
          num_rotamers,
          frag_masks.data(),
          frag_start_indexes.data(),
          frag_stop_indexes.data());
    #endif

    ligand.properties.assign(property_type::SCORE, std::to_string(best_individual_it->score));
  }

} // namespace mudock
