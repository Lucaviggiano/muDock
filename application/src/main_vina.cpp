#include <stdio.h>
#include <vector>

#include <mudock/format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/cpp_implementation/vina.hpp>

int main(int argc, char* argv[]) {

    if (argc != 3) {
        printf("Usage: %s <input_rec.pdb> <input_lig.pdb>\n", argv[0]);
        return 1;
    }

    std::filesystem::path pdb_rec{argv[1]};
    const mudock::ob_mol_wrapper ob_mol_rec = mudock::parser(pdb_rec);
    auto protein_ptr = std::make_unique<mudock::dynamic_molecule>();
    auto& protein    = *protein_ptr;
    mudock::convert<mudock::rotate_check>(protein, ob_mol_rec);

    std::filesystem::path pdb_lig{argv[2]};
    const mudock::ob_mol_wrapper ob_mol_lig = mudock::parser(pdb_lig);
    auto ligand_ptr = std::make_unique<mudock::static_molecule>();
    auto& ligand    = *ligand_ptr;
    mudock::convert<mudock::rotate_check>(ligand, ob_mol_lig);

    std::vector<std::pair<int, int>> interacting_pairs = get_interactive_pairs(ligand);
    size_t num_interacting_pairs = interacting_pairs.size();
    std::vector<int> interacting_pairs_first(num_interacting_pairs);
    std::vector<int> interacting_pairs_second(num_interacting_pairs);
    for(int i = 0; i < num_interacting_pairs; i++){
        interacting_pairs_first[i] = interacting_pairs[i].first;
        interacting_pairs_second[i] = interacting_pairs[i].second;
    }

    std::printf("Score: %f\n", mudock::scoring(  
                                            protein.num_atoms(),
                                            protein.get_x().data(),
                                            protein.get_y().data(),
                                            protein.get_z().data(),
                                            protein.get_is_hbond_acceptor().data(),
                                            protein.get_is_hbond_donor().data(),
                                            protein.get_is_hydrophobic().data(),
                                            protein.get_vdw_radius().data(),
                                            ligand.num_atoms(),
                                            ligand.get_x().data(),
                                            ligand.get_y().data(),
                                            ligand.get_z().data(),
                                            ligand.get_is_hbond_acceptor().data(),
                                            ligand.get_is_hbond_donor().data(),
                                            ligand.get_is_hydrophobic().data(),
                                            ligand.get_vdw_radius().data(),
                                            ligand.num_rotamers(),
                                            interacting_pairs_first.data(),
                                            interacting_pairs_second.data(),
                                            num_interacting_pairs
                                            ));

    return 0;
}
