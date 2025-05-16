#include <stdio.h>
#include <vector>
#include <cmath>
#include <cassert>
#include <utility> 
#include <map>
#include <string>
#include <queue>
#include <set>

#include <filesystem>
#include <iomanip>
#include <iostream>

#include <mudock/format.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/chem/elements.hpp>

#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/parsmart.h>


#define GAUSS1_COEFF        (- 0.035579)
#define GAUSS2_COEFF        (- 0.005156)
#define REPULSION_COEFF     (0.840245)
#define HYDROPHOBIC_COEFF   (- 0.035069)
#define H_BOND_COEFF        (- 0.587439)
#define NROT_COEFF          (0.05846)

#define FLATTENED_2D(x, y, index_x) ((y) * index_x + (x))

namespace mudock {

    fp_type distance(fp_type x1, fp_type y1, fp_type z1, fp_type x2, fp_type y2, fp_type z2) {
        return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
    }

    fp_type gauss1(const std::vector<fp_type> dst_mtx) {
        fp_type gauss1 = 0;
        for( size_t i = 0; i < dst_mtx.size(); i++) {
            if(dst_mtx[i] != 0) gauss1 += exp(- pow(dst_mtx[i] / 0.5, 2));
        }
        return gauss1;
    }

    fp_type gauss2(const std::vector<fp_type> dst_mtx) {
        fp_type gauss2 = 0;
        for( size_t i = 0; i < dst_mtx.size(); i++) {
            if(dst_mtx[i] != 0) gauss2 += exp(- pow((dst_mtx[i] - 3) / 2, 2));
        }
        return gauss2;
    }

    fp_type repulsion(const std::vector<fp_type> dst_mtx) {
        fp_type repulsion = 0;
        for( size_t i = 0; i < dst_mtx.size(); i++) {
            repulsion += pow((dst_mtx[i] < 0) * dst_mtx[i], 2);
        }
        return repulsion;
    }

    fp_type hydrophobic(const std::vector<fp_type> dst_mtx, const std::vector<bool> rec_lig_is_hydrophobic) {
        fp_type hydrophobic = 0;
        for( size_t i = 0; i < dst_mtx.size(); i++) {
            bool hydro_1 = rec_lig_is_hydrophobic[i] && (dst_mtx[i] <= 0.5);
            bool hydro_2_cond = rec_lig_is_hydrophobic[i] && (dst_mtx[i] > 0.5) && (dst_mtx[i] < 1.5);
            fp_type hydro_2 = 1.5 * hydro_2_cond - hydro_2_cond * dst_mtx[i];
            hydrophobic += hydro_1 + hydro_2;
        }
        return hydrophobic;
    }

    fp_type hbonding(const std::vector<fp_type> dst_mtx, const std::vector<bool> rec_lig_is_hb) {
        fp_type h_bonding = 0;
	
        for( size_t i = 0; i < dst_mtx.size(); i++) {
            bool h_bond_1 = rec_lig_is_hb[i] && (dst_mtx[i] <= -0.7);
            bool h_bond_2_cond = rec_lig_is_hb[i] && (dst_mtx[i] < 0) && (dst_mtx[i] > -0.7);
            fp_type h_bond_2 = h_bond_2_cond * (- dst_mtx[i]) / 0.7;
            h_bonding += h_bond_1 + h_bond_2;
        }
        return h_bonding;
    }

    fp_type score_function(
        const std::vector<fp_type> dst_mtx, 
        const std::vector<fp_type> rec_lig_atom_vdw_sum,
        const std::vector<bool> rec_lig_is_hydrophobic,
        const std::vector<bool> rec_lig_is_hbond
    ) {

        std::vector<fp_type> d_ij = std::vector<fp_type>(dst_mtx.size(), 0);

        for( size_t i = 0; i < dst_mtx.size(); i++) {
            d_ij[i] = dst_mtx[i] - rec_lig_atom_vdw_sum[i];
        }

        fp_type g1  = gauss1(d_ij);
        fp_type g2  = gauss2(d_ij);
        fp_type rep = repulsion(d_ij);
        fp_type hydro = hydrophobic(d_ij, rec_lig_is_hydrophobic);
        fp_type hbond = hbonding(d_ij, rec_lig_is_hbond);

        return GAUSS1_COEFF * g1 + GAUSS2_COEFF * g2 + REPULSION_COEFF * rep + HYDROPHOBIC_COEFF * hydro + H_BOND_COEFF * hbond;
    }


    void print_mtx(std::vector<fp_type> mtx){
        for(size_t i = 0; i < mtx.size(); i++){
            printf("%f, ", mtx[i]);
        }
        printf("\n");
    }

    void parse_data(
        /// Protein data
        const size_t num_atoms_protein,
        const fp_type* __restrict__ protein_x,
        const fp_type* __restrict__ protein_y,
        const fp_type* __restrict__ protein_z,
        const int* __restrict__ p_is_hbond_acceptor,
        const int* __restrict__ p_is_hbond_donor,
        const int* __restrict__ p_is_hydrophobic,
        const fp_type* __restrict__ p_vdw_radius,

        ///Ligand data
        const size_t num_atoms_ligand,
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,

        std::vector<fp_type>& dst_mtx,
        std::vector<fp_type>& rec_lig_atom_vdw_sum,     
        std::vector<bool>& rec_lig_is_hbond, 
        std::vector<bool>& rec_lig_is_hydrophobic   
    ){

        for(size_t proteinIdx = 0; proteinIdx < num_atoms_protein; proteinIdx++){
            for(size_t ligandIdx = 0; ligandIdx < num_atoms_ligand; ligandIdx++){
                fp_type dst = sqrt(
                    pow(protein_x[proteinIdx] - ligand_x[ligandIdx], 2) +
                    pow(protein_y[proteinIdx] - ligand_y[ligandIdx], 2) +
                    pow(protein_z[proteinIdx] - ligand_z[ligandIdx], 2)
                );

                if(dst > 8) continue;
                dst_mtx.emplace_back(dst);
                rec_lig_atom_vdw_sum.emplace_back(
                    p_vdw_radius[proteinIdx] + l_vdw_radius[ligandIdx]
                );
                rec_lig_is_hbond.emplace_back(
                    (p_is_hbond_acceptor[proteinIdx] && l_is_hbond_donor[ligandIdx]) || (l_is_hbond_acceptor[ligandIdx] && p_is_hbond_donor[proteinIdx])
                );
                rec_lig_is_hydrophobic.emplace_back(
                    p_is_hydrophobic[proteinIdx] && l_is_hydrophobic[ligandIdx]
                );
            }
        }
    }

    void parse_intra_data(
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const std::vector<std::pair<int, int>> interacting_pairs,

        std::vector<fp_type>& intra_dst_mtx,
        std::vector<fp_type>& intra_rec_lig_atom_vdw_sum,     
        std::vector<bool>& intra_rec_lig_is_hbond, 
        std::vector<bool>& intra_rec_lig_is_hydrophobic   
    ){

        printf("Parsing intra data\n");

        for(size_t i = 0; i < interacting_pairs.size(); i++){
            int atom_1 = interacting_pairs[i].first;
            int atom_2 = interacting_pairs[i].second;

            fp_type dst = sqrt(
                pow(ligand_x[atom_1] - ligand_x[atom_2], 2) +
                pow(ligand_y[atom_1] - ligand_y[atom_2], 2) +
                pow(ligand_z[atom_1] - ligand_z[atom_2], 2)
            );

            if(dst > 8) continue;
            
            intra_dst_mtx.push_back(dst);
            intra_rec_lig_atom_vdw_sum.push_back(
                l_vdw_radius[atom_1] + l_vdw_radius[atom_2]
            );
            intra_rec_lig_is_hbond.push_back(
                (l_is_hbond_acceptor[atom_1] && l_is_hbond_donor[atom_2]) || (l_is_hbond_acceptor[atom_2] && l_is_hbond_donor[atom_1])
            );
            intra_rec_lig_is_hydrophobic.push_back(
                l_is_hydrophobic[atom_1] && l_is_hydrophobic[atom_2]
            );
        }
    }

    fp_type scoring(  
        /// Protein data
        const size_t num_atoms_protein,
        const fp_type* __restrict__ protein_x,
        const fp_type* __restrict__ protein_y,
        const fp_type* __restrict__ protein_z,
        const int* __restrict__ p_is_hbond_acceptor,
        const int* __restrict__ p_is_hbond_donor,
        const int* __restrict__ p_is_hydrophobic,
        const fp_type* __restrict__ p_vdw_radius,

        ///Ligand data
        const size_t num_atoms_ligand,
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const size_t active_torsions,
        const std::vector<std::pair<int, int>> interacting_pairs
    ){

        std::vector<fp_type> dst_mtx                  = std::vector<fp_type>();
        std::vector<fp_type> rec_lig_atom_vdw_sum     = std::vector<fp_type>();
        std::vector<bool> rec_lig_is_hbond            = std::vector<bool>();
        std::vector<bool> rec_lig_is_hydrophobic      = std::vector<bool>();

        std::vector<fp_type> intra_dst_mtx                  = std::vector<fp_type>();
        std::vector<fp_type> intra_rec_lig_atom_vdw_sum     = std::vector<fp_type>();
        std::vector<bool> intra_rec_lig_is_hbond            = std::vector<bool>();
        std::vector<bool> intra_rec_lig_is_hydrophobic      = std::vector<bool>();

        parse_data(
            num_atoms_protein,
            protein_x,
            protein_y,
            protein_z,
            p_is_hbond_acceptor,
            p_is_hbond_donor,
            p_is_hydrophobic,
            p_vdw_radius,

            num_atoms_ligand,
            ligand_x,
            ligand_y,
            ligand_z,
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius,

            dst_mtx, 
            rec_lig_atom_vdw_sum,     
            rec_lig_is_hbond, 
            rec_lig_is_hydrophobic
        );

        parse_intra_data(
            ligand_x,
            ligand_y,
            ligand_z,
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius,
            interacting_pairs,

            intra_dst_mtx,
            intra_rec_lig_atom_vdw_sum,     
            intra_rec_lig_is_hbond, 
            intra_rec_lig_is_hydrophobic   
        );

        printf("dst_mtx len %ld\n", dst_mtx.size());
        printf("intra_dst_mtx len %ld\n", intra_dst_mtx.size());
        fp_type intra_score = score_function(intra_dst_mtx, intra_rec_lig_atom_vdw_sum, intra_rec_lig_is_hydrophobic, intra_rec_lig_is_hbond);
        fp_type inter_score = score_function(dst_mtx, rec_lig_atom_vdw_sum, rec_lig_is_hydrophobic, rec_lig_is_hbond);
        printf("Score inter %f\n", inter_score); 
        printf("Score intra %f\n", intra_score);
        
        return (inter_score + intra_score) / ( 1 + NROT_COEFF * active_torsions);
    }

} 

using namespace OpenBabel;
using namespace std;

// Funzione per trovare la distanza in numero di bond tra due atomi
int shortest_bond_path(OBMol& mol, int startIdx, int endIdx) {
    std::vector<bool> visited(mol.NumAtoms() + 1, false);
    std::queue<std::pair<OBAtom*, int>> q;

    OBAtom* startAtom = mol.GetAtom(startIdx);
    q.push({startAtom, 0});
    visited[startIdx] = true;

    while (!q.empty()) {
        auto [currentAtom, depth] = q.front();
        q.pop();

        if (currentAtom->GetIdx() == endIdx)
            return depth;

        OBBondIterator bondIt;
        OBBond* bond = currentAtom->BeginBond(bondIt);
        while (bond) {
            OBAtom* neighbor = bond->GetNbrAtom(currentAtom);
            int neighborIdx = neighbor->GetIdx();

            if (!visited[neighborIdx]) {
                visited[neighborIdx] = true;
                q.push({neighbor, depth + 1});
            }

            bond = currentAtom->NextBond(bondIt);
        }
    }

    return -1; // Not connected
}

std::vector<std::pair<int, int>> get_distant_atom_pairs(OBMol& mol, int minBondDistance = 4) {
    std::vector<std::pair<int, int>> distant_pairs;
    int numAtoms = mol.NumAtoms();

    for (int i = 1; i <= numAtoms; ++i) {
        for (int j = i + 1; j <= numAtoms; ++j) {
            int bond_distance = shortest_bond_path(mol, i, j);
            if (bond_distance >= minBondDistance) {
                printf("[%d, %d], ", i - 1, j - 1);
                distant_pairs.emplace_back(i - 1, j - 1); /// TODO: riguarda questo
            }
        }
    }

    printf("\n----------------------------------------------------------------\n");

    return distant_pairs;
}

inline auto same_fragment(const std::span<const int>& mask, const int& atom_id1, const int& atom_id2) {
    return mask[atom_id1] == mask[atom_id2] || (mask[atom_id1] == 0 && mask[atom_id2] == 2) ||
           (mask[atom_id2] == 0 && mask[atom_id1] == 2) || (mask[atom_id1] == 1 && mask[atom_id2] == 3) ||
           (mask[atom_id2] == 1 && mask[atom_id1] == 3);
}

/*
Opendock:
    [0, 1, 2]
    [3]
    [4]
    [5, 6, 7, 8, 9, 10, 11]
    [12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    [22]
    [23, 24, 25]
    [26]
    [27]
    [28, 29, 30, 31, 32, 33, 34]
    [35, 36, 37, 38, 39, 40]
    [41, 42, 43]

Mudock:
    Fragment 10: 0 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
    Fragment 9: 0 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 
    Fragment 8: 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
    Fragment 7: 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 
    Fragment 6: 33 35 36 37 38 39 40 41 42 43 
    Fragment 5: 37 41 42 43 
    Fragment 4: 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
    Fragment 3: 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 
    Fragment 2: 9 12 13 14 15 16 17 18 19 20 21 
    Fragment 1: 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 
    Fragment 0: 23 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 
*/

void printAtomsInFragment(const std::unordered_map<int, std::vector<int>>& atoms_in_fragment) {
    for (const auto& [fragment_id, atom_indices] : atoms_in_fragment) {
        std::cout << "Fragment " << fragment_id << ": ";
        for (int atom_index : atom_indices) {
            std::cout << atom_index << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<int> get_lig_frag_mask(const mudock::static_molecule &ligand){

    auto graph = make_graph(ligand.get_bonds(), ligand.num_atoms());
    const auto ligand_fragments =
        std::make_unique<mudock::fragments<mudock::static_containers>>(graph,
                                                                       ligand.get_bonds(),
                                                                       ligand.num_atoms());
    std::vector<int> frag_masks;
    std::vector<int> frag_start_indexes;
    std::vector<int> frag_stop_indexes;
    size_t num_atoms = ligand.num_atoms();
    size_t num_rotamers = ligand.num_rotamers();
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
      printf("Start idx: %d, Stop idx: %d\n", frag_start_indexes[rot], frag_stop_indexes[rot]);
    }

    return frag_masks;
}

std::unordered_map<int, std::vector<int>> get_atoms_in_frag(
    const int* __restrict__ frag_masks,  
    const size_t num_atoms,
    const size_t num_rotamers
){
    std::unordered_map<int, std::vector<int>> atoms_in_fragment;
    
    for (size_t rot = 0; rot < num_rotamers; ++rot) {
        const auto* bitmask_rot = frag_masks + rot * num_atoms;
        std::vector<int> atoms_rot;
        for(size_t atom = 0; atom < num_atoms; ++atom){
            if(bitmask_rot[atom] != 0){
                atoms_rot.push_back(atom);
            }
        }
        atoms_in_fragment[rot] = atoms_rot;
    }

    printAtomsInFragment(atoms_in_fragment);
    return atoms_in_fragment;
}


std::vector<std::pair<int, int>> get_interactive_pairs(
    const int* __restrict__ frag_masks, 
    const size_t num_atoms, 
    const size_t num_rotamers, 
    std::vector<std::pair<int, int>> potential_interacting
){
    
    std::set<std::pair<int, int>> unique_out;

    std::unordered_map<int, std::vector<int>> atoms_in_fragment = get_atoms_in_frag(frag_masks, num_atoms, num_rotamers);
                                                                       
    for (size_t rot1 = 0; rot1 < num_rotamers; ++rot1) {

        const auto* bitmask_rot1 = frag_masks + rot1 * num_atoms;
        std::vector<int> atoms_rot1 = atoms_in_fragment[rot1];
        for(size_t i = 0; i < atoms_rot1.size(); ++i){

            for(size_t rot2 = rot1 + 1; rot2 <num_rotamers; ++rot2){

                const auto* bitmask_rot2 = frag_masks + rot2 * num_atoms;
                std::vector<int> atoms_rot2 = atoms_in_fragment[rot2];
                for(size_t j = 0; j < atoms_rot2.size(); ++j){

                    int atom1 = atoms_rot1[i];
                    int atom2 = atoms_rot2[j];
                    
                    /// Assicurati che il pair non sia già presente
                    bool found = false;
                    for(const auto& pair : potential_interacting) {
                        if(atom1 == pair.first && atom2 == pair.second || atom1 == pair.second && atom2 == pair.first){
                            found = true;
                            break;
                        }
                    }

                    if(!found) continue;

                    int mask_1 = bitmask_rot1[atom1];
                    int mask_2 = bitmask_rot2[atom2];

                    ///if(mask_1 == 3 || mask_2 == 3 || mask_1 == 2 || mask_2 == 2) continue;

                    /// Order pair before adding it to the output
                    int a = std::min(atom1, atom2);
                    int b = std::max(atom1, atom2);
                    if (unique_out.insert({a, b}).second) {
                        printf("[%d, %d], ", a, b);
                    }
                }
            }
        }
    }

    printf("\nOld pairs size: %ld\n", potential_interacting.size());
    printf("New pairs size: %ld\n", unique_out.size());

    std::vector<std::pair<int, int>> out(unique_out.begin(), unique_out.end());
    return out;
}

int main(int argc, char* argv[]) {

    /// print(", ".join(map(str, self.receptor.rec_heavy_atoms_xyz[:, 0].tolist())))

    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_rec.pdb> <input_lig.pdb>" << endl;
        return 1;
    }

    /// TODO: essere sicuri che tutti gli elementi siano diversi dall'idrogeno

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

    std::vector<std::pair<int, int>> p_interactive_pairs = get_distant_atom_pairs(*ob_mol_lig, 4);

    std::vector<std::pair<int, int>> interactive_pairs = get_interactive_pairs(
        get_lig_frag_mask(ligand).data(),
        ligand.num_atoms(),
        ligand.num_rotamers(), 
        p_interactive_pairs);

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
                                            interactive_pairs
                                            ));

    return 0;
}
