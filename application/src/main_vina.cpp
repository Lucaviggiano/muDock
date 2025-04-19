#include <stdio.h>
#include <vector>
#include <cmath>
#include <cassert>
#include <utility> 
#include <map>
#include <string>

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
        printf("G1 result %f\n", g1);
        fp_type g2  = gauss2(d_ij);
        fp_type rep = repulsion(d_ij);
        fp_type hydro = hydrophobic(d_ij, rec_lig_is_hydrophobic);
        printf("Hydro result: %f\n", hydro);
        fp_type hbond = hbonding(d_ij, rec_lig_is_hbond);
        printf("Hbond result: %f\n", hbond);

        return GAUSS1_COEFF * g1 + GAUSS2_COEFF * g2 + REPULSION_COEFF * rep + HYDROPHOBIC_COEFF * hydro + H_BOND_COEFF * hbond;
    }

    std::vector<fp_type> ignore_distant_atoms(std::vector<fp_type> dist_mtx) {
        std::vector<fp_type> vina_dist_mtx = std::vector<fp_type>();
        for(size_t i = 0; i < dist_mtx.size(); i++){
            if(dist_mtx[i] <= 8) vina_dist_mtx.push_back(dist_mtx[i]);
        }
        // printf("vina_dist_mtx len %ld\n", vina_dist_mtx.size());
        return vina_dist_mtx;
    }

    void print_mtx(std::vector<fp_type> mtx){
        for(size_t i = 0; i < mtx.size(); i++){
            printf("%f, ", mtx[i]);
        }
        printf("\n");
    }

    std::vector<fp_type> generate_pldist_mtrx(
        const std::vector<fp_type> rec_heavy_atoms_x, 
        const std::vector<fp_type> rec_heavy_atoms_y,
        const std::vector<fp_type> rec_heavy_atoms_z,
        const std::vector<fp_type> lig_pose_heavy_atoms_coords_x,
        const std::vector<fp_type> lig_pose_heavy_atoms_coords_y,
        const std::vector<fp_type> lig_pose_heavy_atoms_coords_z
    ){
        std::vector<fp_type> dst_mtx = std::vector<fp_type>(rec_heavy_atoms_x.size() * lig_pose_heavy_atoms_coords_x.size(), 0);

        for(size_t i = 0; i < rec_heavy_atoms_x.size(); i++){
                for(size_t j = 0; j < lig_pose_heavy_atoms_coords_x.size(); j++){
                dst_mtx[FLATTENED_2D(j, i, lig_pose_heavy_atoms_coords_x.size())] = sqrt(
                    pow(rec_heavy_atoms_x[i] - lig_pose_heavy_atoms_coords_x[j], 2) +
                    pow(rec_heavy_atoms_y[i] - lig_pose_heavy_atoms_coords_y[j], 2) +
                    pow(rec_heavy_atoms_z[i] - lig_pose_heavy_atoms_coords_z[j], 2)
                );
            }
        }
        return dst_mtx;
    }

    std::vector<fp_type> generate_intra_mtrx(
        const std::vector<fp_type> lig_pose_heavy_atoms_coords_x,
        const std::vector<fp_type> lig_pose_heavy_atoms_coords_y,
        const std::vector<fp_type> lig_pose_heavy_atoms_coords_z,
        const std::vector<std::pair<size_t, size_t>> lig_intra_interacting_pairs
    ){
        std::vector<fp_type> intra_dst_mtx = std::vector<fp_type>();

        for(size_t i = 0; i < lig_intra_interacting_pairs.size(); i++){
            
            size_t atom_i = lig_intra_interacting_pairs[i].first;
            size_t atom_j = lig_intra_interacting_pairs[i].second; 

            intra_dst_mtx.push_back(sqrt(
                pow(lig_pose_heavy_atoms_coords_x[atom_i] - lig_pose_heavy_atoms_coords_x[atom_j], 2) +
                pow(lig_pose_heavy_atoms_coords_y[atom_i] - lig_pose_heavy_atoms_coords_y[atom_j], 2) +
                pow(lig_pose_heavy_atoms_coords_z[atom_i] - lig_pose_heavy_atoms_coords_z[atom_j], 2)
            ));
        }
        return intra_dst_mtx;
    }

    fp_type scoring(  
        const std::vector<fp_type> dst_mtx,
        const std::vector<fp_type> intra_dst_mtx,

        const std::vector<fp_type> rec_lig_atom_vdw_sum,
        const std::vector<fp_type> intra_rec_lig_atom_vdw_sum,
                        
        const size_t active_torsion,
        const size_t inactive_torsion,

        /// Inter
        const std::vector<bool> rec_lig_is_hydrophobic,
        const std::vector<bool> rec_lig_is_hbond,
        /// Intra
        const std::vector<bool> intra_rec_lig_is_hydrophobic,
        const std::vector<bool> intra_rec_lig_is_hbond
    ){

        printf("dst_mtx len %ld\n", dst_mtx.size());
        printf("intra_dst_mtx len %ld\n", intra_dst_mtx.size());
        
        fp_type inter_score = score_function(dst_mtx, rec_lig_atom_vdw_sum, rec_lig_is_hydrophobic, rec_lig_is_hbond);
        fp_type intra_score = score_function(intra_dst_mtx, intra_rec_lig_atom_vdw_sum, intra_rec_lig_is_hydrophobic, intra_rec_lig_is_hbond);
        printf("Score inter %f\n", inter_score); 
        printf("Score intra %f\n", intra_score);
        
        return (inter_score + intra_score) / ( 1 + NROT_COEFF * (active_torsion + 0.5 * inactive_torsion));
    }

} 

using namespace OpenBabel;
using namespace std;

bool isHydrophobicAtom(OBMol &mol, OBAtom *atom) {
    // Define a simple hydrophobic SMARTS: non-polar carbon, sp3
    OBSmartsPattern smarts;
    smarts.Init(
        "[c,s,Br,I,S&H0&v2,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6X4H0]);+0]"); // aliphatic carbon not bonded to N/O/S
  
    if (!smarts.Match(mol))
      return false;
  
    // Check if the atom is part of any match
    for (const auto &match: smarts.GetMapList()) {
      for (uint idx: match) {
        if (idx == atom->GetIdx())
          return true;
      }
    }
  
    return false;
}

int main(int argc, char* argv[]) {

    /// print(", ".join(map(str, self.receptor.rec_heavy_atoms_xyz[:, 0].tolist())))

    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_rec.pdb> <input_lig.pdb>" << endl;
        return 1;
    }

    std::filesystem::path pdb_rec{argv[1]};
    const mudock::ob_mol_wrapper ob_mol_rec = mudock::parser(pdb_rec);

    std::filesystem::path pdb_lig{argv[2]};
    const mudock::ob_mol_wrapper ob_mol_lig = mudock::parser(pdb_lig);

    size_t active_torsion = 11;
    size_t inactive_torsion = 0;

    std::vector<std::pair<size_t, size_t>> lig_intra_interacting_pairs;
    std::vector<mudock::fp_type> intra_dst_mtx;
    std::vector<bool> intra_is_hbond;
    std::vector<bool> intra_is_hydro;
    std::vector<mudock::fp_type> intra_rec_lig_atom_vdw_sum;

    for (auto atom_it = ob_mol_lig->BeginBonds(); atom_it < ob_mol_lig->EndBonds(); ++atom_it){
        const auto bond = *atom_it;

        const auto atom1 = bond->GetBeginAtom();
        const auto atom2 = bond->GetEndAtom();

        mudock::fp_type dst = bond->GetLength() /// Tested is equal to the distance

        if(dst > 8) continue;

        intra_dst_mtx.push_back(dst);

        bool isA1Donor    = atom1->IsHbondDonor();
        bool isA1Acceptor = atom1->IsHbondAcceptor();

        bool isA2Donor    = atom2->IsHbondDonor();
        bool isA2Acceptor = atom2->IsHbondAcceptor();

        bool isA1Hydro    = isHydrophobicAtom(*ob_mol_lig, atom1);
        bool isA2Hydro    = isHydrophobicAtom(*ob_mol_lig, atom2); 

        mudock::fp_type vdwA1 = OBElements::GetVdwRad(atom1->GetAtomicNum());
        mudock::fp_type vdwA2 = OBElements::GetVdwRad(atom2->GetAtomicNum());

        intra_is_hbond.push_back((isA1Donor && isA2Acceptor) || (isA2Donor && isA1Acceptor));
        intra_is_hydro.push_back(isA1Hydro && isA2Hydro);
        intra_rec_lig_atom_vdw_sum.push_back(vdwA1 + vdwA2);

        const auto atom1_idx = atom1->GetIdx() - 1;
        const auto atom2_idx = atom2->GetIdx() - 1;
        lig_intra_interacting_pairs.push_back(std::make_pair(atom1_idx, atom2_idx));
    }

    std::vector<mudock::fp_type> dst_mtx;
    std::vector<bool> is_hydro;
    std::vector<bool> is_hbond;
    std::vector<mudock::fp_type> rec_lig_atom_vdw_sum;

    for (auto atom_it_rec = ob_mol_rec->BeginAtoms(); atom_it_rec < ob_mol_rec->EndAtoms(); ++atom_it_rec) {
        for (auto atom_it_lig = ob_mol_lig->BeginAtoms(); atom_it_lig < ob_mol_lig->EndAtoms(); ++atom_it_lig) {
            const auto atom_rec = *atom_it_rec;
            const auto atom_lig = *atom_it_lig;

            mudock::fp_type xr = atom_rec->GetX();
            mudock::fp_type yr = atom_rec->GetY();
            mudock::fp_type zr = atom_rec->GetZ();

            mudock::fp_type xl = atom_lig->GetX();
            mudock::fp_type yl = atom_lig->GetY();
            mudock::fp_type zl = atom_lig->GetZ();

            mudock::fp_type dst =  sqrt(
                    pow(xr - xl, 2) +
                    pow(yr - yl, 2) +
                    pow(zr - zl, 2)
                );

            if(dst > 8) continue;

            dst_mtx.push_back(dst);

            bool isRecHydro    = isHydrophobicAtom(*ob_mol_rec, atom_rec);
            bool isRecDonor    = atom_rec->IsHbondDonor();
            bool isRecAcceptor = atom_rec->IsHbondAcceptor();
            mudock::fp_type vdwRec = OBElements::GetVdwRad(atom_rec->GetAtomicNum());
        
            bool isLigHydro    = isHydrophobicAtom(*ob_mol_lig, atom_lig);
            bool isLigDonor    = atom_lig->IsHbondDonor();
            bool isLigAcceptor = atom_lig->IsHbondAcceptor();
            mudock::fp_type vdwLig = OBElements::GetVdwRad(atom_lig->GetAtomicNum());

            is_hydro.push_back(isRecHydro && isLigHydro);
            is_hbond.push_back((isRecDonor && isLigAcceptor) || (isRecAcceptor && isLigDonor));
            rec_lig_atom_vdw_sum.push_back(vdwRec + vdwLig); /// 0.6 max difference from opendock
        }
    }

    /// current results
    /*
        G1 result 111.071415
        Hydro result: 2.000000
        Hbond result: 0.000000
        G1 result 0.000110 dovrebbe essere ~21
        Hydro result: 0.000000
        Hbond result: 0.000000
        Score inter -15.313697
        Score intra 137.946701
    */

    /// Opendock -12, -1, -8
    std::printf("Score: %f\n", mudock::scoring(  dst_mtx,
                                            intra_dst_mtx,
                                            rec_lig_atom_vdw_sum,
                                            intra_rec_lig_atom_vdw_sum,
                                            active_torsion, 
                                            inactive_torsion,
                                            is_hydro, 
                                            is_hbond, 
                                            intra_is_hydro, 
                                            intra_is_hbond
                                            ));

    return 0;
}
