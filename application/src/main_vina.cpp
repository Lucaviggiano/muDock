#include <stdio.h>
#include <mudock/type_alias.hpp>
#include <vector>
#include <cmath>

// #include "mudock/src/cpp_implementation/vina.cpp"


#define GAUSS1_COEFF        (- 0.035579)
#define GAUSS2_COEFF        (- 0.005156)
#define REPULSION_COEFF     (0.840245)
#define HYDROPHOBIC_COEFF   (- 0.035069)
#define H_BOND_COEFF        (- 0.587439)
#define NROT_COEFF          (0.05846)


#define FLATTENED_2D(x, y, index_x) ((y) * index_x + (x))

namespace mudock
{
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
    
        fp_type hydrophobic(const std::vector<fp_type> dst_mtx, const std::vector<fp_type> rec_lig_is_hydrophobic) {
            fp_type hydrophobic = 0;
            for( size_t i = 0; i < dst_mtx.size(); i++) {
                fp_type hydro_1 = rec_lig_is_hydrophobic[i] * (dst_mtx[i] <= 0.5);
                fp_type hydro_2_cond = rec_lig_is_hydrophobic[i] * (dst_mtx[i] > 0.5) * (dst_mtx[i] < 1.5);
                fp_type hydro_2 = 1.5 * hydro_2_cond - hydro_2_cond * dst_mtx[i];
                hydrophobic += hydro_1 + hydro_2;
            }
            return hydrophobic;
        }
    
        fp_type hbonding(const std::vector<fp_type> dst_mtx, const std::vector<fp_type> rec_lig_is_hb) {
            fp_type h_bonding = 0;
            for( size_t i = 0; i < dst_mtx.size(); i++) {
                fp_type h_bond_1 = rec_lig_is_hb[i] * (dst_mtx[i] <= -0.7);
                fp_type h_bond_2 = rec_lig_is_hb[i] * (dst_mtx[i] < 0) * (dst_mtx[i] > -0.7) * 1.0 * (- dst_mtx[i]) / 0.7;
                h_bonding += h_bond_1 + h_bond_2;
            }
            return h_bonding;
        }
    
        fp_type score_function(
            const std::vector<fp_type> dst_mtx, 
            const std::vector<fp_type> rec_lig_atom_vdw_sum,
            const std::vector<fp_type> rec_lig_is_hydrophobic,
            const std::vector<fp_type> rec_lig_is_hbond
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

    /*
    *    Args:
    *        dist_matrix [N, M]: the distance matrix with less than 8 angstroms.
    *        N is the number of poses,
    *    M is the number of rec-lig atom pairs less than 8 Angstroms in each pose.
    *
    *    Returns:
    *        final_inter_score [N, 1]
    */

    fp_type scoring(  
                        /// Inter
                        const std::vector<fp_type> dst_mtx, 
                        const std::vector<fp_type> rec_lig_atom_vdw_sum,
                        const std::vector<fp_type> rec_lig_is_hydrophobic,
                        const std::vector<fp_type> rec_lig_is_hbond,
                        /// Intra
                        const std::vector<fp_type> intra_dst_mtx, 
                        const std::vector<fp_type> intra_rec_lig_atom_vdw_sum,
                        const std::vector<fp_type> intra_rec_lig_is_hydrophobic,
                        const std::vector<fp_type> intra_rec_lig_is_hbond,
                        /// Others
                        const fp_type active_torsion,
                        const fp_type inactive_torsion
    ){
        fp_type inter_score = score_function(dst_mtx, rec_lig_atom_vdw_sum, rec_lig_is_hydrophobic, rec_lig_is_hbond);
        fp_type intra_score = score_function(intra_dst_mtx, intra_rec_lig_atom_vdw_sum, intra_rec_lig_is_hydrophobic, intra_rec_lig_is_hbond);
        printf("Score inter %f\n", inter_score); 
        printf("Score inter %f\n", intra_score);
        
        return (inter_score + intra_score) / ( 1 + NROT_COEFF * (active_torsion + 0.5 * inactive_torsion));
    }
}
