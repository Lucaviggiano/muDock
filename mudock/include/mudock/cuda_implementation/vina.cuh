#pragma once

#include <mudock/type_alias.hpp>
#include <thrust/device_ptr.h>

/// -10.219824: risultato finale variazione solo di +0.000009       <- with parallelisation
/// Score inter -15.313681, Score intra -1.478089, Score -10.219815 <- with no parallelisation
/// Score inter -15.313679, Score intra -1.478089, Score -10.219813 <- real

namespace mudock {

    __device__ static constexpr fp_type GAUSS1_COEFF_CUDA{- 0.035579};
    __device__ static constexpr fp_type GAUSS2_COEFF_CUDA{- 0.005156};
    __device__ static constexpr fp_type REPULSION_COEFF_CUDA{0.840245};
    __device__ static constexpr fp_type HYDROPHOBIC_COEFF_CUDA{- 0.035069};
    __device__ static constexpr fp_type H_BOND_COEFF_CUDA{- 0.587439};
    __device__ static constexpr fp_type NROT_COEFF_CUDA{0.05846};

    __device__ inline fp_type distance(fp_type x, fp_type y, fp_type z) {
        return sqrt( x*x + y*y + z*z );
    }

    __device__ inline fp_type gauss1(const size_t idx, const fp_type* __restrict__ dst_mtx) {
        fp_type gauss1 = 0;
        if(dst_mtx[idx] != 0) gauss1 = exp(- pow(dst_mtx[idx] / 0.5, 2));
        return gauss1;
    }

    __device__ inline fp_type gauss2(const size_t idx, const fp_type* __restrict__ dst_mtx) {
        fp_type gauss2 = 0;
        if(dst_mtx[idx] != 0) gauss2 = exp(- pow((dst_mtx[idx] - 3) / 2, 2));
        return gauss2;
    }

    __device__ inline fp_type repulsion(const size_t idx, const fp_type* __restrict__ dst_mtx) {
        return pow((dst_mtx[idx] < 0) * dst_mtx[idx], 2);
    }

    __device__ inline fp_type hydrophobic(const size_t idx, const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hydrophobic) {
        if(rec_lig_is_hydrophobic[idx] < 0) return 0; // Sentinel value, no interaction
        bool hydro_1 = rec_lig_is_hydrophobic[idx] && (dst_mtx[idx] <= 0.5);
        bool hydro_2_cond = rec_lig_is_hydrophobic[idx] && (dst_mtx[idx] > 0.5) && (dst_mtx[idx] < 1.5);
        fp_type hydro_2 = 1.5 * hydro_2_cond - hydro_2_cond * dst_mtx[idx];
        return hydro_1 + hydro_2;
    }

    __device__ inline fp_type hbonding(const size_t idx, const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hb) {
        if(rec_lig_is_hb[idx] < 0) return 0; // Sentinel value, no interaction
        bool h_bond_1 = rec_lig_is_hb[idx] && (dst_mtx[idx] <= -0.7);
        bool h_bond_2_cond = rec_lig_is_hb[idx] && (dst_mtx[idx] < 0) && (dst_mtx[idx] > -0.7);
        fp_type h_bond_2 = h_bond_2_cond * (- dst_mtx[idx]) / 0.7;
        return h_bond_1 + h_bond_2;
    }


    __device__ inline fp_type score_function(
        const fp_type* __restrict__ dst_mtx, 
        const int* __restrict__ ij_is_hydrophobic,
        const int* __restrict__ ij_is_hbond,
        const size_t size
    ) {

        fp_type g1 = 0; 
        fp_type g2 = 0; 
        fp_type rep = 0;
        fp_type hydro = 0;
        fp_type hbond = 0;

        for(size_t i = threadIdx.x; i < size; i += blockDim.x) {
            if(ij_is_hydrophobic[i] < 0) continue; // Sentinel value, no interaction
            g1  += gauss1(i, dst_mtx);
            g2  += gauss2(i, dst_mtx);
            rep += repulsion(i, dst_mtx);
            hydro += hydrophobic(i, dst_mtx, ij_is_hydrophobic);
            hbond += hbonding(i, dst_mtx, ij_is_hbond);
        }

        return GAUSS1_COEFF_CUDA * g1 + GAUSS2_COEFF_CUDA * g2 + REPULSION_COEFF_CUDA * rep + HYDROPHOBIC_COEFF_CUDA * hydro + H_BOND_COEFF_CUDA * hbond;
    }

    __device__ inline void parse_data(
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

        /// Output buffers
        fp_type* __restrict__ dst_mtx,
        int* __restrict__ is_hbond,
        int* __restrict__ is_hydrophobic
    ){

        fp_type vdw_sum;
        for (size_t proteinIdx = 0; proteinIdx < num_atoms_protein; proteinIdx++) {
            for(size_t ligandIdx = threadIdx.x; ligandIdx < num_atoms_ligand; ligandIdx += blockDim.x) {
                fp_type dst = distance(
                    (protein_x[proteinIdx] - ligand_x[ligandIdx]),
                    (protein_y[proteinIdx] - ligand_y[ligandIdx]),
                    (protein_z[proteinIdx] - ligand_z[ligandIdx])
                );

                size_t idx = ligandIdx + (proteinIdx * num_atoms_ligand);
                
                if(dst > 8){
                    is_hydrophobic[idx] = -1; /// Sentinel value. No interaction.
                }
                else{
                    vdw_sum = p_vdw_radius[proteinIdx] + l_vdw_radius[ligandIdx];
                    dst_mtx[idx] = dst - vdw_sum;
                    is_hbond[idx] = (p_is_hbond_acceptor[proteinIdx] && l_is_hbond_donor[ligandIdx]) || (l_is_hbond_acceptor[ligandIdx] && p_is_hbond_donor[proteinIdx]);
                    is_hydrophobic[idx] = p_is_hydrophobic[proteinIdx] && l_is_hydrophobic[ligandIdx];
                }

            }
        }
    }

    __device__ inline void parse_intra_data(
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const size_t num_interacting_pairs,

        /// Output buffers
        fp_type* __restrict__ intra_dst_mtx,
        int* __restrict__ intra_is_hbond,
        int* __restrict__ intra_is_hydrophobic
    ){

        fp_type vdw_sum;
        for(size_t i = threadIdx.x; i < num_interacting_pairs; i += blockDim.x) {
            int atom_1 = interacting_pairs_first[i];
            int atom_2 = interacting_pairs_second[i];

            fp_type dst = distance(
                (ligand_x[atom_1] - ligand_x[atom_2]),
                (ligand_y[atom_1] - ligand_y[atom_2]),
                (ligand_z[atom_1] - ligand_z[atom_2])
            );
            
            if(dst > 8){
                intra_is_hydrophobic[i] = -1; /// Sentinel value. No interaction.
            }
            else {
                vdw_sum = l_vdw_radius[atom_1] + l_vdw_radius[atom_2];
                intra_dst_mtx[i] = dst - vdw_sum;
                intra_is_hbond[i] = (l_is_hbond_acceptor[atom_1] && l_is_hbond_donor[atom_2]) || (l_is_hbond_acceptor[atom_2] && l_is_hbond_donor[atom_1]);
                intra_is_hydrophobic[i] = l_is_hydrophobic[atom_1] && l_is_hydrophobic[atom_2];
            }
        }
    }

    __device__ inline fp_type scoring_cuda(  
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
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const size_t num_interacting_pairs,

        /// Buffers (worst dim: num_atoms_ligand * num_atoms_protein)
        fp_type* __restrict__ dst_mtx,
        int* __restrict__ is_hbond,
        int* __restrict__ is_hydrophobic

    ){
        /// TODO: essere sicuri che tutti gli elementi siano diversi dall'idrogeno
        size_t mtx_size = 0;

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
            is_hbond, 
            is_hydrophobic
        );
        
        mtx_size = num_atoms_ligand * num_atoms_protein;
        fp_type inter_score = score_function(dst_mtx, is_hydrophobic, is_hbond, mtx_size);
        
        parse_intra_data(
            ligand_x,
            ligand_y,
            ligand_z,
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius,
            interacting_pairs_first,
            interacting_pairs_second,
            num_interacting_pairs,
            
            dst_mtx,
            is_hbond, 
            is_hydrophobic
        );

        mtx_size = num_interacting_pairs;
        fp_type intra_score = score_function(dst_mtx, is_hydrophobic, is_hbond, mtx_size);
        
        fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF_CUDA * active_torsions);
        
        // printf("dst_mtx len %ld\n", dst_mtx.size());
        // printf("intra_dst_mtx len %ld\n", intra_dst_mtx.size());
        // printf("Score inter %f, Score intra %f, Score %f\n", inter_score, intra_score, score); 
        
        return score;
    }
}
