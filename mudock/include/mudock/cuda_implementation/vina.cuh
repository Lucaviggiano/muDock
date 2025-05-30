#pragma once

#include <mudock/type_alias.hpp>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

namespace mudock {

    __device__ static constexpr fp_type GAUSS1_COEFF_CUDA{- 0.035579};
    __device__ static constexpr fp_type GAUSS2_COEFF_CUDA{- 0.005156};
    __device__ static constexpr fp_type REPULSION_COEFF_CUDA{0.840245};
    __device__ static constexpr fp_type HYDROPHOBIC_COEFF_CUDA{- 0.035069};
    __device__ static constexpr fp_type H_BOND_COEFF_CUDA{- 0.587439};
    __device__ static constexpr fp_type NROT_COEFF_CUDA{0.05846};

    __device__ inline fp_type gauss1(const fp_type* __restrict__ dst_mtx, int size) {
        fp_type gauss1 = 0;
        for( size_t i = 0; i < size; i++) {
            if(dst_mtx[i] != 0) gauss1 += exp(- pow(dst_mtx[i] / 0.5, 2));
        }
        return gauss1;
    }

    __device__ inline fp_type gauss2(const fp_type* __restrict__ dst_mtx, int size) {
        fp_type gauss2 = 0;
        for( size_t i = 0; i < size; i++) {
            if(dst_mtx[i] != 0) gauss2 += exp(- pow((dst_mtx[i] - 3) / 2, 2));
        }
        return gauss2;
    }

    __device__ inline fp_type repulsion(const fp_type* __restrict__ dst_mtx, int size) {
        fp_type repulsion = 0;
        for( size_t i = 0; i < size; i++) {
            repulsion += pow((dst_mtx[i] < 0) * dst_mtx[i], 2);
        }
        return repulsion;
    }

    __device__ inline fp_type hydrophobic(const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hydrophobic, int size) {
        fp_type hydrophobic = 0;
        for( size_t i = 0; i < size; i++) {
            bool hydro_1 = rec_lig_is_hydrophobic[i] && (dst_mtx[i] <= 0.5);
            bool hydro_2_cond = rec_lig_is_hydrophobic[i] && (dst_mtx[i] > 0.5) && (dst_mtx[i] < 1.5);
            fp_type hydro_2 = 1.5 * hydro_2_cond - hydro_2_cond * dst_mtx[i];
            hydrophobic += hydro_1 + hydro_2;
        }
        return hydrophobic;
    }

    __device__ inline fp_type hbonding(const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hb, int size) {
        fp_type h_bonding = 0;
        for( size_t i = 0; i < size; i++) {
            bool h_bond_1 = rec_lig_is_hb[i] && (dst_mtx[i] <= -0.7);
            bool h_bond_2_cond = rec_lig_is_hb[i] && (dst_mtx[i] < 0) && (dst_mtx[i] > -0.7);
            fp_type h_bond_2 = h_bond_2_cond * (- dst_mtx[i]) / 0.7;
            h_bonding += h_bond_1 + h_bond_2;
        }
        return h_bonding;
    }


    __device__ inline fp_type score_function(
        fp_type* __restrict__ dst_mtx, 
        const fp_type* __restrict__ rec_lig_atom_vdw_sum,
        const int* __restrict__ rec_lig_is_hydrophobic,
        const int* __restrict__ rec_lig_is_hbond,
        const size_t size
    ) {

        for( size_t i = 0; i < size; i++) {
            dst_mtx[i] = dst_mtx[i] - rec_lig_atom_vdw_sum[i];
        }

        fp_type g1  = gauss1(dst_mtx, size);
        fp_type g2  = gauss2(dst_mtx, size);
        fp_type rep = repulsion(dst_mtx, size);
        fp_type hydro = hydrophobic(dst_mtx, rec_lig_is_hydrophobic, size);
        fp_type hbond = hbonding(dst_mtx, rec_lig_is_hbond, size);

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
        fp_type* __restrict__ atom_vdw_sum,
        int* __restrict__ is_hbond,
        int* __restrict__ is_hydrophobic,
        int* mtx_size 
    ){

        int idx = 0;
        for(size_t proteinIdx = 0; proteinIdx < num_atoms_protein; proteinIdx++){
            for(size_t ligandIdx = 0; ligandIdx < num_atoms_ligand; ligandIdx++){
                fp_type dst = sqrt(
                    pow(protein_x[proteinIdx] - ligand_x[ligandIdx], 2) +
                    pow(protein_y[proteinIdx] - ligand_y[ligandIdx], 2) +
                    pow(protein_z[proteinIdx] - ligand_z[ligandIdx], 2)
                );

                if(dst > 8) continue;

                dst_mtx[idx] = dst;
                atom_vdw_sum[idx] = p_vdw_radius[proteinIdx] + l_vdw_radius[ligandIdx];
                is_hbond[idx] = (p_is_hbond_acceptor[proteinIdx] && l_is_hbond_donor[ligandIdx]) || (l_is_hbond_acceptor[ligandIdx] && p_is_hbond_donor[proteinIdx]);
                is_hydrophobic[idx] = p_is_hydrophobic[proteinIdx] && l_is_hydrophobic[ligandIdx];
                idx++;
            }
        }
        *mtx_size = idx;
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
        fp_type* __restrict__ intra_atom_vdw_sum,
        int* __restrict__ intra_is_hbond,
        int* __restrict__ intra_is_hydrophobic,
        int* intra_mtx_size
    ){

        int idx = 0;
        for(size_t i = 0; i < num_interacting_pairs; i++){
            int atom_1 = interacting_pairs_first[i];
            int atom_2 = interacting_pairs_second[i];

            fp_type dst = sqrt(
                pow(ligand_x[atom_1] - ligand_x[atom_2], 2) +
                pow(ligand_y[atom_1] - ligand_y[atom_2], 2) +
                pow(ligand_z[atom_1] - ligand_z[atom_2], 2)
            );

            if(dst > 8) continue;
            
            intra_dst_mtx[idx] = dst;
            intra_atom_vdw_sum[idx] = l_vdw_radius[atom_1] + l_vdw_radius[atom_2];
            intra_is_hbond[idx] = (l_is_hbond_acceptor[atom_1] && l_is_hbond_donor[atom_2]) || (l_is_hbond_acceptor[atom_2] && l_is_hbond_donor[atom_1]);
            intra_is_hydrophobic[idx] = l_is_hydrophobic[atom_1] && l_is_hydrophobic[atom_2];
            idx++;
        }
        *intra_mtx_size = idx;
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
        fp_type* __restrict__ atom_vdw_sum,
        int* __restrict__ is_hbond,
        int* __restrict__ is_hydrophobic

    ){
        /// TODO: essere sicuri che tutti gli elementi siano diversi dall'idrogeno
        int mtx_size;
        
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
            atom_vdw_sum,     
            is_hbond, 
            is_hydrophobic,
            &mtx_size
        );
        
        fp_type inter_score = score_function(dst_mtx, atom_vdw_sum, is_hydrophobic, is_hbond, mtx_size);
        
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
            atom_vdw_sum,     
            is_hbond, 
            is_hydrophobic,
            &mtx_size  
        );

        fp_type intra_score = score_function(dst_mtx, atom_vdw_sum, is_hydrophobic, is_hbond, mtx_size);
        
        fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF_CUDA * active_torsions);
        
        // printf("dst_mtx len %ld\n", dst_mtx.size());
        // printf("intra_dst_mtx len %ld\n", intra_dst_mtx.size());
        printf("Score inter %f, Score intra %f, Score %f\n", inter_score, intra_score, score); 
        
        return score;
    }
}