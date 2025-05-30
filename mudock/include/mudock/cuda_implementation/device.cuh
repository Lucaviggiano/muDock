#pragma once

#include <cuda_runtime.h>
#include <memory>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/grid.hpp>

namespace mudock {

  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point3D center_maps;
    const cudaStream_t stream;
    cudaTextureObject_t electro_tex, desolv_tex;
    cuda_wrapper<std::vector, cudaTextureObject_t> atom_texs;

    /// Vina protein data
    const int num_atoms;
    cuda_wrapper<std::vector, fp_type> protein_x;
    cuda_wrapper<std::vector, fp_type> protein_y;
    cuda_wrapper<std::vector, fp_type> protein_z;
    cuda_wrapper<std::vector, fp_type> p_vdw_radius;
    cuda_wrapper<std::vector, int> p_is_hbond_acceptor;
    cuda_wrapper<std::vector, int> p_is_hbond_donor;
    cuda_wrapper<std::vector, int> p_is_hydrophobic;

    device(const std::size_t gpu_id,
           std::shared_ptr<dynamic_molecule>& protein,
           std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
           std::shared_ptr<const grid_map>& electro_map,
           std::shared_ptr<const grid_map>& desolv_map);

    cudaStream_t get_stream() const;
  };
} // namespace mudock
