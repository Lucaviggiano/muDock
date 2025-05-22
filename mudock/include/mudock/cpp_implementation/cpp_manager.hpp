#pragma once

#include <memory>
#include <mudock/compute.hpp>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  // this function will configure and create (if needed) cpp workers to the threadpool
  void manage_cpp(const std::vector<std::string>& configurations,
                  threadpool& pool,
                  std::shared_ptr<dynamic_molecule>& protein,
                  const knobs knobs,
                  std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                  std::shared_ptr<safe_stack<static_molecule>>& output_molecules);
} // namespace mudock
