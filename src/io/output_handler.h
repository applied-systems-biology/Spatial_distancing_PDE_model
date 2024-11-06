//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef ABM_OUTPUT_OUTPUTHANDLER_H_
#define ABM_OUTPUT_OUTPUTHANDLER_H_

#include <string>
#include <unordered_map>
#include "utils/io_util.h"

class Site;
class SimulationTime;

class OutputHandler {
 public:
  OutputHandler() = default;
  OutputHandler(const std::string &config_path, std::string project_dir);
  void concludeSimulation(int runs, int seed, double runtime) const;
  void setupOutputForRun(int run) const;
  void outputCurrentConfiguration(const Site &site,
                                  const SimulationTime &time,
                                  int run,
                                  int seed = 0,
                                  bool simulation_end = false,
                                  std::unordered_map<std::string, std::string> cmd_input_args = {}) const;

 private:
  void concludeSimulationRun(int current_run, int seed, const std::string &hash, std::unordered_map<std::string, std::string> cmd_input_args) const;
  std::string current_project_string_{};
  std::string current_project_dir_{};
  abm::util::OutputParameters parameters_{};
  abm::util::ConfigParameters config_param_{};
  mutable std::unordered_map<int, std::string> toc_dir_cache_;
};
#endif //ABM_OUTPUT_OUTPUTHANDLER_H_
