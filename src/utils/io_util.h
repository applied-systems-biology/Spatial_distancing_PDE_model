//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef ABM_UTILS_IO_UTIL_H_
#define ABM_UTILS_IO_UTIL_H_

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <array>

#include <type_traits>


#include "basic/Coordinate3D.h"
namespace abm::util {

struct ConfigParameters {
  int runs{};
  int number_of_threads{};
  int system_seed{};
  std::string config_path{};
  std::string output_dir{};
  std::string input_dir{};
  std::string conda_dir{};
  std::unordered_map<std::string, std::vector<std::string>> screening_parameters{};
  int screen_start_idx{};
};

struct OutputParameters {
  bool activated{};
  bool output_graphs{};
  bool output_ovito{};
};

struct SimulationParameters {
  struct MorphologyParameters {
    double radius{};
    std::string type{};
  };
  struct MoleculeMoleculeInteractionParameters{
    std::string first{};
    std::string second{};
    double binding{};
    double unbinding{};
  };

  struct Flow{
      double amplitude{};
      double frequence{};
      double target_conc{};
  };

  struct MoleculeParameters{
    double initial_concentration{};
    double diffusion_coefficient{};
    double decay{};
    Flow boundary_flow{};
    double spontaneous_rate{};
    std::string name{};
  };

  struct MoleculeInteractionParameters{
    double secretion{};
    double secretion_delay{};
    bool secretion_variable_over_time;
    double uptake{};
    double inside_conc_decay{}; // decay for conc inside the cell
  };


  struct AgentParameters {
    int number{};
    int initial_distribution{};
    double pos_dist{};
    std::string input_distribution_path{};
    std::string type{};
    MorphologyParameters morphology_parameters{};
    std::vector<std::pair<std::string,MoleculeInteractionParameters>> molecule_interactions{};
    bool track_inside_conc{};
  };

  struct AgentManagerParameters {
    std::string site_identifier{};
    std::vector<std::shared_ptr<AgentParameters>> agents;
  };

  struct MoleculeManagerParameters{
      std::vector<std::shared_ptr<MoleculeParameters>> molecules;
      std::vector<std::pair<std::string, MoleculeMoleculeInteractionParameters>> molecule_molecule_interactions{};
      std::string boundary_condition{};
  };

struct StoppingCriteria{
    std::string molecule;
    double threshold;
    double time;
    double init_time;
};
  struct SiteParameters {
    std::string identifier{};
    std::string type{};
    AgentManagerParameters agent_manager_parameters{};
    MoleculeManagerParameters molecule_manager_parameters{};
    bool molecular_layer_{};
   StoppingCriteria stopping_criteria{};

  };
  struct CuboidSiteParameters : public SiteParameters {
    Coordinate3D upper_bound{};
    Coordinate3D lower_bound{};
    std::tuple<int, int, int> molecules_grid_size{};
  };

  int dimensions{};
  double max_time{};
  double time_stepping{};
  bool dt_auto{};

  std::string topic{};
  std::unique_ptr<SiteParameters> site_parameters;
  std::unordered_map<std::string, std::string> cmd_input_args{};
};

struct AnalyserParameters {
  std::string species{};
  std::string analyser_path{};

  std::unordered_set<std::string> active_measurements{};
  std::vector<std::string> cell_state_count{};
  double output_every_X_time_unit{};
};


struct InputParameters {
  struct DefaultRateParameters {
    double rate{};
    std::string type{};
    std::string key{};
    // for conditional rate
    std::string condition{};
  };
    struct MoleculeRateParameters : public DefaultRateParameters {
        std::string molecule_name{};
    };
  std::vector<std::unique_ptr<DefaultRateParameters>> rates{};
};

InputParameters getInputParameters(const std::string &input_config);
AnalyserParameters getAnalyserParameters(const std::string &analyser_config);
ConfigParameters getMainConfigParameters(const std::string &config_path);
OutputParameters getOutputParameters(const std::string &output_config);
SimulationParameters getSimulationParameters(const std::string &simulator_config);

void executeShellCommand(const std::string &command, bool suppress_output = true);

} // namespace abm::util
#endif//ABM_UTILS_IO_UTIL_H_
