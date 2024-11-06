//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "utils/io_util.h"

#include "external/json.hpp"
#include "utils/macros.h"
#include <optional>
#include <random>
#include <optional>

namespace abm::util {

using json = nlohmann::json;

ConfigParameters getMainConfigParameters(const std::string &main_config) {
  ConfigParameters parameters{};
  std::ifstream json_file(main_config);
  json json_parameters;
  json_file >> json_parameters;
  parameters.runs = json_parameters["Spatial_Distancing"].value("runs", 1);
  parameters.config_path = json_parameters["Spatial_Distancing"].value("config_path", "");
  parameters.output_dir = json_parameters["Spatial_Distancing"].value("output_dir", "/tmp/results");
  parameters.conda_dir = json_parameters["Spatial_Distancing"].value("conda_dir", "");
  parameters.number_of_threads = json_parameters["Spatial_Distancing"].value("number_of_threads", 1);
  std::mt19937 gen(std::random_device{}());
  std::uniform_int_distribution dis;
  parameters.system_seed = json_parameters["Spatial_Distancing"].value("seed", dis(gen));
  for (const auto &item: json_parameters["Spatial_Distancing"]["parameter_screening"].items()) {
      std::vector<std::string> values;
      for (float value: item.value()) {
          std::string valuestr = std::to_string(value);
          if (std::stoi(valuestr.substr(valuestr.find_first_of('.') + 1, -1)) == 0) {
              values.push_back(std::to_string(static_cast<int>(value)));
          } else {
              values.push_back(std::to_string(value));
          }
      }
      parameters.screening_parameters[item.key()] = values;
  }
  parameters.screen_start_idx = json_parameters["Spatial_Distancing"].value("screen_start_idx", 1);
  json_file.close();
  return parameters;
}


OutputParameters getOutputParameters(const std::string &output_config) {
  OutputParameters parameters{};
  std::ifstream json_file(output_config);
  json json_parameters;
  json_file >> json_parameters;
  parameters.activated = json_parameters["Spatial_Distancing"].value("activated", false);
  parameters.output_graphs = json_parameters["Spatial_Distancing"].value("output_graphs", false);
  parameters.output_ovito = json_parameters["Spatial_Distancing"].value("output_ovito", false);
  json_file.close();
  return parameters;
}
SimulationParameters getSimulationParameters(const std::string &simulator_config) {
  SimulationParameters parameters{};
  std::ifstream json_file(simulator_config);
  json json_parameters;
  json_file >> json_parameters;
  try {
    parameters.topic = json_parameters["Spatial_Distancing"]["topic"];
    parameters.dimensions = json_parameters["Spatial_Distancing"]["dimensions"];
    parameters.max_time = json_parameters["Spatial_Distancing"]["max_time"];
    parameters.time_stepping = json_parameters["Spatial_Distancing"]["timestepping"];
    parameters.dt_auto = json_parameters["Spatial_Distancing"].value("dt_auto", false);
  } catch (const std::bad_optional_access &e) {
    ERROR_STDERR(e.what());
    throw;
  }
  for (const auto &site: json_parameters["Spatial_Distancing"]["Sites"]) {
    std::unique_ptr<SimulationParameters::SiteParameters> site_para;
    const auto molecular_layer = site.value("molecular_layer", false);

    const auto type = site["type"];
    if ("CuboidSite" == type) {
      auto cs_para = std::make_unique<SimulationParameters::CuboidSiteParameters>();
      cs_para->upper_bound = {site["CuboidSite"]["x_range"][0], site["CuboidSite"]["y_range"][0], site["CuboidSite"]["z_range"][0]};
      cs_para->lower_bound = {site["CuboidSite"]["x_range"][1], site["CuboidSite"]["y_range"][1], site["CuboidSite"]["z_range"][1]};
      if (molecular_layer) {
          cs_para->molecules_grid_size = {site["CuboidSite"]["molecules_grid_size"][0], site["CuboidSite"]["molecules_grid_size"][1], site["CuboidSite"]["molecules_grid_size"][2]};
      }
      site_para = std::move(cs_para);
    }
    site_para->type = type;
    site_para->molecular_layer_ = molecular_layer;
    site_para->identifier = site["identifier"];
    if(auto tmp = site.find("stopping_criteria"); tmp!= site.end()){
        SimulationParameters::StoppingCriteria stop_param;
        stop_param.molecule = site["stopping_criteria"].value("molecule", "");
        stop_param.threshold = site["stopping_criteria"].value("threshold", 0.0);
        stop_param.time = site["stopping_criteria"].value("time", parameters.max_time);
        stop_param.init_time = stop_param.time;
        site_para->stopping_criteria = stop_param;
    }

    //load Molecule manager
      if (auto mol = site.find("MoleculeManager"); mol != site.end()) {
          std::shared_ptr<SimulationParameters::MoleculeManagerParameters> molecule_manager_parameters = std::make_unique<SimulationParameters::MoleculeManagerParameters>();
          site_para->molecule_manager_parameters.boundary_condition = site["MoleculeManager"]["boundary_condition"];
          for (const auto& molecule_type : site["MoleculeManager"]["Types"]) {
              auto molecule = site["MoleculeManager"]["Molecules"].at(std::string(molecule_type));
              std::shared_ptr<SimulationParameters::MoleculeParameters> molecule_parameters = std::make_unique<SimulationParameters::MoleculeParameters>();
              molecule_parameters->name = molecule_type;
              molecule_parameters->diffusion_coefficient = molecule.value("DiffusionCoefficient", 0.0);
              molecule_parameters->decay = molecule.value("Decay", 0.0);
              molecule_parameters->initial_concentration = molecule.value("InitialConcentration", 0.0);
              if(const auto& tmp = molecule.find("Boundaryflow"); tmp != molecule.end()){
                  molecule_parameters->boundary_flow.amplitude = molecule["Boundaryflow"].value("amplitude", 0.0);
                  molecule_parameters->boundary_flow.frequence= molecule["Boundaryflow"].value("frequence", 0.0);
                  molecule_parameters->boundary_flow.target_conc= molecule["Boundaryflow"].value("target_conc", 0.0);
              }
              else{
                  molecule_parameters->boundary_flow.amplitude = 0.0;
                  molecule_parameters->boundary_flow.frequence = 0.0;
              }
              molecule_parameters->spontaneous_rate = molecule.value("SpontaneousRate", 0.0);

              site_para->molecule_manager_parameters.molecules.emplace_back(std::move(molecule_parameters));
          }
          for (const auto& interaction :site["MoleculeManager"]["MoleculeMoleculeInteractions"].items()){
              SimulationParameters::MoleculeMoleculeInteractionParameters molecule_molecule_parameters;
              molecule_molecule_parameters.first = interaction.value().value("first", "");
              molecule_molecule_parameters.second = interaction.value().value("second", "");
              molecule_molecule_parameters.binding = interaction.value().value("binding", 0.0);
              molecule_molecule_parameters.unbinding = interaction.value().value("unbinding", 0.0);
              site_para->molecule_manager_parameters.molecule_molecule_interactions.emplace_back(std::make_pair(interaction.key(), std::move(molecule_molecule_parameters)));
          }
      }

    // load agent manager
    // we need to keep an ordering for reproducing simulations, changing it leads to agents get differently initialized due to random values
    for (const auto& agent_type:site["AgentManager"]["Types"]) {
      std::shared_ptr<SimulationParameters::AgentParameters> agent_parameters{};
      auto agent = site["AgentManager"]["Agents"].at(std::string(agent_type));
      agent_parameters = std::make_shared<SimulationParameters::AgentParameters>();
      agent_parameters->track_inside_conc = agent.value("inside_conc", false);
      agent_parameters->initial_distribution = agent.value("initial_distribution", 0);
      agent_parameters->pos_dist = agent.value("pos_dist", -1.0);
      agent_parameters->morphology_parameters.type = "SphericalMorphology";
      agent_parameters->morphology_parameters.radius = agent["Morphology"]["SphericalMorphology"]["radius"];

      agent_parameters->type = agent_type;
      agent_parameters->number = agent.value("number", 0);

      for (const auto &molecule:agent["Molecule Interactions"].items()) {
            SimulationParameters::MoleculeInteractionParameters molecule_interaction_parameters{};
            molecule_interaction_parameters.secretion  = molecule.value().value("SecretionRate",0.0);
            molecule_interaction_parameters.secretion_delay = molecule.value().value("Secretion_delay", 0.0);
            molecule_interaction_parameters.uptake = molecule.value().value("UptakeRate", 0.0);
            molecule_interaction_parameters.inside_conc_decay = molecule.value().value("InsideConcDecay", 0.0);
            molecule_interaction_parameters.secretion_variable_over_time = molecule.value().value("secretion_variable_over_time", false);
            agent_parameters->molecule_interactions.emplace_back(std::make_pair(molecule.key(),std::move(molecule_interaction_parameters)));
      }
      site_para->agent_manager_parameters.agents.emplace_back(std::move(agent_parameters));
    }



    parameters.site_parameters = std::move(site_para);
  }
  json_file.close();
  return parameters;
}

AnalyserParameters getAnalyserParameters(const std::string &analyser_config) {
  AnalyserParameters parameters{};
  std::ifstream json_file(analyser_config);
  json json_parameters;
  json_file >> json_parameters;
  parameters.active_measurements = json_parameters["Spatial_Distancing"]["active_measurements"].get<std::unordered_set<std::string>>();
  parameters.cell_state_count = json_parameters["Spatial_Distancing"]["cell_state_count"].get<std::vector<std::string>>();
  parameters.output_every_X_time_unit = json_parameters["Spatial_Distancing"].value("output_every_X_time_unit", 0.0);
  json_file.close();
  return parameters;
}


InputParameters getInputParameters(const std::string &input_config) {
  InputParameters parameters{};
  std::ifstream json_file(input_config);
  json json_parameters;
  json_file >> json_parameters;

  for (auto &rate:json_parameters["Spatial_Distancing"]["Rates"].items()) {
    std::unique_ptr<InputParameters::DefaultRateParameters> para_rate;
    const auto type = rate.value().value("type", "ConstantRate");
   if (type == "MoleculeRate") {
        auto molecule_rate = std::make_unique<InputParameters::MoleculeRateParameters>();
        molecule_rate->molecule_name = rate.value()["Molecule Name"];
        para_rate = std::move(molecule_rate);

    }else {
      para_rate = std::make_unique<InputParameters::DefaultRateParameters>();
    }
    para_rate->type = type;
    para_rate->key = rate.key();
    para_rate->rate = rate.value().value("rate", 0.0);
    para_rate->condition = rate.value().value("condition", "");
    parameters.rates.emplace_back(std::move(para_rate));
  }
  json_file.close();
  return parameters;
}

void executeShellCommand(const std::string &command, const bool suppress_output) {
  if (auto system_return = std::system(command.c_str());!suppress_output) {
    DEBUG_STDOUT("shell.command: " + command + "with return " + std::to_string(system_return));
  }
}

}