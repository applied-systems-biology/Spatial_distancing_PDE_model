//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "simulation/Site.h"

#include <numeric>
#include <algorithm>
#include "external/json.hpp"

#include "simulation/Algorithms.h"

#include "simulation/AgentManager.h"
#include "utils/macros.h"

using json = nlohmann::json;

Site::Site(Randomizer *random_generator,
           std::shared_ptr<InSituMeasurements> measurements,
           std::string config_path, std::unordered_map<std::string, std::string> cmd_input_args,
           std::string output_path) : random_generator_(
        random_generator), measurements_(std::move(measurements)) {

    measurements_->setSite(this);
}

void Site::doAgentDynamics(Randomizer *random_generator, const SimulationTime &time) {
    measurements_->observeMeasurements(time);

    const auto current_time = time.getCurrentTime();
    const auto dt = time.getCurrentDeltaT();

    const auto &all_agents = agent_manager_->getAllAgents();
    if (!all_agents.empty()) {
        // Generate a random order for executing agent and molecule tasks
        std::vector<unsigned int> task_order = Algorithms::generateRandomPermutation(random_generator, 2);

        for (const auto &task : task_order) {
            if (task == 0) {  // Molecule task
                if (molecular_layer_) {
                    double radius = 0.0;
                    std::vector<std::vector<std::tuple<int, int, int>>> covered_point, surface_point;

                    // Gather surface and covered points for all agents
                    for (const auto &ag : all_agents) {
                        covered_point.emplace_back(ag->get_inside());
                        surface_point.emplace_back(ag->get_membrane());
                        radius = ag->getAgentProperties()->getMorphology()->getBasicSphereOfThis()->getRadius();
                    }

                    // Flatten vectors
                    auto cov = std::accumulate(covered_point.begin(), covered_point.end(), decltype(covered_point)::value_type{},
                                               [](auto &x, auto &y) {
                                                   x.insert(x.end(), y.begin(), y.end());
                                                   return x;
                                               });
                    auto surf = std::accumulate(surface_point.begin(), surface_point.end(), decltype(surface_point)::value_type{},
                                                [](auto &x, auto &y) {
                                                    x.insert(x.end(), y.begin(), y.end());
                                                    return x;
                                                });

                    // Randomize order of diffusion and reaction tasks
                    std::vector<unsigned int> molecule_order = Algorithms::generateRandomPermutation(random_generator, 2);
                    for (const auto &molecule_task : molecule_order) {
                        switch (molecule_task) {
                        case 0:
                            molecule_manager_->doDiffusion(dt, current_time, cov, surf, radius);
                            break;
                        case 1:
                            molecule_manager_->do_reactions(dt, radius);
                            break;
                        }
                    }
                }
            } else {  // Agent task
                const auto current_permutation = Algorithms::generateRandomPermutation(random_generator, all_agents.size());
                for (auto agent_idx = current_permutation.begin(); agent_idx < current_permutation.end(); ++agent_idx) {
                    auto curr_agent = all_agents[*agent_idx];
                    if (nullptr != curr_agent) {
                        curr_agent->doAllActionsForTimestep(dt, current_time);
                        if (curr_agent->isDeleted()) {
                            for (const auto &sphere : curr_agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis()) {
                                agent_manager_->removeSphereRepresentation(sphere);
                            }
                            curr_agent = nullptr;
                        }
                    }
                }
                agent_manager_->cleanUpAgents();
            }
        }
    }

    // Check if steady state reached
    if (!steady_state_) {
        steady_state_ = molecule_manager_->steady_state_reached(stopping_criteria_, dt);
    }
}

void Site::handleCmdInputArgs(std::unordered_map<std::string, std::string>cmd_input_args) {
    // Parameters to be screened or from cmd input
    // You can add your parameters you want to screen below here
    if (cmd_input_args.size() > 0) {
        for (const auto &[key, value]: cmd_input_args) {
            if (key == "binding" || key=="unbinding") {
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecule_molecule_interactions) {
                    if (mol.first == "Complex") {
                        if ("binding" == key) {
                            mol.second.binding = std::stod(value);
                            SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                        } else if ("unbinding" == key) {
                            mol.second.unbinding = std::stod(value);
                            SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                        }
                    }
                }
            }
            else if (key=="uptake" || key == "secretion" || key == "secretion_delay"){
                for (auto& agent: parameters_.site_parameters->agent_manager_parameters.agents){
                    for (auto& mol_interact: agent->molecule_interactions) {
                        if ("uptake" == key) {
                            if (mol_interact.first == "AMP") {
                                mol_interact.second.uptake = std::stod(value);
                                SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                            }
                        }
                        if("secretion" == key){
                            if (mol_interact.first == "Defensive") {
                                //mol_interact.second.secretion_variable_over_time = false;
                                mol_interact.second.secretion = std::stod(value);
                                SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                            }
                        }
                        if (key== "secretion_delay"){
                            if (mol_interact.first == "Defensive"){
                                mol_interact.second.secretion_delay = std::stod(value);
                                SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                            }
                        }
                    }
                }
            }
            else if (key == "degrad_AMP" || key == "degrad_DEF" || key == "D_AMP" || key == "D_DEF"){
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecules){
                    if (mol->name == "AMP" && key == "degrad_AMP"){
                        mol->decay = std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                    if (mol->name == "Defensive" && key == "degrad_DEF"){
                        mol->decay = std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                    if (mol->name == "AMP" && key == "D_AMP"){
                        mol->diffusion_coefficient = std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                    if (mol->name == "Defensive" && key == "D_DEF"){
                        mol->diffusion_coefficient = std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                }
            }
            else if (key == "D_COM"){
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecules){
                    if (mol->name == "Complex" && key == "D_COM"){
                        mol->diffusion_coefficient = std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                }
            }
            else if (key == "D"){
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecules){
                    if (mol->name == "Complex"){
                        mol->diffusion_coefficient = std::stod(value)/2;
                    }
                    else{
                        mol->diffusion_coefficient = std::stod(value);
                    }
                }
                SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
            }
            else if(key == "Kd"){
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecule_molecule_interactions) {
                    if (mol.first == "Complex") {
                        mol.second.binding = std::stod(value);
                        mol.second.unbinding = std::stod(value)*44;
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                }
            }
            else if(key == "d"){
                for (auto &agent: parameters_.site_parameters->agent_manager_parameters.agents){
                    if (agent->initial_distribution == 5){
                        agent->pos_dist= std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                }
            }
            else if (key == "flow"){
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecules) {
                    if (mol->name == "AMP") {
                        mol->boundary_flow.amplitude = std::stod(value);
                        SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                    }
                }
            }
            else if (key == "degrad"){
                for (auto &mol: parameters_.site_parameters->molecule_manager_parameters.molecules){
                    mol->decay = std::stod(value);
                    SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                }
            }
        }
    }
}