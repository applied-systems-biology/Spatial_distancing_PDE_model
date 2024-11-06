//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "simulation/Cell.h"

#include "simulation/Agent.h"
#include "simulation/AgentProperties.h"
#include "analyser/Analyser.h"

#include "simulation/morphology/SphericalMorphology.h"

#include "analyser/InSituMeasurements.h"
#include "utils/macros.h"

#include "simulation/site/CuboidSite.h"

Cell::Cell(std::unique_ptr<Coordinate3D> c,
           int id,
           Site *site,
           double time_delta,
           double current_time)
    : Agent(std::move(c), id, site) {
}

void Cell::doAllActionsForTimestep(double timestep, double current_time) {
    std::vector<unsigned int>::iterator it;
    interactWithMolecules(timestep, current_time);
}

std::string Cell::getTypeName() {
  return "Cell";
}

void Cell::interactWithMolecules(double timestep, double current_time) {
    bool mol_layer = site->get_molecular_layer();

    if (mol_layer){
        for (const auto sphere : agentProps->getMorphology()->getAllSpheresOfThis()) {

            // SECRETION
            for (const auto& [name, rate] : secretion_rate) {
                if (secretion_start) {
                    site->molecule_manager_->secreteMolecules(name, rate, membrane_grid_points, timestep, current_time, this->getSurface()->getBasicSphereOfThis()->getRadius());
                }
            }

            // UPTAKE
            for (const auto& [name, rate] : uptake_rate) {
                current_uptake_[name] = site->molecule_manager_->uptakeMolecules(name, rate, membrane_grid_points, timestep);
                if (this->track_inside_conc){
                    this->update_inside_conc(name, current_uptake_[name]);
                }
                for (auto& [mol_name, tmp] : secretion_variable_over_time){
                    if (secretion_variable_over_time[mol_name]) {
                        if (current_time >= secretion_delay[mol_name]){
                            secretion_start= true;
                            previous_uptake_[mol_name].push_back(current_uptake_[name]);
                            previous_uptake_[mol_name].erase(previous_uptake_[mol_name].begin());
                        }
                        else{
                            previous_uptake_[mol_name].push_back(current_uptake_[name]);
                        }
                        this->update_secretion_rate(mol_name, name, previous_uptake_[mol_name].front(), timestep);
                    }
                }
            }
            // INSIDE CONCENTRATION UPDATE
            for (const auto& [name, rate] : inside_conc_decay) {
                if (this->track_inside_conc) {
                    this->decay(name, rate, timestep);
                }
            }
        }
    }
}

Morphology *Cell::getSurface() {
  return surface.get();
}

void Cell::setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) {

    // Initialize the variable with a negative value, so that in the first run the
    // cell is not marked as checked
    surface = std::make_shared<Morphology>(this);
    if (parameters->morphology_parameters.type == "SphericalMorphology") {
    auto radius = site->getRandomGenerator()->generateNormalDistributedValue(parameters->morphology_parameters.radius,
                                                                             0.0);

    surface->appendAssociatedCellpart(std::make_unique<SphericalMorphology>(surface.get(),
                                                                            position,
                                                                            radius,
                                                                            0.0,
                                                                            "basic"));
    }
    agentProps->setMorphology(surface);
    secretion_start = true;

    for(auto &[name, data]:parameters->molecule_interactions) {
        if (data.secretion != 0) {
            secretion_rate.emplace(std::pair<std::string, double>(name, data.secretion));
            initial_rate_.emplace(std::pair<std::string, double>(name, data.secretion));
        }
        if (data.secretion_delay != 0){
            secretion_delay.emplace(std::pair<std::string, double>(name, data.secretion_delay));
            if (data.secretion_delay){
                secretion_start = false;
            }
        }
        if (data.secretion_variable_over_time != 0){
            secretion_variable_over_time.emplace(std::pair<std::string, bool>(name, data.secretion_variable_over_time));
        }
        if (data.uptake != 0) {
            uptake_rate.emplace(std::pair<std::string, double>(name, data.uptake));
        }
        if(data.inside_conc_decay != 0){
            inside_conc_decay.emplace(std::pair<std::string, double>(name, data.inside_conc_decay));
        }
    }
    if (site->molecule_manager_ != nullptr){
        for (const auto &name: site->molecule_manager_->get_molecules_names()){
          if (parameters->track_inside_conc) {
              track_inside_conc = parameters->track_inside_conc;
              inside_conc[name] = 0.0;
        }
          current_uptake_[name] = 0.0;
          previous_uptake_[name].push_back(0.0);
      }
    }

    auto cell_info = site->cell_grid_info(agentProps->getMorphology()->getBasicSphereOfThis()->getRadius(), getCurrentPosition());
    membrane_grid_points = std::get<1>(cell_info);
    inside_grid_points = std::get<0>(cell_info);
    //std::cout << "Nb of inside points = " << inside_grid_points.size();

}

void Cell::update_inside_conc(std::string molecule, double value) {
    inside_conc[molecule] += value;
}

void Cell::update_secretion_rate(const std::string& secreted_molecule, const std::string& uptaken_molecule, const double& value, double dt){
    secretion_rate[secreted_molecule] = value * initial_rate_[secreted_molecule];
    if (secretion_rate[secreted_molecule] > 100000){ //max secretion rate
        secretion_rate[secreted_molecule] = 100000;
    }
    secretion_rate[secreted_molecule] = secretion_rate[secreted_molecule]/dt;
}

void Cell::decay(std::string molecule, double rate, double timestep) {
    inside_conc[molecule] -= timestep * rate*inside_conc[molecule];
}


std::string Cell::getAgentCSVTagAoc(double current_time) {
    std::ostringstream csvTag;
    csvTag << current_time << "\t" <<
           position->x << "\t" <<
           position->y << "\t" <<
           position->z << "\t" <<
           getTypeName();
    return csvTag.str();
}

std::string Cell::getAgentCSVTagToc() {

    std::ostringstream csvTag;
    csvTag << id << "\t" <<
           position->x << "\t" <<
           position->y << "\t" <<
           position->z << "\t" <<
           getTypeName();

    return csvTag.str();
}
