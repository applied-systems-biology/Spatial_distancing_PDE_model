//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef CELL_H
#define    CELL_H
#include <memory>

#include "simulation/Agent.h"
#include "utils/io_util.h"
#include "simulation/Site.h"

class AgentProperties;
class Analyser;

class Cell : public Agent {
 public:
    /// Class for an abstract cell that is further specified by inherited classes
  Cell(std::unique_ptr<Coordinate3D>, int, Site *, double time_delta, double current_time);

  void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters);

    /*!
     * Central function: Performs all actions for one timestep for one cell
     * @param timestep Double for courrent timestep
     * @param current_time Double for current time
     */
  void doAllActionsForTimestep(double timestep, double current_time) final;

    /*!
     * Central function: Performs all actions for one timestep for one cell related to molecular interactions
     * @param timestep Double for courrent timestep
     * @param current_time Double for current time
     */
  void interactWithMolecules(double timestep, double current_time);

  std::string getAgentCSVTagAoc(double current_time) final;
  std::string getAgentCSVTagToc() final;

  std::string getTypeName() override;

  Morphology *getSurface();

  std::map<std::string, double> get_molecule_uptake() {return uptake_rate;}

  std::map<std::string, double> get_inside_conc() const {return inside_conc;}
  void update_inside_conc(std::string molecule, double value);

  void update_secretion_rate(const std::string& secreted_molecule, const std::string& uptaken_molecule, const double& value, double dt);

  void decay(std::string molecule, double rate, double timestep);

  bool get_track_inside_conc() {return track_inside_conc;};

  std::vector<std::tuple<int, int, int>> get_membrane() const {return membrane_grid_points;}
  virtual   std::vector<std::tuple<int, int, int>> get_inside() const {return inside_grid_points;}

  double get_current_uptake(std::string name){return current_uptake_.at(name);}

protected:

  std::shared_ptr<Morphology> surface;

  //Molecule Parameters
  std::map<std::string, double> secretion_rate;
  std::map<std::string, double> initial_rate_;
  std::map<std::string, bool> secretion_variable_over_time;
  std::map<std::string, double> secretion_delay;
  std::map<std::string, double> uptake_rate;
  std::map<std::string, double> inside_conc_decay;
  std::map<std::string , double> inside_conc;
  bool secretion_start;
  bool track_inside_conc;
  std::map<std::string, double> current_uptake_;
  std::map<std::string, std::vector<double>> previous_uptake_;

  std::vector<std::tuple<int, int, int>> membrane_grid_points;
  std::vector<std::tuple<int, int, int>> inside_grid_points;
};

#endif    /* CELL_H */

