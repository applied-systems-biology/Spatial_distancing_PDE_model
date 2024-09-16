//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef _SITE_H
#define _SITE_H

#include <memory>
#include <utility>
#include <vector>
#include <algorithm>

#include "utils/time_util.h"

#include "basic/Coordinate3D.h"

#include "simulation/Cell.h"
#include "simulation/Algorithms.h"

#include "analyser/InSituMeasurements.h"
#include "simulation/morphology/SphereRepresentation.h"

#include "simulation/AgentManager.h"

#include "MoleculeManager.h"
#include "io/output_handler.h"

#include "CellFactory.h"

class Agent; //forward declaration
class Analyser;
class Randomizer;
class CellFactory;


class Site {

 public:
    /// Class for handling environment interactions and wrapping functionality of all aspects that are happening during the simulation inside of the environment (i.e. site)

    Site(Randomizer *random_generator,
       std::shared_ptr<InSituMeasurements> measurements,
       std::string config_path,
       std::unordered_map<std::string, std::string> cmd_input_args,
       std::string output_path);

  virtual ~Site() = default;

  friend void OutputHandler::outputCurrentConfiguration(const Site &site,
                                                        const SimulationTime &time,
                                                        int run,
                                                        int seed,
                                                        bool simulation_end,
                                                        std::unordered_map<std::string, std::string> cmd_input_args) const;

  friend void InSituMeasurements::observeMeasurements(const SimulationTime &time);

  virtual bool containsPosition(Coordinate3D position) = 0;

  [[nodiscard]] virtual std::string getType() const = 0;

  virtual Coordinate3D getRandomPosition(double radius) = 0;
  virtual Coordinate3D getCenterPosition() = 0;
  virtual Coordinate3D getCloseToCenterPosition() = 0;
  virtual Coordinate3D getZEqualsZeroPosition(double radius) = 0;
  virtual Coordinate3D getZandYequalsZeroPosition(double radius) = 0;
  virtual Coordinate3D getLowerLimits() = 0;
  virtual Coordinate3D getUpperLimits() = 0;

  virtual std::vector<Coordinate3D> getSystemBoundaries() {return {}; };

  void doAgentDynamics(Randomizer *random_generator, const SimulationTime &time);

  Randomizer *getRandomGenerator() { return random_generator_; }
  AgentManager *getAgentManager() const { return agent_manager_.get(); }
  InSituMeasurements *getMeasurments() const { return measurements_.get(); }

  [[nodiscard]] unsigned int getNumberOfSpatialDimensions() const { return dimensions; }

  [[nodiscard]] std::string getIdentifier() const { return identifier_; }
  virtual   std::pair<std::vector<std::tuple<int, int, int>>, std::vector<std::tuple<int, int, int>>> cell_grid_info(double radius, Coordinate3D position) = 0;

  virtual std::vector<std::tuple<int, int, int>> covered_points(double radius, Coordinate3D position) = 0;
  virtual std::vector<std::tuple<int, int, int>> surface_points(double radius, Coordinate3D position) = 0;


  virtual float average_conc_in_agent(const std::string& name, double radius, Coordinate3D position) = 0;

  std::shared_ptr<MoleculeManager> molecule_manager_;

  bool get_molecular_layer() const {return molecular_layer_;};


  void handleCmdInputArgs(std::unordered_map<std::string, std::string> cmd_input_args);

    double getMaxTime() {return parameters_.max_time;}
    double getTimeStepping() {return parameters_.time_stepping;}
    CellFactory *getCellFactory() const { return cell_factory_.get(); }

    bool get_steady_state() {return steady_state_;}


private:
protected:
  void setBoundaryCondition();


  unsigned int dimensions{};
  std::string identifier_{};


  Randomizer *random_generator_;
  std::shared_ptr<InSituMeasurements> measurements_;
  std::unique_ptr<AgentManager> agent_manager_;

    bool steady_state_{};

    bool molecular_layer_{};

    std::string config_path_{};

    abm::util::SimulationParameters parameters_{};

    std::unique_ptr<CellFactory> cell_factory_;

    abm::util::SimulationParameters::StoppingCriteria stopping_criteria_;

};


#endif /* _SITE_H */

