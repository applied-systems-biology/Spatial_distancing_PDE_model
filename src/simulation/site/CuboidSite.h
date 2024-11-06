//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef _CUBOIDSITE_H
#define	_CUBOIDSITE_H

#include <string>
#include <iostream>

#include "simulation/Site.h"


class CuboidSite : public Site {
 public:
    // Class for a cuboid site
    CuboidSite(Randomizer *random_generator,
               std::shared_ptr<InSituMeasurements> measurements,
               std::string config_path, std::unordered_map<std::string,
            std::string> cmd_input_args, std::string output_path);

  [[nodiscard]] bool containsPosition(Coordinate3D) final;
  [[nodiscard]] std::string getType() const override{return "CuboidSite";}
  Coordinate3D getRandomPosition(double radius) final;
  Coordinate3D getCenterPosition() final;
  Coordinate3D getCloseToCenterPosition() final;
  Coordinate3D getZEqualsZeroPosition(double radius) final;
  Coordinate3D getZandYequalsZeroPosition(double radius) final;
  Coordinate3D getLowerLimits() final;
  Coordinate3D getUpperLimits() final;

  std::vector<Coordinate3D> getSystemBoundaries() {return {lower_bound_, upper_bound_};};

  std::pair<std::vector<std::tuple<int, int, int>>, std::vector<std::tuple<int, int, int>>> cell_grid_info(double radius, Coordinate3D position) final;

  std::vector<std::tuple<int, int, int>> covered_points(double radius, Coordinate3D position) final;
  std::vector<std::tuple<int, int, int>> surface_points(double radius, Coordinate3D position) final;

    void initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                                      double current_time,
                                      double time_delta);
 protected:

  float average_conc_in_agent(const std::string& name, double radius, Coordinate3D position) final;


  Coordinate3D upper_bound_;
  Coordinate3D lower_bound_;
  std::tuple<int, int, int> molecules_grid_size_;

    bool check_collision(Coordinate3D pos_to_check);

    Coordinate3D two_cells_layout(double radius, int cell, double pos_dist);
};

#endif	/* _CUBOIDSITE_H */