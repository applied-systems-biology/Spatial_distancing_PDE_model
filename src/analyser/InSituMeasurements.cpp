//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "analyser/InSituMeasurements.h"

#include "simulation/Site.h"
#include "simulation/Cell.h"

#include "utils/macros.h"
#include "cmath"

InSituMeasurements::InSituMeasurements(std::unordered_set<std::string> active_measurements, const std::string &id, double output_every_X_time_unit)
    : active_measurements_(std::move(active_measurements)) {
    every_X_time_unit_ = output_every_X_time_unit;
    for (const auto &active_:active_measurements_) {
        // Remove everything after "%" in the active measurement string
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
          active = active_.substr(0, char_pos);
        }
        if ("agent-statistics" == active) {
          pair_measurements_["agent-statistics"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                              "radius", "x", "y", "z",
                                                                                              "cellpart", "cellpart_id");
        } else if ("environment" == active) {
            pair_measurements_["environment"] = std::make_unique<PairMeasurement>(id, "x", "y", "z", "radius_or_length", "type", "additional");
        } else if ("record-random-numbers" == active) {
          pair_measurements_["record-random-numbers"] = std::make_unique<PairMeasurement>(id, "origin", "value");
        }
        else if ("molecules" == active) {
            pair_measurements_["molecules"] = std::make_unique<PairMeasurement>(id, "name", "time", "x", "y", "z", "value");
        }
        else if ("molecules_test" == active){
            // auto mol_names = site_->molecule_manager_->get_molecules_names();
            pair_measurements_["molecules_test"] = std::make_unique<PairMeasurement>(id, "time", "x", "y", "z", "AMP", "Complex", "Defensive");
        }
        else if ("mol" == active){
            pair_measurements_["mol"] = std::make_unique<PairMeasurement>(id, "name", "time", "grid_point_id", "value");
        }
        else if("mol_environment" == active){
          pair_measurements_["mol_environment"]= std::make_unique<PairMeasurement>(id, "grid_point_id", "x", "y", "z");
        }
        else if ("inside_conc" == active){
            pair_measurements_["inside_conc"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid", "mol_type", "inside_concentration");
            //histogram_measurements_["conc_kill"] = std::make_unique<HistogramMeasurement>(id, "concentration");
        }
        else if ("membrane" == active){
            pair_measurements_["membrane"] = std::make_unique<PairMeasurement>(id, "agentID", "TypeName", "x", "y", "z", "distance");
        }
        else if ("inside" == active){
            pair_measurements_["inside"] = std::make_unique<PairMeasurement>(id, "agentID", "TypeName", "x", "y", "z");
        }
        else if("spatial_dist_center" == active){
            pair_measurements_["spatial_dist_center"] = std::make_unique<PairMeasurement>(id, "time", "dist_to_cell", "AMP", "Complex", "Defensive");
        }
        else if("spatial_dist_center_updated" == active){
            pair_measurements_["spatial_dist_center_updated"] = std::make_unique<PairMeasurement>(id, "time", "dist_to_cell", "AMP", "Complex", "Defensive", "AMP_uptaken");
        }
        else if("spatial_dist_cell" == active){
            pair_measurements_["spatial_dist_cell"] = std::make_unique<PairMeasurement>(id, "time", "dist_to_cell", "AMP", "Complex", "Defensive");
        }
        else if("final_conc" == active){
            pair_measurements_["final_conc"] = std::make_unique<PairMeasurement>(id, "AMP", "Complex", "Defensive");
        }
        else if("dist_2_cells" == active){
            pair_measurements_["dist_2_cells"] = std::make_unique<PairMeasurement>(id, "time", "x", "AMP", "Complex", "Defensive", "center_cell", "radius");
        }
        else if ("mod_distrib" == active){
            pair_measurements_["mod_distrib"] = std::make_unique<PairMeasurement>(id, "time", "AMP", "Complex", "Defensive");
        }
        else if ("uptake" == active){
            pair_measurements_["uptake"] = std::make_unique<PairMeasurement>(id, "time", "molecule", "uptake");
        }
    }
}

void InSituMeasurements::observeMeasurements(const SimulationTime &time) {
  using std::make_pair;
  const auto &current_time = time.getCurrentTime();

  for (const auto &active_:active_measurements_) {

    // Code snippet to let your measurement only execute ever X-th simulation step
    // Name your measurement as "measurment%X" -> example: "agent-statistics%25" ... every 25th step
    // Default if no value is given: Every 10th step
    int every_x_step = 10;
    auto active = active_;
    int char_pos = active_.find("%");
    if (char_pos > 0) {
      every_x_step = stoi(active_.substr(char_pos+1, active_.size()));
      active = active_.substr(0, char_pos);
    }
    bool do_measurement = false;

    if (every_X_time_unit_ != 0) {
        int xstep = round(every_X_time_unit_ / time.getCurrentDeltaT()); // for every X time unit output
        do_measurement = time.getCurrentTimeStep() % xstep == 0;
    }
    else{
        do_measurement = (time.getCurrentTimeStep() % every_x_step) == 0;
    }

    // Measurements ever X timestep
    if ("molecules" == active && do_measurement) {
      double value;
      auto const location = site_->molecule_manager_->get_location();
      const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
      for (const auto& [name, data] : *site_->molecule_manager_) {
          for (int i= 0; i< n_x*n_y*n_z;++i){
              value = data[i];
              auto const loc = location[i];
              pair_measurements_["molecules"]->addValuePairs(name, current_time, loc.x, loc.y, loc.z, value);
          }
      }
    } else if ("molecules_test" == active && do_measurement){
        auto const location = site_->molecule_manager_->get_location();
        const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
        for (int i= 0; i< n_x*n_y*n_z;++i){
            std::vector<double> value;
            for (const auto& [name, data] : *site_->molecule_manager_) {
                value.push_back(data[i]);
            }
            auto const loc = location[i];
            //SYSTEM_STDOUT("time: " << current_time << "| location: " << loc.x << ", " << loc.y << ", " << loc.z << "| values: " << value[0] << ", " << value[1] << ", " << value[2]);
            pair_measurements_["molecules_test"]->addValuePairs(current_time, loc.x, loc.y, loc.z, value[0], value[1], value[2]);
        }
    }else if("spatial_dist_center" == active && do_measurement) {
        const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
        auto dist = 0.0;
        for (int i = round(n_z / 2); i < n_z; ++i) {
            std::vector<double> value;
            for (const auto& [name, data] : *site_->molecule_manager_) {
                value.push_back(data[n_x / 2 + n_x * n_y / 2 + n_y * n_x * i]); // center of the simulation site
                Coordinate3D position = {0, 0, 0};
                dist = position.calculateEuclidianDistance(site_->molecule_manager_->get_loc(i, n_y/2, n_x/2));
            }
            pair_measurements_["spatial_dist_center"]->addValuePairs(time.getCurrentTime(), dist, value[0], value[1], value[2]);
        }
    }else if("spatial_dist_center_updated" == active && do_measurement) {
        const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
        auto dist = 0.0;
        for (int i = round(n_z / 2); i < n_z; ++i) {
            std::vector<double> value;
            for (const auto& [name, data] : *site_->molecule_manager_) {
                value.push_back(data[n_x / 2 + n_x * n_y / 2 + n_y * n_x * i]); // center of the simulation site
                Coordinate3D position = {0, 0, 0};
                dist = position.calculateEuclidianDistance(site_->molecule_manager_->get_loc(i, n_y/2, n_x/2));
            }
            auto inside = 0.0;
            if (dist < 3.4){
                for (const auto ag : site_->getAgentManager()->getAllAgents()) {
                    if (ag->get_track_inside_conc()) {
                        inside = ag->get_inside_conc().at("AMP");
                    }
                }
            }
            pair_measurements_["spatial_dist_center_updated"]->addValuePairs(time.getCurrentTime(), dist, value[0], value[1], value[2], inside);
        }
    }else if ("mol" == active && do_measurement){
        double value;
        const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
        for (auto& [name, data] : *site_->molecule_manager_) {
            for (int i = 0; i < n_x * n_y * n_z; ++i) {
                value = data[i];
                pair_measurements_["mol"]->addValuePairs(name, current_time, i, value);
            }
        }
    }else if ("agent-statistics" == active && do_measurement) {
      for (const auto agent: site_->getAgentManager()->getAllAgents()) {
        std::string cellpart = "Mothercell";
        for (const auto cellparts: agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis()) {
          pair_measurements_["agent-statistics"]->addValuePairs(current_time,
                                                                agent->getTypeName(), agent->getId(),
                                                                cellparts->getRadius(),
                                                                cellparts->getPosition().x,
                                                                cellparts->getPosition().y,
                                                                cellparts->getPosition().z,
                                                                cellpart, cellparts->getDescription());
        }
      }
    }
    else if ("inside_conc" == active && do_measurement) {
        for (const auto ag : site_->getAgentManager()->getAllAgents()) {
            if (ag->get_track_inside_conc()){
                for (const auto& [name, data] : *site_->molecule_manager_) {
                    pair_measurements_["inside_conc"]->addValuePairs(current_time, ag->getTypeName(), ag->getId(), name, ag->get_inside_conc().at(name));
                }
            }
        }
    }
    else if("dist_2_cells" == active && do_measurement) {
        const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
        auto pos_cell1 = site_->agent_manager_->getAllAgents()[0]->getPosition();
        auto radius = site_->agent_manager_->getAllAgents()[0]->getSurface()->getBasicSphereOfThis()->getRadius();
        for (int k = 0; k < n_x; ++k) {
            std::vector<double> value;
            for (const auto &[name, data]: *site_->molecule_manager_) {
                value.push_back(data[k + n_x * n_y/2 + n_y * n_x * n_z/2]);
            }
            pair_measurements_["dist_2_cells"]->addValuePairs(time.getCurrentTime(),
                                                                      site_->molecule_manager_->get_loc(n_x/2, n_y/2, k).x, value[0], value[1],
                                                              value[2], pos_cell1.x, radius);
        }
    }
    else if ("uptake" == active && do_measurement){
        for (const auto &ag : site_->getAgentManager()->getAllAgents()) {
            if (ag->getTypeName() == "Cell"){
                auto cell = std::dynamic_pointer_cast<Cell>(ag);
                for (const auto& [name, data] : *site_->molecule_manager_) {
                    pair_measurements_["uptake"]->addValuePairs(time.getCurrentTime(), name, cell->get_current_uptake(name));
                }
            }
        }
    }
    // End of the simulation
    if (time.getMaxTime() - time.getCurrentTime() <= time.getCurrentDeltaT()) {
      if ("environment" == active) {
          std::vector<Coordinate3D> sysBound = site_->getSystemBoundaries();
          Coordinate3D pointA = sysBound[0];
          Coordinate3D pointB = sysBound[1];
          pair_measurements_["environment"]->addValuePairs(pointA.x, pointA.y, pointA.z, "0.0",
                                                           "minBounds", 0.0);
          pair_measurements_["environment"]->addValuePairs(pointB.x, pointB.y, pointB.z, "0.0",
                                                           "maxBounds", 0.0);
      }
      else if ("mol_environment" == active){
          const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
          auto const location = site_->molecule_manager_->get_location();
          for (int i= 0; i< n_x*n_y*n_z;++i) {
              Coordinate3D const loc = location[i];
              pair_measurements_["mol_environment"]->addValuePairs(i, loc.x, loc.y, loc.z);
          }
      }
    }
  }
}

void InSituMeasurements::writeToFiles(const std::string &output_dir) const {
   for (const auto&[name, measurement]: histogram_measurements_) {
       const auto file_name = static_cast<boost::filesystem::path>(output_dir).append(name + ".csv").string();
       if (!boost::filesystem::exists(file_name)) {
           std::ofstream file{file_name};

           const auto &keys = measurement->getKeys();
           file << "id" << HistogramMeasurement::delimeter;
           for (const auto &key: keys) {
               file << key << HistogramMeasurement::delimeter;
           }
           file << '\n';
           file << *measurement;
           file.close();
           measurement->clearData();
       } else {
           std::ofstream file{file_name, std::ios_base::app};
           file << *measurement;
           file.close();
           measurement->clearData();
       }
   }
   for (const auto&[name, measurement]: pair_measurements_) {
       const auto file_name = static_cast<boost::filesystem::path>(output_dir).append(name + ".csv").string();
       //DEBUG_STDOUT("file name: " << file_name);
       if (!boost::filesystem::exists(file_name)) {
           std::ofstream file{file_name};
           const auto &keys = measurement->getKeys();
           file << "id" << PairMeasurement::delimeter;
           for (const auto &key: keys) { file << key << PairMeasurement::delimeter;
           }
           file << '\n';
           file << *measurement;
           file.close();
           measurement->clearData();

       } else {
           std::ofstream file{file_name, std::ios_base::app};
           file << *measurement;
           file.close();
           measurement->clearData();
       }
   }
}

void InSituMeasurements::post_simulations_measurements(SimulationTime &time) {
    for (const auto &active_:active_measurements_) {
        // Remove everything after "%" in the active measurement string
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            active = active_.substr(0, char_pos);
        }
        if ("membrane" == active){
            for (auto& agent: site_->getAgentManager()->getAllAgents()){
                for (auto& point: agent->get_membrane()){
                    auto loc = site_->molecule_manager_->get_loc(std::get<0>(point), std::get<1>(point), std::get<2>(point));
                    auto dist = loc.calculateEuclidianDistance(Coordinate3D{0, 0, 0});
                    pair_measurements_["membrane"]->addValuePairs(agent->getId(), agent->getTypeName(), loc.x, loc.y, loc.z, dist);
                }
            }
        }
        else if ("inside" == active) {
            for (auto& agent : site_->getAgentManager()->getAllAgents()) {
                for (auto& point : agent->get_inside()) {
                    auto loc = site_->molecule_manager_->get_loc(std::get<0>(point), std::get<1>(point), std::get<2>(point));
                    pair_measurements_["inside"]->addValuePairs(agent->getId(), agent->getTypeName(), loc.x, loc.y, loc.z);
                }
            }
        }
        else if ("final_conc" == active){
            std::map<std::string, double> avg; // Create a map to store average values
            for (auto &[name, conc]: site_->molecule_manager_->get_conc()){
                double sum = 0.0;
                for (const auto& value : conc) {
                    sum += value;
                }
                double average = sum / conc.size();
                avg.emplace(name, average);
            }
            pair_measurements_["final_conc"]->addValuePairs(avg["AMP"], avg["Complex"], avg["Defensive"]);
        }
        else if ("molecules" == active) {
            double value;
            auto const location = site_->molecule_manager_->get_location();
            const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
            for (const auto& [name, data] : *site_->molecule_manager_) {
                for (int i= 0; i< n_x*n_y*n_z;++i){
                    value = data[i];
                    auto const loc = location[i];
                    pair_measurements_["molecules"]->addValuePairs(name, time.getCurrentTime(), loc.x, loc.y, loc.z, value);
                }
            }
        }
        else if ("molecules_test" == active){
            auto const location = site_->molecule_manager_->get_location();
            const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
            for (int i= 0; i< n_x*n_y*n_z;++i){
                std::vector<double> value;
                for (const auto& [name, data] : *site_->molecule_manager_) {
                    value.push_back(data[i]);
                }
                auto const loc = location[i];
                //SYSTEM_STDOUT("time: " << current_time << "| location: " << loc.x << ", " << loc.y << ", " << loc.z << "| values: " << value[0] << ", " << value[1] << ", " << value[2]);
                pair_measurements_["molecules_test"]->addValuePairs(time.getCurrentTime(), loc.x, loc.y, loc.z, value[0], value[1], value[2]);
            }
        }
        else if("spatial_dist_center" == active) {
            const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
            auto dist = 0.0;
            for (int i = (int)round(n_z / 2); i < n_z; ++i) {
                std::vector<double> value;
                for (const auto& [name, data] : *site_->molecule_manager_) {
                    value.push_back(data[n_x / 2 + n_x * n_y / 2 + n_y * n_x * i]); // center of the simulation site
                    Coordinate3D position = {0, 0, 0};
                    dist = position.calculateEuclidianDistance(site_->molecule_manager_->get_loc(i, n_y/2, n_x/2));
                }
                pair_measurements_["spatial_dist_center"]->addValuePairs(time.getCurrentTime(), dist, value[0], value[1], value[2]);
            }
        }
        else if("spatial_dist_center_updated" == active) {
            const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
            auto dist = 0.0;
            for (int i = round(n_z / 2); i < n_z; ++i) {
                std::vector<double> value;
                for (const auto& [name, data] : *site_->molecule_manager_) {
                    value.push_back(data[n_x / 2 + n_x * n_y / 2 + n_y * n_x * i]); // center of the simulation site
                    Coordinate3D position = {0, 0, 0};
                    dist = position.calculateEuclidianDistance(site_->molecule_manager_->get_loc(i, n_y/2, n_x/2));
                }
                auto inside = 0.0;
                if (dist < 3.4){
                    for (const auto ag : site_->getAgentManager()->getAllAgents()) {
                        if (ag->get_track_inside_conc()) {
                            inside = ag->get_inside_conc().at("AMP");
                        }
                    }
                }
                pair_measurements_["spatial_dist_center_updated"]->addValuePairs(time.getCurrentTime(), dist, value[0], value[1], value[2], inside);
            }
        }
        else if ("mol" == active){
            double value;
            const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
            for (auto& [name, data] : *site_->molecule_manager_) {
                for (int i = 0; i < n_x * n_y * n_z; ++i) {
                    value = data[i];
                    pair_measurements_["mol"]->addValuePairs(name, time.getCurrentTime(), i, value);
                }
            }
        }
        else if ("agent-statistics" == active) {
            for (const auto agent: site_->getAgentManager()->getAllAgents()) {
                std::string cellpart = "Mothercell";
                for (const auto cellparts: agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis()) {
                    pair_measurements_["agent-statistics"]->addValuePairs(time.getCurrentTime(),
                                                                          agent->getTypeName(), agent->getId(),
                                                                          cellparts->getRadius(),
                                                                          cellparts->getPosition().x,
                                                                          cellparts->getPosition().y,
                                                                          cellparts->getPosition().z,
                                                                          cellpart, cellparts->getDescription());
                }
            }
        }
        else if ("inside_conc" == active) {
            for (const auto ag : site_->getAgentManager()->getAllAgents()) {
                if (ag->get_track_inside_conc()){
                    for (const auto& [name, data] : *site_->molecule_manager_) {
                        pair_measurements_["inside_conc"]->addValuePairs(time.getCurrentTime(), ag->getTypeName(), ag->getId(), name, ag->get_inside_conc().at(name));
                    }
                }
            }
        }

        else if ("mod_distrib" == active) {
            const auto& [n_x, n_y, n_z] = site_->molecule_manager_->getGridSize();
            std::vector<double> dist;
            for (const auto& [name, data] : *site_->molecule_manager_) {
                auto max_value = 0.0;
                auto tmp = 0.0;
                auto max_dist = 0;
                for (int i = (int) round(n_z / 2); i < n_z; ++i) {
                    tmp = data[n_x / 2 + n_x * n_y / 2 + n_y * n_x * i];
                    if (tmp > max_value) {
                        max_value = tmp;
                        max_dist = i;
                    }
                }
                Coordinate3D position = {0, 0, 0};
                dist.push_back(position.calculateEuclidianDistance(site_->molecule_manager_->get_loc(max_dist, n_y/2, n_x/2)));
            }
            pair_measurements_["mod_distrib"]->addValuePairs(time.getCurrentTime(), dist[0], dist[1], dist[2]);
        }
    }
}
