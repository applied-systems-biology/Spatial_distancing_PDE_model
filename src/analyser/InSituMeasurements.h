//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef INSITUMEASUREMENTS_H
#define INSITUMEASUREMENTS_H

#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <memory>
#include <set>

#include "utils/time_util.h"
#include "analyser/histogram_measurment.h"
#include "analyser/pair_measurement.h"
#include "utils/macros.h"

class Site;
class Collision;
class Cell;
class Coordinate3D;
namespace abm::test{void test_insitu_measurements();}

class InSituMeasurements {
 public:
  // Class for providing the main functionality of taking measurements during a simulation run
  InSituMeasurements() = default;
  InSituMeasurements(std::unordered_set<std::string> active_measurements, const std::string &, double  every_X_time_unit);
  void observeMeasurements(const SimulationTime &time);
  void setSite(Site *site) { site_ = site; }
  void writeToFiles(const std::string &output_dir) const;

  void post_simulations_measurements(SimulationTime &time);

  template<typename T, typename U>
  void increment(const std::string &measurement_name, const std::string &curve_name, U value) {
    if constexpr(std::is_same_v<T, HistogramMeasurement>) {
      const auto &measurement = histogram_measurements_.find(measurement_name);
      if (measurement != histogram_measurements_.end()) {
        auto cache_access = abm::util::overloaded{
            [&](std::monostate /*unused*/)->HistogramMeasurement::cache_type { return value; },
            [&](auto &arg)->HistogramMeasurement::cache_type {
              return value + arg;
            },};
        measurement->second->cache[curve_name] = std::visit(cache_access, measurement->second->cache[curve_name]);
      }
    } else if constexpr(std::is_same_v<T, PairMeasurement>) {
      const auto &measurement = pair_measurements_.find(measurement_name);
      if (measurement != pair_measurements_.end()) {
        auto cache_access = abm::util::overloaded{
            [&](std::monostate /*unused*/)->PairMeasurement::cache_type { return value; },
            [&](auto &arg)->PairMeasurement::cache_type { return value + arg;},
            [&](std::vector<int> &arg)->PairMeasurement::cache_type { arg.emplace_back(value);return std::move(arg); },
            [&](std::vector<double> &arg)->PairMeasurement::cache_type { arg.emplace_back(value);return std::move(arg); }
            ,};
        measurement->second->cache[curve_name] = std::visit(cache_access, measurement->second->cache[curve_name]);
      }
    }
  }
  template<typename T, typename ...ARGS>
  void addValues(const std::string &measurement_name, ARGS &&...args) {
    if constexpr(std::is_same_v<T, HistogramMeasurement>) {
      const auto &measurement = histogram_measurements_.find(measurement_name);
      if (measurement != histogram_measurements_.end()) {
        measurement->second->addValues(std::forward<ARGS>(args)...);
      }
    } else if constexpr(std::is_same_v<T, PairMeasurement>) {
      const auto &measurement = pair_measurements_.find(measurement_name);
      if (measurement != pair_measurements_.end()) {
        measurement->second->addValuePairs(std::forward<ARGS>(args)...);
      }
    }
  }
 friend  void abm::test::test_insitu_measurements();
 private:
  Site *site_{};
  std::unordered_set<std::string> active_measurements_;
  std::unordered_map<std::string, std::unique_ptr<HistogramMeasurement>> histogram_measurements_{};
  std::unordered_map<std::string, std::unique_ptr<PairMeasurement>> pair_measurements_{};
  double every_X_time_unit_{};

};

#endif	/* INSITUMEASUREMENTS_H */

