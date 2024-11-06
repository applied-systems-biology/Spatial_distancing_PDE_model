//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "analyser/histogram_measurment.h"

#include <sstream>

std::ostream &operator<<(std::ostream &out, HistogramMeasurement &measurement) {
  if (!measurement.data_.empty()) {
    auto cache_access = abm::util::overloaded{
        [&](std::monostate /*unused*/)->std::string { return ""; }, [&](auto arg)->std::string { return std::to_string(arg); }
        , [&](std::string arg)->std::string { return arg; },};

    std::uint64_t max_size = 0;
    for (const auto &key: measurement.keys_) {
      max_size = std::max(max_size, measurement.data_[key].size());
    }
    for (std::uint64_t i = 0; i < max_size; ++i) {
      out << measurement.owner_id_<< HistogramMeasurement::delimeter;
      for (const auto &key: measurement.keys_) {
        if (measurement.data_[key].size() > i) {
          out << std::visit(cache_access, measurement.data_[key][i]) << HistogramMeasurement::delimeter;
        } else {
          out << HistogramMeasurement::delimeter;
        }
      }
      out << '\n';
    }
  }
  return out;
}
