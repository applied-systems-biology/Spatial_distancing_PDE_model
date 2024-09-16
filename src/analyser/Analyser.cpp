//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "analyser/Analyser.h"

#include "analyser/InSituMeasurements.h"

Analyser::Analyser(const std::string &config_path, const std::string &project_dir) {
  auto const analyser_config = static_cast<boost::filesystem::path>(config_path).append("analyser-config.json");
  parameters_ = abm::util::getAnalyserParameters(analyser_config.string());
  if (!parameters_.active_measurements.empty()) {
    measurement_path_ = static_cast<boost::filesystem::path>(project_dir).append("measurements").string();
    boost::filesystem::create_directories(measurement_path_);
  }
}

std::shared_ptr<InSituMeasurements> Analyser::generateMeasurement(const std::string& id) const{
  std::shared_ptr<InSituMeasurements> new_measurement = nullptr;
#pragma omp critical
  {
    new_measurement = measurements_.emplace_back(std::make_shared<InSituMeasurements>(parameters_.active_measurements, id, parameters_.output_every_X_time_unit));
  }
  return new_measurement;
}

void Analyser::outputAllMeasurements() const {
  for (const auto &measurement: measurements_) {
   measurement->writeToFiles(measurement_path_);
  }
}
