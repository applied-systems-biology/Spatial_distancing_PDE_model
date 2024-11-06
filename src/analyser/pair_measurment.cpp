//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "analyser/pair_measurement.h"

std::ostream &operator<<(std::ostream &out, PairMeasurement &measurement) {
  if (!measurement.data_.empty()) {
    for (const auto &line: measurement.data_) {
      out << measurement.owner_id_ << PairMeasurement::delimeter << line << '\n';
    }
  }
  return out;
}

