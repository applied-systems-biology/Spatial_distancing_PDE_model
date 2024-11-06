//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef ANALYSER_H
#define ANALYSER_H

#include <boost/filesystem.hpp>
#include <memory>

#include "utils/io_util.h"

class InSituMeasurements;

class Analyser {

 public:
    // Class for wrapping analyzer functionality to take measurements during a simulation run
    Analyser() = default;
    ~Analyser() = default;
    Analyser(const std::string &config_path, const std::string &project_dir);
    Analyser(const Analyser &) = delete;
    Analyser &operator=(const Analyser &) = delete;
    Analyser(Analyser &&) = delete;
    Analyser &operator=(Analyser &&) = delete;

    std::shared_ptr<InSituMeasurements> generateMeasurement(const std::string& id) const;
    void outputAllMeasurements() const;

private:
    mutable std::vector<std::shared_ptr<InSituMeasurements>> measurements_;
    std::string measurement_path_;
    abm::util::AnalyserParameters parameters_;
};

#endif /* ANALYSER_H */
