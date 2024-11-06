//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "testConfigurations.h"

#include <memory>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "external/doctest/doctest.h"
#include "simulation/Site.h"
#include "simulation/Simulator.h"
#include "analyser/Analyser.h"

using boost::filesystem::path;
using boost::filesystem::exists;

std::string abm::test::test_simulation(const std::string &config) {
  const auto parameters = abm::util::getMainConfigParameters(config);
    const auto simulator = std::make_unique<Simulator>();
    simulator->setConfigPath(parameters.config_path);
    auto run_seed = parameters.system_seed;
    const auto analyser = std::make_unique<const Analyser>();
    const auto random_generator = std::make_unique<Randomizer>(run_seed);
    const auto sites = simulator->createSites(run_seed, random_generator.get(), analyser.get());
    const auto &site = sites;
    SimulationTime time{site->getTimeStepping(), site->getMaxTime()};
    for (time.updateTimestep(0); !time.endReached(); ++time) {    //inner loop: t -> t + dt
        site->doAgentDynamics(random_generator.get(), time);
    }

  return abm::util::generateHashFromMolecules(time.getCurrentTime(), site->molecule_manager_->get_conc());
}
// Spatial distancing Model Test
TEST_CASE ("Check Spatial Distancing Test") {
  path config("../../test/configurations/testSpatialDistancing/config.json");
        CHECK(exists(config) == true);
  const auto string_return = abm::test::test_simulation(config.string());
        CHECK(string_return == "9120002391968631136");
}

