//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "simulation/CellFactory.h"

#include "simulation/Site.h"

CellFactory::CellFactory(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters) {
    for (const auto &agent: site_parameters->agent_manager_parameters.agents) {
        agent_configurations_.emplace(std::make_pair(site_parameters->identifier + agent->type, agent));
    }
}


std::shared_ptr<Cell> CellFactory::createCell(const std::string &agenttype,
                                              std::unique_ptr<Coordinate3D> c,
                                              int id,
                                              Site *site,
                                              double time_delta,
                                              double current_time) {
  std::shared_ptr<Cell> agent{};
  auto site_tag = site->getIdentifier();
  if (agenttype == "Cell") {
      agent = std::make_shared<Cell>(std::move(c), id, site, time_delta, current_time);
  }
  agent->setup(time_delta, current_time,agent_configurations_[site_tag + agenttype].get());
  return agent;
}

