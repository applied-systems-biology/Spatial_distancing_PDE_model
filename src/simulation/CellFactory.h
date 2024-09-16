//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef CELLFACTORY_H
#define	CELLFACTORY_H

#include <string>
#include <array>
#include <map>

#include "simulation/Cell.h"
#include "utils/io_util.h"
#include "Site.h"


/**
 * class that generates Cells from given input data
 */

class CellFactory {
public:
    // Factory class for initializing new agents according to the simulator configuration.
    CellFactory(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters);

    virtual std::shared_ptr<Cell> createCell(const std::string &agenttype,
                                             std::unique_ptr<Coordinate3D> c,
                                             int id,
                                             Site *site,
                                             double time_delta,
                                             double current_time);


protected:
    std::map<std::string, std::shared_ptr<abm::util::SimulationParameters::AgentParameters>> agent_configurations_;
};


#endif	/* CELLFACTORY_H */