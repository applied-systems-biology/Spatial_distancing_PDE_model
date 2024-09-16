//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.



#include "simulation/AgentManager.h"

#include <boost/algorithm/string.hpp>
#include <cmath>

#include "simulation/CellFactory.h"
#include "analyser/InSituMeasurements.h"
#include "utils/macros.h"
#include "simulation/Site.h"


AgentManager::AgentManager(double time_delta, Site *site) {
  this->site = site;
  time_delta_ = time_delta;
  idHandling = 0;
  idHandlingSphereRepresentation = 0;
}

Agent *AgentManager::createAgent(Site *site,
                                 std::string agenttype,
                                 Coordinate3D c,
                                 double current_time) {
  auto agent = emplace_back(site->getCellFactory()->createCell(agenttype,
                                                    std::make_unique<Coordinate3D>(c),
                                                    generateNewID(),
                                                    site,
                                                    time_delta_,
                                                    current_time));
  return agent.get();
}

void AgentManager::replaceAgent(Site *site,
                                Agent *agent,
                                std::unique_ptr<Coordinate3D> newCoord,
                                Coordinate3D *prevMove,
                                double current_time) {
  std::string agentType = agent->getTypeName();
  auto newAgent = site->getCellFactory()->createCell(agentType,
                                          std::move(newCoord),
                                          idHandling,
                                          site,
                                          time_delta_,
                                          current_time);
  if (newAgent != 0) {
    std::replace_if(allAgents.begin(), allAgents.end(), [agent](const auto &a) { return agent == a.get(); }, newAgent);
    idHandling++;
  }
  agent->setDeleted();
}

void AgentManager::replaceAgent(Site *site, Agent *agentToReplace, std::shared_ptr<Agent> newAgent, double current_time) {
  std::string name = agentToReplace->getTypeName();
  std::string siteName = site->getType();

  std::replace_if(allAgents.begin(),
                  allAgents.end(),
                  [agent = agentToReplace](const auto &a) { return agent == a.get(); },
                  newAgent);

  if (newAgent == 0) {
  }
  agentToReplace->setDeleted();
}

void AgentManager::removeAgent(Site *site, Agent *agent, double current_time) {

  agent->setDeleted();

  allAgents.erase(std::remove_if(allAgents.begin(),
                                 allAgents.end(),
                                 [agent](const auto &a) { return agent == a.get(); }), allAgents.end());
}

int AgentManager::getAgentQuantity(std::string agenttype) {
  int count = 0;
  for (auto agent: allAgents){
    if (agent != 0){
      if (agent->getTypeName() == agenttype){
        count++;
      }
    }
  }
  return count;
}

const std::vector<std::shared_ptr<Agent>> &AgentManager::getAllAgents() {
  return allAgents;
}

std::vector<Agent *> AgentManager::getAllAgentsOfOneCelltype(std::string agenttype) {
  std::vector<Agent *> agents(getAgentQuantity(agenttype));
  const auto &all = getAllAgents();

  size_t i = 0;
  while (i < agents.size()) {
    for (size_t j = 0; j < all.size(); j++) {
      if (all[j] != 0) {
        if (all[j]->getTypeName().compare(agenttype) == 0) {
          agents[i] = all[j].get();
          i++;
        }
      }

    }
  }
  return agents;
}

void AgentManager::cleanUpAgents() {
  for (auto it = allAgents.begin(); it != allAgents.end();) {
    if (*it == 0) {
      it = allAgents.erase(it);

    } else {
      if ((*it)->isDeleted()) {
        for (const auto &sphere: (*it)->getAgentProperties()->getMorphology()->getAllSpheresOfThis()) {
          removeSphereRepresentation(sphere);
        }
        it = allAgents.erase(it);
      } else {
        it++;
      }
    }
  }

}

int AgentManager::getNextSphereRepresentationId(SphereRepresentation *sphereRep) {
  sphereIdToCell[idHandlingSphereRepresentation] =
      sphereRep->getMorphologyElementThisBelongsTo()->getMorphologyThisBelongsTo()->getCellThisBelongsTo();
  sphereIdToSphereRep[idHandlingSphereRepresentation] = sphereRep;

  //tell the neighbourhood locator that here is a new sphere included
  allSphereRepresentations.insert(sphereRep);

  return idHandlingSphereRepresentation++;
}
//cout << "all agent: " << count << " not moved:" << notMoved << '\n';

Cell *AgentManager::getCellBySphereRepId(int sphereRepId) {
  if (sphereIdToCell.find(sphereRepId) != sphereIdToCell.end()) {
    return sphereIdToCell[sphereRepId];
  }
  return NULL;
}

SphereRepresentation *AgentManager::getSphereRepBySphereRepId(int sphereRepId) {
  if (sphereIdToSphereRep.find(sphereRepId) != sphereIdToSphereRep.end()) {
    return sphereIdToSphereRep[sphereRepId];
  }
  return NULL;
}

void AgentManager::removeSphereRepresentation(SphereRepresentation *sphereRep) {
  allSphereRepresentations.erase(sphereRep);
  sphereIdToSphereRep.erase(sphereRep->getId());
  sphereIdToCell.erase(sphereRep->getId());
}

int AgentManager::getIdHandling() const {
  return idHandling;
}

void AgentManager::incrementIdHandling() {
  idHandling++;
}
