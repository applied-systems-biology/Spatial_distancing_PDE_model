//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef AGENTMANAGER_H
#define AGENTMANAGER_H

#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>


#include "simulation/morphology/SphereRepresentation.h"

class Site;
class Analyser;
class Agent;
class Cell;

class AgentManager {
 public:
    // Class for managing all the agents in data container. It provides functionality for input and replacement of agents during a simulation.
  AgentManager(double time_delta, Site* site);

  // Iterator support for AgentManager
  using iterator = std::vector<std::shared_ptr<Agent>>::iterator;
  using const_iterator = std::vector<std::shared_ptr<Agent>>::const_iterator;
  iterator begin() { return allAgents.begin(); }
  [[nodiscard]] const_iterator begin() const { return allAgents.begin(); }
  iterator end() { return allAgents.end(); }
  [[nodiscard]] const_iterator end() const { return allAgents.end(); }
  std::shared_ptr<Agent> &emplace_back(std::shared_ptr<Agent> &&value) {
    return allAgents.emplace_back(std::forward<std::shared_ptr<Agent>>(value));
  }


  int generateNewID(){
    return idHandling++;
  }

  /**
   * create an agent at the initial coordinate
   * @XMLNode information about the agent that should be created
   * @string name of agenttype that should be created
   * @Coordinate3D initial coordinate/position
   */
  Agent *createAgent(Site *, std::string, Coordinate3D,double current_time);

  /**
     * replace an agent in the old array by a new one
     * @param agent to remove
     * @param coordinate of the new position
     * @param vector of the previous movement which was done
     */
  void replaceAgent(Site *site, Agent *, std::unique_ptr<Coordinate3D>, Coordinate3D *, double current_time);

  /**
   * replace an agent in the old array by a new one
   * @param agent to replace
   * @param agent that replaces the old one
   */
  void replaceAgent(Site *site, Agent *agentToReplace, std::shared_ptr<Agent> newAgent, double current_time);

  /**
   * removes one specified agent by its pointer from this site
   * @param agent that should be removed from this site
   */
  void removeAgent(Site *site, Agent *agent, double current_time);

  /**
   * retrieve the current quantity of a given agenttype
   * @param agenttype name of the agent
   * @return quantity of agenttype in this site
   */
  int getAgentQuantity(std::string agenttype);

  const std::vector<std::shared_ptr<Agent>> &getAllAgents();


  std::vector<Agent *> getAllAgentsOfOneCelltype(std::string agenttype);

  void cleanUpAgents();

  int getNextSphereRepresentationId(SphereRepresentation *sphereRep);

  Cell *getCellBySphereRepId(int sphereRepId);

  SphereRepresentation *getSphereRepBySphereRepId(int sphereRepId);

  std::set<SphereRepresentation *> *getAllSphereRepresentations() { return &allSphereRepresentations; };

  void removeSphereRepresentation(SphereRepresentation *sphereRep);

  int getIdHandling() const;

  void incrementIdHandling();

  int system{};
 private:

  /**
   * number of agents of a specific type in the current site
   */

  std::vector<std::shared_ptr<Agent>> allAgents;

  std::vector<Agent *> deathAgents;

  std::map<int, Cell *> sphereIdToCell;

  std::map<int, SphereRepresentation *> sphereIdToSphereRep;

  std::set<SphereRepresentation *> allSphereRepresentations;

  int idHandling;

  int idHandlingSphereRepresentation;

  Site *site{};
  double time_delta_;
};

#endif    /* AGENTMANAGER_H */

