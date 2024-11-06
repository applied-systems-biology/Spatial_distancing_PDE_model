//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef _AGENT_H
#define	_AGENT_H

#include <iostream>

#include "basic/Coordinate3D.h"

#include "map"
#include "simulation/Morphology.h"
#include "simulation/AgentProperties.h"



class Site; //forward declaration

class Agent {
 public:
    // Abstract class for agents in the hABM. This class provides the cell class with it main functionality.

    Agent();
  Agent(std::unique_ptr<Coordinate3D>, int, Site *);
  virtual ~Agent() = default;

  virtual void doAllActionsForTimestep(double timestep, double current_time) = 0;
  virtual std::string getTypeName() = 0;

  virtual std::string getAgentCSVTagAoc(double current_time) = 0;
  virtual std::string getAgentCSVTagToc() = 0;

  virtual Morphology *getSurface() = 0;

  Coordinate3D getPosition();

  void setDeleted();
  [[nodiscard]] int getId() const;
  [[nodiscard]] bool isDeleted() const { return is_deleted_;}

  Site *getSite();
  AgentProperties *getAgentProperties();

  Coordinate3D getInitialPosition() { return *initialPosition; };
  Coordinate3D getCurrentPosition() { return *position; };

  std::map<std::string, double> molecule_uptake;

  virtual   std::map<std::string, double> get_inside_conc() const = 0;

  virtual bool get_track_inside_conc() = 0;

  virtual   std::vector<std::tuple<int, int, int>> get_membrane() const = 0;
  virtual   std::vector<std::tuple<int, int, int>> get_inside() const = 0;

private:

 protected:
  void setInitialPosition(Coordinate3D initPos);

  unsigned int id{};
  bool is_deleted_;

  Site *site;

  std::unique_ptr<AgentProperties> agentProps;
  std::unique_ptr<Coordinate3D> initialPosition;
  std::shared_ptr<Coordinate3D> position;

};

#endif	/* _AGENT_H */

