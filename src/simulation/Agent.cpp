//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "simulation/Agent.h"

#include "simulation/Site.h"
#include <cmath>
#include "simulation/AgentManager.h"

Agent::Agent() {
  agentProps = std::make_unique<AgentProperties>();

  id = 0;
  is_deleted_ = false;
  initialPosition = std::make_unique<Coordinate3D>();
  setInitialPosition(getPosition());
}

Agent::Agent(std::unique_ptr<Coordinate3D> c , int id , Site* site){

  agentProps = std::make_unique<AgentProperties>();

  this->id = id;
  is_deleted_ = false;
  this->site = site;

  initialPosition = std::make_unique<Coordinate3D>(Coordinate3D{c->x, c->y, c->z});
  position = std::move(c);

  setInitialPosition(getPosition());

}


int Agent::getId() const{
    return id;
}

Coordinate3D Agent::getPosition(){
    
    return *position;
}

Site* Agent::getSite(){
    return site;
}

void Agent::setInitialPosition(Coordinate3D initPos){
    *initialPosition = initPos;
}

AgentProperties* Agent::getAgentProperties(){
    return agentProps.get();
}

void Agent::setDeleted(){
  is_deleted_ = true;
    
}


std::string Agent::getAgentCSVTagAoc(double current_time){
    std::ostringstream csvTag;
    csvTag << current_time << "\t" <<
           position->x << "\t" <<
           position->y << "\t" <<
           position->z<< "\t" <<
           getTypeName();
    return csvTag.str();
}

std::string Agent::getAgentCSVTagToc(){

    std::ostringstream csvTag;
    csvTag << id << "\t" <<
           position->x << "\t" <<
           position->y << "\t" <<
           position->z << "\t" <<
           getTypeName();

    return csvTag.str();
}
