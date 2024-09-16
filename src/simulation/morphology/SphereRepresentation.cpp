//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "SphereRepresentation.h"

#include "MorphologyElement.h"
#include "simulation/Cell.h"
#include "simulation/Site.h"
#include "simulation/AgentManager.h"

SphereRepresentation::SphereRepresentation() {
}

SphereRepresentation::SphereRepresentation(std::shared_ptr<Coordinate3D> coord, double radius,  MorphologyElement* morphologyElement, std::string description){
    this->position = coord;
    this->morphologyElementThisBelongsTo = morphologyElement;

    this->radius = radius;
    this->virtualRadiusExtension = 0;
    this->description_ = description;

    id = morphologyElementThisBelongsTo->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getSite()->getAgentManager()->getNextSphereRepresentationId(this);

    if(!morphologyElementThisBelongsTo->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getSite()->containsPosition(*position)){
      std::cout << "[SphereRepresentation] Warning: Sphere at " << position->x << " is not within site!; t=" << "\n";
    }

}

SphereRepresentation::SphereRepresentation(std::shared_ptr<Coordinate3D> coord, double radius, double virtualRadiusExtension, MorphologyElement* morphologyElement, std::string description){
    this->position = coord;
    this->morphologyElementThisBelongsTo = morphologyElement;

    this->radius = radius;
    this->virtualRadiusExtension = virtualRadiusExtension;
    this->description_ = description;
    id = morphologyElementThisBelongsTo->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getSite()->getAgentManager()->getNextSphereRepresentationId(this);

    if(!morphologyElementThisBelongsTo->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getSite()->containsPosition(*position)){
      std::cout << "[SphereRepresentation] Warning: Sphere at " << position->x << " is not within site!; "<< "\n";
    }

}

SphereRepresentation::~SphereRepresentation() {
}

Coordinate3D SphereRepresentation::getEffectiveConnection(SphereRepresentation* sphereRep){
  return Coordinate3D(sphereRep->getPosition()-*position);
}

void SphereRepresentation::shiftPosition(Coordinate3D* shifter){
  *position += *shifter;
}

void SphereRepresentation::setRandomizedRadius(double time){
  //radius = r0 + rDev*sin(2*M_PI*time/tPeriod);
}