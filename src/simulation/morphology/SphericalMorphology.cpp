//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include <list>

#include "SphericalMorphology.h"
#include "SphereRepresentation.h"
#include "simulation/Site.h"

SphericalMorphology::SphericalMorphology() : MorphologyElement(){
}

SphericalMorphology::SphericalMorphology(Morphology* refMorphology, std::shared_ptr<Coordinate3D> position, double radius, std::string description) : MorphologyElement(refMorphology, description){
  center_of_mass_ = position;
    this->radius = radius;
    this->virtualRadiusExtension = 0;
    this->description_ = description;
    
    generateSphereRepresentation();
}

SphericalMorphology::SphericalMorphology(Morphology* refMorphology, std::shared_ptr<Coordinate3D> position, double radius, double virtualRadiusExtension, std::string description) : MorphologyElement(refMorphology, description){
  center_of_mass_ = position;
    this->radius = radius;
    this->virtualRadiusExtension = virtualRadiusExtension;
    this->description_ = description;
    
    generateSphereRepresentation();
}


void SphericalMorphology::generateSphereRepresentation(){
    
    auto sR = std::make_unique<SphereRepresentation>(center_of_mass_, radius, virtualRadiusExtension, this, description_);

    sphere_rep_.push_back(std::move(sR));

}