//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef MORPHOLOGYELEMENT_H
#define	MORPHOLOGYELEMENT_H

#include <utility>
#include <vector>



#include "basic/Coordinate3D.h"
#include "simulation/morphology/SphereRepresentation.h"

class Morphology;

class MorphologyElement {
 public:
    // Abstract class for morphology representation
    MorphologyElement() = default;
  MorphologyElement(Morphology *morphologyRef, std::string description)
      : morphology_this_belongs_to_(morphologyRef), description_(std::move(description)) {};
  virtual ~MorphologyElement() = default;
  virtual void generateSphereRepresentation() = 0;
  const std::vector<std::shared_ptr<SphereRepresentation>> &getSphereRepresentation();
  std::string getDescription();
  Morphology *getMorphologyThisBelongsTo();
  double getVolume() const;

 protected:
  std::shared_ptr<Coordinate3D> center_of_mass_{};
  double volume_{};
  Morphology *morphology_this_belongs_to_{};
  std::string description_;
  std::vector<std::shared_ptr<SphereRepresentation>> sphere_rep_;

};

#endif	/* MORPHOLOGYELEMENT_H */

