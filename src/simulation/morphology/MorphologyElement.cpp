//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "MorphologyElement.h"

const std::vector<std::shared_ptr<SphereRepresentation>> &MorphologyElement::getSphereRepresentation() {
  return sphere_rep_;
}

std::string MorphologyElement::getDescription() {
  return description_;
}

Morphology *MorphologyElement::getMorphologyThisBelongsTo() {
  return
      morphology_this_belongs_to_;
}

double MorphologyElement::getVolume() const {
  return volume_;
}