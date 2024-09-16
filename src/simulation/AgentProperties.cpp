//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "simulation/AgentProperties.h"

void AgentProperties::setMorphology(std::shared_ptr<Morphology>  morphology){
  morphology_ = std::move(morphology);
}

Morphology *AgentProperties::getMorphology() {
  if (morphology_ == nullptr) {
    morphology_ = std::make_unique<Morphology>();
  }

  return morphology_.get();
}