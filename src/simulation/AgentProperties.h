//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef AGENTPROPERTIES_H
#define	AGENTPROPERTIES_H
#include <iostream>

#include "simulation/Morphology.h"

class Interactions;
/**
 * this is a kind of an interface class, to promote access to specialised elements of subclasses 
 * via the original agent-class, just one compendium of getters and setters towards the corresponding variables
 */
class AgentProperties {

 public:
  AgentProperties() = default;

  void setMorphology(std::shared_ptr<Morphology> morphology);

  Morphology *getMorphology();
 private:

  std::shared_ptr<Morphology> morphology_{};
};

#endif	/* AGENTPROPERTIES_H */

