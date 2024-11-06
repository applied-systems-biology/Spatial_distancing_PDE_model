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

#include "simulation/Morphology.h"
#include "simulation/Cell.h"
#include "simulation/morphology/SphereRepresentation.h"
#include "simulation/Site.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif


Morphology::Morphology(Cell *cell) {
  this->cell_this_belongs_to_ = cell;
}

void Morphology::appendAssociatedCellpart(std::unique_ptr<MorphologyElement> morphElement){
    morphologyElements.emplace_back(std::move(morphElement));
    //cout << "[Morphology] elements of cell : " << cellThisBelongsTo->getId() << " " << morphologyElements.size() << '\n';;
}

std::vector<SphereRepresentation*> Morphology::getAllSpheresOfThis(){
  std::vector<SphereRepresentation*> allSpheresOfThisCell;
    
    auto it = morphologyElements.begin();
    
    while(it != morphologyElements.end()){
      std::vector<std::shared_ptr<SphereRepresentation>> sphereReps = (*it)->getSphereRepresentation();
        //cout << "[Morphology] size: " << sphereReps.size() << '\n';;
        for(const auto& ptr: sphereReps){
          allSpheresOfThisCell.emplace_back(ptr.get());
        }
        it++;
    }
    
    return allSpheresOfThisCell;
}

SphereRepresentation* Morphology::getBasicSphereOfThis(){
  std::vector<SphereRepresentation*> allSpheresOfThisCell = getAllSpheresOfThis();
  auto it = allSpheresOfThisCell.begin();
    SphereRepresentation *basicSphere = 0;
    while(it != allSpheresOfThisCell.end()){
        if((*it)->getMorphologyElementThisBelongsTo()->getDescription() == "basic"){
            basicSphere = *it;
            break;
        }
        it++;
    }
    
    return basicSphere;
}
