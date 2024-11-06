//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef SURFACE_H
#define	SURFACE_H

#include <list>

#include "basic/Coordinate3D.h"
#include "simulation/morphology/MorphologyElement.h"


class Cell;

class Morphology {
 public:
    // Class for handling the morphology of cells

    Morphology() = default;

  Morphology(Cell *cell);
  virtual ~Morphology() = default;

  Cell *getCellThisBelongsTo() { return cell_this_belongs_to_; };
  void appendAssociatedCellpart(std::unique_ptr<MorphologyElement> morphElement);
  std::vector<SphereRepresentation *> getAllSpheresOfThis();
  SphereRepresentation *getBasicSphereOfThis();

 protected:
  std::list<std::unique_ptr<MorphologyElement>> morphologyElements;
  Cell *cell_this_belongs_to_ {};

};

#endif	/* SURFACE_H */

