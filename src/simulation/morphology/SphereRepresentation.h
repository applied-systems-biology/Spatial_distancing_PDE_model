//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef SPHEREREPRESENTATION_H
#define	SPHEREREPRESENTATION_H

#include <memory>
#include <sstream>
#include "basic/Coordinate3D.h"

class MorphologyElement;

class SphereRepresentation {
public:
    // Class for a spherical representation based on a spherical morphology.
    SphereRepresentation();
    
    SphereRepresentation(std::shared_ptr<Coordinate3D> coord, double radius, MorphologyElement* morphologyElement, std::string description="");
    SphereRepresentation(std::shared_ptr<Coordinate3D> coord, double radius, double virtualRadiusExtension, MorphologyElement* morphologyElement, std::string description="");
    
    virtual ~SphereRepresentation();

    Coordinate3D getPosition(){return *position;};
    
    double getRadius(){return radius;};
    void   setRadius(double r){ radius=r;};
    double getVirtualRadiusExtension(){return virtualRadiusExtension;};

    double getTotalRadius(){return radius + virtualRadiusExtension;};

    int getId(){return id;};

    std::string getDescription(){return description_;};

    MorphologyElement* getMorphologyElementThisBelongsTo(){return morphologyElementThisBelongsTo;};

    Coordinate3D getEffectiveConnection(SphereRepresentation* sphereRep);

    void shiftPosition(Coordinate3D* shifter);

    void setRadiusToOrigin(double r);

    void setRandomizedRadius(double time);

private:

    std::shared_ptr<Coordinate3D> position;

    double radius;

    double virtualRadiusExtension;

    std::string description_;

    /**
     * the morphological element this sphere is a representation of
     */
    MorphologyElement* morphologyElementThisBelongsTo;

    /**
     * a unique identifier of this sphere
     */
    int id;


    //=========HERE COMES EXPERIMENTAL STUFF===============//
    /*
    double tPeriod;
    double rDev;
    double r0;
     * */
};

#endif	/* SPHEREREPRESENTATION_H */

