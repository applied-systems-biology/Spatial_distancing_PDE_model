//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "basic/Coordinate3D.h"
#include "utils/macros.h"
#include "basic/Randomizer.h"

double Coordinate3D::calculateEuclidianDistance(const Coordinate3D &coord) const noexcept{
    const double dx = x - coord.x;
    const double dy = y - coord.y;
    const double dz = z - coord.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double Coordinate3D::calculateEuclidianDistancePeriodic(const Coordinate3D& coord,  Coordinate3D& site_limit) const noexcept {
    double dx = x - coord.x;
    if (dx > site_limit.x * 0.5){
        dx = dx - site_limit.x;
    }
    else if (dx <= -site_limit.x * 0.5){
        dx = dx + site_limit.x;
    }
    double dy = y - coord.y;
    if (dy > site_limit.y * 0.5){
        dy = dy - site_limit.y;
    }
    else if (dy <= -site_limit.y * 0.5){
        dy = dy + site_limit.y;
    }
    double dz = z - coord.z;
    if (dz > site_limit.z * 0.5){
        dz = dx - site_limit.z;
    }
    else if (dz <= -site_limit.z * 0.5){
        dz = dz + site_limit.z;
    }
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double Coordinate3D::getMagnitude() const {
    double a;
    if (x == 0 && y == 0 && z == 0){
        a = 0;
    }
    else{
        a = sqrt(x*x + y*y + z*z);
    }
    return a;
}
void Coordinate3D::setMagnitude(double length) {
    double magnitude = getMagnitude();
    if (magnitude > 0) {
      *this *= length/magnitude;
    }
}
Coordinate3D Coordinate3D::operator+(const Coordinate3D& vec) const {
    return {vec.x+x,vec.y +y, vec.z+z};
}
Coordinate3D& Coordinate3D::operator+=(const Coordinate3D& vec) {
    x += vec.x;
    y += vec.y;
    z += vec.z;
    return *this;
}
Coordinate3D Coordinate3D::operator-(const Coordinate3D& vec) const {
    return {x-vec.x,y-vec.y, z-vec.z};
}
Coordinate3D& Coordinate3D::operator-=(const Coordinate3D& vec) {
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
    return *this;
}
Coordinate3D Coordinate3D::operator*(double value) const {
    return {x*value,y*value, z*value};
}
Coordinate3D& Coordinate3D::operator*=(double value) {
    x *= value;
    y *= value;
    z *= value;
    return *this;
}
double Coordinate3D::scalarProduct(const Coordinate3D& vec) const {
    return x*vec.x + y*vec.y + z*vec.z;
}
Coordinate3D Coordinate3D::crossProduct(const Coordinate3D& vec) const {
    Coordinate3D new_vec{};
    new_vec.x = y*vec.z - z*vec.y;
    new_vec.y = z*vec.x  - x*vec.z ;
    new_vec.z = x*vec.y  - y*vec.x;
    return new_vec;

}
