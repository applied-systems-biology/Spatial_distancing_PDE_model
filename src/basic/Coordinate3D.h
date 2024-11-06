//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef COORDINATE3D_H
#define	COORDINATE3D_H


struct Coordinate3D {
    double x;
    double y;
    double z;

    Coordinate3D operator+(const Coordinate3D &vec) const;
    Coordinate3D &operator+=(const Coordinate3D &vec);
    Coordinate3D operator-(const Coordinate3D &vec) const;
    Coordinate3D &operator-=(const Coordinate3D &vec);
    Coordinate3D operator*(double value) const;
    Coordinate3D &operator*=(double value);

    [[nodiscard]] double getMagnitude() const;
    void setMagnitude(double length);
    [[nodiscard]] double scalarProduct(const Coordinate3D& vec) const;
    [[nodiscard]] double calculateEuclidianDistance(const Coordinate3D& coord) const noexcept;
    [[nodiscard]] Coordinate3D crossProduct(const Coordinate3D& vec) const;
    double calculateEuclidianDistancePeriodic(const Coordinate3D& coord, Coordinate3D& site_limit) const noexcept;
};

#endif	/* COORDINATE3D_H */


