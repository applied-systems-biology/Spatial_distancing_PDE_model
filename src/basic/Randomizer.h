//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.



#ifndef RANDOMIZER_H
#define	RANDOMIZER_H

#include <cstdlib>
#include <iostream>

#include "basic/Coordinate3D.h"

#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>


class Randomizer {
public:
    // Class for wrapping C++ random generator
    explicit Randomizer(int seed);
    double generateDouble();
    double generateDouble(double);
    double generateDouble(double, double);
    unsigned int generateInt(unsigned int maxVal);
    unsigned int generateInt(unsigned int minVal, unsigned int maxVal);

    Coordinate3D generateRandomDirection(unsigned int spatialDims, double length);
    double generateReighlayDistributedValue(double sigma = 1);
    double generateNormalDistributedValue(double mean = 0, double stddev = 1, bool box_muller_method = true);

private:
    int seed_{};
    boost::mt19937 random_mt_{};
    bool second_gauss_available_{false};
    double second_gauss_value_{0.0};
    double getNormalDistributed01Value(bool box_muller_method = true);
};
#endif	/* RANDOMIZER_H */