//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef ALGORITHMS_H
#define	ALGORITHMS_H
#include <boost/numeric/ublas/matrix.hpp>

#include "basic/Randomizer.h"


class Algorithms {
public:

    
    /**
     *
     * @param size permutation size
     * @return a random permutation
     */
    static std::vector<unsigned int> generateRandomPermutation(Randomizer* randomizer, unsigned int size);
    
    /**
     * n über k berechnen
     * @param n
     * @param k
     * @return 
     */
    static double nChoosek( unsigned long n, unsigned long k );
    
    static double bernoulliProbability(unsigned long n, unsigned long k, double p);


    /*
     *  Matrix inversion routine.
     *  Uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
     */
    bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);

private:

    /**
     * swaps a pair of positions in a given vector
     * @param vectorToSwap vector in which the swap should take place
     * @param i position 1 in the array
     * @param j position 2 in the array
     */
    static void swap(std::vector<unsigned int>* vectorToSwap, unsigned int i, unsigned int j);

};

#endif	/* ALGORITHMS_H */

