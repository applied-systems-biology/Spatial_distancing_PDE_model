//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include <vector>
#include <cmath>
#include <algorithm>

#include <stdio.h>
#include <math.h>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <utils/misc_util.h>

#include "simulation/Algorithms.h"
#include "utils/macros.h"

std::vector<unsigned int> Algorithms::generateRandomPermutation(Randomizer* randomizer, unsigned int size){
  std::vector<unsigned int> permutationVector;
  if (size > 0) {
    //initialize the vector
    for (unsigned int k = 0; k < size; k++) {
      permutationVector.push_back(k);
    }


    for (unsigned int k = 0; k < size - 1; k++) {
      swap(&permutationVector, k, randomizer->generateInt(k, size - 1));
    }
  }
  return permutationVector;
}

void Algorithms::swap(std::vector<unsigned int>* vectorToSwap, unsigned int i, unsigned int j){
  unsigned int temp;
  temp = (*vectorToSwap)[i];
  (*vectorToSwap)[i] = (*vectorToSwap)[j];
  (*vectorToSwap)[j] = temp;
}

double Algorithms::nChoosek( unsigned long n, unsigned long k )
{
  if (k * 2 > n){ k = n-k;}

  if (k > n) return 0;
  if (k == 0) return 1;


  double result = 1.0;
  for(  unsigned long i = 1; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

double Algorithms::bernoulliProbability(unsigned long n, unsigned long k, double p){

  return nChoosek(n,k)*pow(p,k)*pow(1-p,n-k);
}


using namespace boost::numeric::ublas;
bool Algorithms::InvertMatrix(const matrix<double>& input, matrix<double>& inverse)
{
  typedef permutation_matrix<std::size_t> pmatrix;

  // create a working copy of the input
  matrix<double> A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(identity_matrix<double> (A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}