//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "basic/Randomizer.h"


#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

Randomizer::Randomizer(int seed) : seed_(seed), random_mt_(seed) {}

double Randomizer::generateDouble() {

  return 1.0 * (random_mt_)() / boost::mt19937::max();
}

double Randomizer::generateDouble(double max_val) {

  return (1.0 * (random_mt_)() / boost::mt19937::max()) * max_val;
}

double Randomizer::generateDouble(double min_val, double max_val) {

  return min_val + (1.0 * (random_mt_)() / boost::mt19937::max()) * (max_val - min_val);
}

unsigned int Randomizer::generateInt(unsigned int maxVal){
  return (unsigned int) ((random_mt_)() % (maxVal + 1));
}

unsigned int Randomizer::generateInt(unsigned int minVal, unsigned int maxVal){
  return (unsigned int) minVal + ((random_mt_)() % ((maxVal - minVal) + 1));
}

Coordinate3D Randomizer::generateRandomDirection(unsigned int spatialDims, double length){
  double u,phi,r,x,y,z,subst;

  switch(spatialDims){
    case 2:
      phi = generateDouble(M_PI*2.0);
      r = length;

      x = r*cos(phi);
      y = r*sin(phi);
      z = 0.0;
      break;

    case 3:
      //get a sin() sampled value in [0,PI] for theta via a uniform value of u
      u = generateDouble();//sampler->sample();

      phi = generateDouble(M_PI*2.0);
      r = length;

      subst = 2*r*sqrt(u*(1-u));
      x = subst*cos(phi);
      y = subst*sin(phi);
      z = r*(1-2*u);

      break;


    default:
      //get a sin() sampled value in [0,PI] for theta via a uniform value of u
      u = generateDouble();//sampler->sample();

      phi = generateDouble(M_PI*2.0);
      r = length;

      subst = 2*r*sqrt(u*(1-u));
      x = subst*cos(phi);
      y = subst*sin(phi);
      z = r*(1-2*u);


      break;
  }
  return Coordinate3D{x,y,z};
}


double Randomizer::generateReighlayDistributedValue(double sigma) {
  double normX, normY;
  normX = generateNormalDistributedValue(0, sigma);
  normY = generateNormalDistributedValue(0, sigma);

  return sqrt(normX * normX + normY * normY);
}

double Randomizer::generateNormalDistributedValue(double mean, double stddev, bool box_muller_method){
  double returnVal;
  if(stddev == 0){
    returnVal = mean;
  }else{
    returnVal = mean + stddev * getNormalDistributed01Value(box_muller_method);
  }
  return returnVal;
}

double Randomizer::getNormalDistributed01Value(bool box_muller_method) {
  double value;
  double u1, u2;

  if (second_gauss_available_) {
    value = second_gauss_value_;
    second_gauss_available_ = false;
  } else {
    double interim;
    if (box_muller_method) {
      u1 = 1 - generateDouble();
      u2 = generateDouble();
      interim = sqrt(-2.0 * log(u1));
      value = interim * cos(2 * M_PI * u2);
      second_gauss_value_ = interim * sin(2 * M_PI * u2);
    } else {
      //application of polar method
      double r;
      do {
        u1 = 2 * generateDouble() - 1;
        u2 = 2 * generateDouble() - 1;
        r = u1 * u1 + u2 * u2;
      } while (r > 1 || r == 0);
      interim = sqrt(-2.0 * log(r) / r);
      value = interim * u1;
      second_gauss_value_ = interim * u2;
    }
    second_gauss_available_ = true;
  }

  return value;
}