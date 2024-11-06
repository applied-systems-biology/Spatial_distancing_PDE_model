//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

//
// Created by prudolph on 10.08.20.
//


#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "external/doctest/doctest.h"

#include <memory>
#include <unordered_set>

#include "analyser/pair_measurement.h"
#include "analyser/histogram_measurment.h"


TEST_CASE ("Check Pair Measurements") {
  const auto measurement = std::make_unique<PairMeasurement>("TEST", "Case1", "Case2");
  measurement->addValuePairs(1, 2);
  std::ostringstream out;
  out << *measurement;
  CHECK(out.str() == "TEST;1;2;\n");
  measurement->cache["Case2"] = 5;
  measurement->addValuesFromCache(1, "Case2");
  out.str("");
  out << *measurement;
  CHECK(out.str() == "TEST;1;2;\nTEST;1;5;\n");
  measurement->cache["Case2"] = std::vector<int>{5,2};
  measurement->addValuesFromCache(3, "Case2");
  out.str("");
  out << *measurement;
  CHECK(out.str() == "TEST;1;2;\nTEST;1;5;\nTEST;3;5,2,;\n");
}

TEST_CASE ("Check Histogram Measurements") {
  const auto measurement = std::make_unique<HistogramMeasurement>("TEST", "Case1", "Case2");
  measurement->addValues(std::make_pair("Case1",1), std::make_pair("Case2",2));
  measurement->addValues(std::make_pair("Case2",5), std::make_pair("Case2",2));
  std::ostringstream out;
  out << *measurement;
  CHECK(out.str() == "TEST;1;2;\nTEST;;5;\nTEST;;2;\n");
  out.str("");
  measurement->cache["Case1"] = 2;
  measurement->addValuesFromCache("Case1");
  out << *measurement;
  CHECK(out.str() == "TEST;1;2;\nTEST;2;5;\nTEST;;2;\n");
}

TEST_CASE("Check distance function") {
  Coordinate3D test1{1.0, 2.0, 4.0};
  Coordinate3D test2{4.0, 6.0, 4.0};
  double result = 5.0;
          CHECK(test1.calculateEuclidianDistance(test2) == result);
}