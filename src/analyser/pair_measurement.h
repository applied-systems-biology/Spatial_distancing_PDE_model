//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef ABM_PAIR_MEASUREMENT_H_
#define ABM_PAIR_MEASUREMENT_H_

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <variant>
#include <utility>

#include "utils/macros.h"
#include "utils/misc_util.h"

class PairMeasurement {
 public:
    // Class for pair measurements. Pair measurements are measurements which are coupled with dynamic information in simulation such as the current time point
  using cache_type = std::variant<std::monostate, int, double, std::vector<int>, std::vector<double>>;
  template<typename... Args>
  explicit PairMeasurement(std::string id, Args &&... args)
      : owner_id_(std::move(id)), number_of_elements_(sizeof...(Args)), keys_({std::forward<Args>(args)...}) {}
  template<typename... Args>
  void addValuePairs(Args &&... args) {
    std::ostringstream output;
    ((output << std::forward<Args>(args) << delimeter), ...);
    data_.emplace_back(output.str());
  }
  template<typename T, typename... Args>
  typename std::enable_if<(std::is_convertible<Args, std::string>::value && ...), void>::type
  addValuesFromCache(T value, Args &&...args) {
    auto cache_access = abm::util::overloaded{
        [&](std::monostate /*unused*/)->std::string { return ""; }, [&](auto arg)->std::string { return std::to_string(arg); },
        [&](const std::vector< int> &vec)->std::string {
          std::ostringstream out;
          for (const auto &arg: vec) {
            out << arg << ",";
          }
          return out.str();
        }, [&](const std::vector<double> &vec)->std::string {
          std::ostringstream out;
          for (const auto &arg: vec) {
            out << arg << ",";
          }
          return out.str();
        }
    };
    auto cache_reset = abm::util::overloaded{
        [&](auto&  /*unused*/)->cache_type { return {}; },
        [&](std::vector<int>& /*unused*/)->cache_type { return std::vector<int>{}; },
        [&](std::vector<double>& /*unused*/)->cache_type { return std::vector<double>{}; },
    };
    std::ostringstream output;
    output << value << delimeter;
    ((output << std::visit(cache_access, cache[std::forward<Args>(args)]) << delimeter), ...);
    data_.emplace_back(output.str());
    ((std::visit(cache_reset, cache[std::forward<Args>(args)])), ...);
  }
  const std::vector<std::string> &getKeys() {
    return keys_;
  }
  const void clearData() {
      data_ = std::vector<std::string> {};
  }
  friend std::ostream &operator<<(std::ostream &out, PairMeasurement &measurement);
  static constexpr auto delimeter = ';';
  std::map<std::string, cache_type> cache{};
 private:
  std::string owner_id_{};
  int number_of_elements_;
  std::vector<std::string> keys_{};
  std::vector<std::string> data_;
};
std::ostream &operator<<(std::ostream &out, PairMeasurement &measurement);
#endif //ABM_PAIR_MEASUREMENT_H_
