//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "utils/time_util.h"

#include "utils/macros.h"
#include <algorithm>
#include <cassert>
#include <iomanip>

namespace abm::util {
std::string getCurrentLocalTimeAsString() {
    std::ostringstream output;
    auto t = std::time(nullptr);
    output << std::put_time(std::localtime(&t), "%F_%OH-%OM-%OS");
    return output.str();
}

std::string generateUniformString(double time, const int fill_to_size) {
    std::string ret_string;
    constexpr auto precision = 6;
    std::ostringstream buffer;
    buffer << std::right << std::setfill('0') << std::setw(fill_to_size) << std::setprecision(precision) << std::fixed
           << time;
    ret_string = buffer.str();
    std::replace(ret_string.begin(), ret_string.end(), '.', 'k');
    return ret_string;
}
} // namespace abm::util

SimulationTime& SimulationTime::operator++() noexcept {
    ++time_step_;
    current_time_ = time_step_ * delta_t_;
    return *this;
}
SimulationTime& SimulationTime::operator--() noexcept {
    --time_step_;
    current_time_ = time_step_ * delta_t_;
    return *this;
}

void SimulationTime::updateDeltaT(double new_delta) noexcept {
    assert(new_delta != 0);
    last_delta_t_ = delta_t_;
    delta_t_ = new_delta;
    time_step_ = static_cast<int>(std::ceil(time_step_ * (last_delta_t_ / delta_t_)));
    max_steps_ = static_cast<int>(std::ceil(max_time_ / delta_t_));
}
double SimulationTime::getMaxTime() const noexcept {
    return max_time_;
}
double SimulationTime::getCurrentDeltaT() const noexcept {
    return delta_t_;
}
double SimulationTime::getLastDeltaT() const noexcept {
    return last_delta_t_;
}
double SimulationTime::getCurrentTime() const noexcept {
    return current_time_;
}
int SimulationTime::getCurrentTimeStep() const noexcept {
    return time_step_;
}
bool SimulationTime::endReached() const noexcept {
    return time_step_ >= max_steps_;
}
bool SimulationTime::checkForNumberOfExecutions(int number_of_executions, bool execute_first) const {
    if (execute_first) {
        if (0 == time_step_) {
            return true;
        }
        if (1 == number_of_executions) {
            return false;
        }
        number_of_executions -= 1;
    }

    assert(number_of_executions > 0);
    const auto step_size = static_cast<int>(floor(static_cast<double>(max_steps_) / number_of_executions));
    assert(step_size > 0);
    return 0 == (time_step_ + 1) % step_size; //time step needs an increment here, since first step is 0
}

int operator%(const SimulationTime& time, int value) {
    assert(value != 0);
    return time.time_step_ % value;
}
void SimulationTime::updateTimestep(int new_step) noexcept {
    time_step_ = new_step;
    current_time_ = delta_t_ * time_step_;
}
