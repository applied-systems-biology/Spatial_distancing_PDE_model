//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#ifndef ABM_UTILS_MISC_UTIL_H_
#define ABM_UTILS_MISC_UTIL_H_

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <basic/Coordinate3D.h>
#include <map>
#include "simulation/MoleculeManager.h"

class Agent;
class Randomizer;
namespace abm::util {

//little helper for std::visit
template<typename... Ts>
struct overloaded : Ts ... { using Ts::operator()...; };
template<typename... Ts> overloaded(Ts...) -> overloaded<Ts...>;

std::string generateHashFromAgents(const double &current_time, const std::vector<std::shared_ptr<Agent>> &agents);
std::string generateHashFromMolecules(const double &current_time, const std::map<std::string, std::vector<double>> &concentrations);
template<char T = ';', typename... Args>
std::string concatenate(Args &&... args);

std::vector<std::string> getFileNamesFromDirectory(const std::string& path, const std::string& fileMask="");

std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    calculateCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para);



    void read3DCoordinatesFromFile(std::vector<Coordinate3D>& AMpos, const std::string& inputString);

double readLambdaValueFromFile(const std::string& input_string);

std::unordered_map<std::string, std::string> handleCmdInputs(int argc, char **argv);

std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>> getCartesianProd(const std::unordered_map<std::string, std::vector<std::string>>& para);

bool approxEqual(double d1, double d2, double epsilon=1e-8);


std::pair<double, double> getRandomDirection(const Coordinate3D & coordinate, Randomizer* randomizer, double mean, double std, std::pair<double, double> prevDir, std::string origin="");
std::pair<double, double> getRandomDirectionOfLength(const Coordinate3D & coordinate, Randomizer* randomizer, std::string origin="");
}
#endif //ABM_UTILS_MISC_UTIL_H_
