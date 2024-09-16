//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "utils/misc_util.h"

#include "simulation/Agent.h"
#include "utils/macros.h"
#include <bits/basic_string.h>
#include <boost/algorithm/string.hpp>
#include <dirent.h>
#include <set>
#include <sstream>
#include <fstream>

namespace abm::util {

    std::string generateHashFromAgents(const double &current_time, const std::vector<std::shared_ptr<Agent>> &agents) {
        std::set<std::string> sorted_content;
        std::transform(agents.begin(),
                       agents.end(),
                       std::inserter(sorted_content, sorted_content.begin()),
                       [current_time](const auto &agent) {
                           return std::to_string(agent->getId()) + "-" +
                                  std::to_string(agent->getSurface()->getAllSpheresOfThis().size());// +
                                  //"-" + agent->getAgentCSVTagAoc(current_time);
                       });
        std::ostringstream buff_string;
        std::copy(sorted_content.begin(), sorted_content.end(), std::ostream_iterator<std::string>(buff_string, " "));
        return std::to_string(std::hash<std::string>{}(buff_string.str()));
    }

    std::string generateHashFromMolecules(const double &current_time, const std::map<std::string,std::vector<double>> &concentrations){
        std::vector<double> mol;
        for (const auto &[name, val]: concentrations){
            mol.insert(mol.end(), val.begin(), val.end());
        }
        std::set<std::string> sorted_content;
        std::transform(mol.begin(),
                       mol.end(),
                       std::inserter(sorted_content, sorted_content.begin()),
                       [current_time](const auto &conc) {
                           return std::to_string(conc);
                       });
        std::ostringstream buff_string;
        std::copy(sorted_content.begin(), sorted_content.end(), std::ostream_iterator<std::string>(buff_string, " "));
        return std::to_string(std::hash<std::string>{}(buff_string.str()));
    }


    template<char T, typename... Args>
    std::string concatenate(Args &&... args) {
        std::ostringstream output;
        ((output << std::forward<Args>(args) << T), ...);
        return output.str();
    }

    std::vector<std::string> getFileNamesFromDirectory(const std::string &path, const std::string &fileMask) {
        std::vector<std::string> fileList;
        DIR *dirp = opendir(path.c_str());
        struct dirent *dp;
        while ((dp = readdir(dirp)) != nullptr) {
            if (static_cast<std::string>(dp->d_name).find(fileMask) != std::string::npos) {
                fileList.emplace_back(dp->d_name);
            }
        }
        std::sort(fileList.begin(),
                  fileList.end()); // apply some order on files to have equal vector on different machines
        closedir(dirp);
        return fileList;
    }

    void read3DCoordinatesFromFile(std::vector<Coordinate3D> &AMpos, const std::string &inputString) {
        std::ostringstream am_dist_path;
        am_dist_path << inputString;
        std::ifstream file(am_dist_path.str().c_str());
        std::string currentLine;
        getline(file, currentLine);
        while (getline(file, currentLine)) {
            std::vector<std::string> tokens;
            boost::algorithm::split(tokens, currentLine, boost::is_any_of(","));
            double x = atof(tokens.at(1).c_str());
            double y = atof(tokens.at(2).c_str());
            double z = atof(tokens.at(3).c_str());

            Coordinate3D receivedCoordinate{x, y, z};
            AMpos.emplace_back(receivedCoordinate);
        }
        file.close();
    }

    double readLambdaValueFromFile(const std::string &inputString) {
        std::ostringstream AM_dist_path;
        AM_dist_path << inputString;
        std::ifstream file(AM_dist_path.str().c_str());
        std::string currentLine;
        getline(file, currentLine);
        return std::stod(currentLine);
    }

    std::unordered_map<std::string, std::string> handleCmdInputs(int argc, char **argv) {

        std::unordered_map<std::string, std::string> inputArgs;

        if (argc > 2) {
            if (argc % 2 == 0) {
                for (int i = 2; i < argc; i = i + 2) {
                    std::ostringstream ss, ssVal;
                    ss << argv[i];
                    std::string curArg = ss.str();
                    std::string firstCharArg = curArg.substr(0, 1);
                    if (firstCharArg == "-") {
                        ssVal << argv[i + 1];
                        std::string curVal = ssVal.str();
                        std::string withoutFirstCharArg = curArg.substr(1, -1);
                        inputArgs[withoutFirstCharArg] = curVal;
                    } else {
                        ERROR_STDERR(
                                R"(usage hint: "./Spatial_Distancing path/to/config.json -id1 value1 -id2 value2 ... ")");
                        ERROR_STDERR("now stopping execution");
                        exit(1);
                    }
                }
            } else {
                ERROR_STDERR(
                        R"(usage hint: "./Spatial_Distancing path/to/config.json -id1 value1 -id2 value2 ... ")");
                ERROR_STDERR("now stopping execution");
                exit(1);
            }
        }

        return inputArgs;
    }

    std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    getCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para) {
        std::vector<std::string> keys;
        std::vector<std::vector<std::string>> values;
        std::vector<int> its;
        std::vector<std::vector<std::string>> tuples;
        for (const auto &x: para) {
            keys.emplace_back(x.first);
            values.emplace_back(x.second);
            its.emplace_back(0);
        }

        bool keep_running = true;
        while (keep_running) {
            std::vector<std::string> valtuple;
            for (int i = 0; i < keys.size(); ++i) {
                valtuple.emplace_back(values[i][its[i]]);
            }
            tuples.emplace_back(valtuple);

            for (int j = 0; j < keys.size();) {
                if (its[j] < values[j].size() - 1) {
                    its[j] += 1;
                    j = keys.size();
                } else {
                    if (j == keys.size() - 1) {
                        keep_running = false;
                    }
                    its[j] = 0;
                    ++j;
                }
            }

        }
        return {keys, tuples};
    }

    std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    calculateCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para) {
        std::vector<std::string> keys;
        std::vector<std::vector<std::string>> values;
        std::vector<int> its;
        std::vector<std::vector<std::string>> tuples;
        for (const auto &x: para) {
            keys.emplace_back(x.first);
            values.emplace_back(x.second);
            its.emplace_back(0);
        }

        bool keep_running = true;
        while (keep_running) {
            std::vector<std::string> valtuple;
            for (int i = 0; i < keys.size(); ++i) {
                valtuple.emplace_back(values[i][its[i]]);
            }
            tuples.emplace_back(valtuple);

            for (int j = 0; j < keys.size();) {
                if (its[j] < values[j].size() - 1) {
                    its[j] += 1;
                    j = keys.size();
                } else {
                    if (j == keys.size() - 1) {
                        keep_running = false;
                    }
                    its[j] = 0;
                    ++j;
                }
            }

        }
        return {keys, tuples};
    }
}