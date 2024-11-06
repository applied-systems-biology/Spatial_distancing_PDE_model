//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include <boost/filesystem.hpp>
#include <iostream>
#include <sstream>
#include "simulation/Simulator.h"
#include "utils/io_util.h"
#include "utils/macros.h"
#include "utils/misc_util.h"
#include <omp.h>
#include "simulation/Site.h"

int main(int argc, char** argv) {

    // Default configuration location (if no parameter is specified)
    boost::filesystem::path config_xml("../../config.json");
    if (argc > 1 && !(std::istringstream(argv[1]) >> config_xml)) {
        ERROR_STDERR("usage: " << argv[0] << " <config.json>");
        return 1;
    }
    if (!(boost::filesystem::exists(config_xml))) {
        ERROR_STDERR("Configuration File in " << config_xml.string() << " does not exist!");
        ERROR_STDERR("usage: " << argv[0] << " <config.json>");
        return 2;
    }
    // change root directory for simulation to configuration location
    chdir(&config_xml.parent_path().c_str()[0]);
    const auto parameters = abm::util::getMainConfigParameters(config_xml.filename().string());
    std::unordered_map<std::string, std::string> input_args = abm::util::handleCmdInputs(argc, argv);

#if defined(_OPENMP)
    omp_set_num_threads(parameters.number_of_threads);
    DEBUG_STDOUT("OpenMP activated with " << parameters.number_of_threads << " Thread(s).");
#endif

    if (parameters.screening_parameters.empty()) {
        const auto simulator = std::make_unique<Simulator>();
        simulator->setConfigPath(parameters.config_path);
        simulator->setCmdInputArgs(input_args);
        simulator->setOutputPath(parameters.output_dir);
        simulator.get()->executeRuns(parameters.runs, parameters.system_seed, parameters.output_dir,
                                     parameters.input_dir);
    } else {
        /// Screening over all parameter combinations specified as sets in the configuration file <config.json>
        // For screening, the cartesian product of all the single parameters sets is generated
        // If you want to resume a screening from a certain index, you can change screen_start_idx in main config
        const auto [parameter_names, value_combinations] = abm::util::calculateCartesianProd(
                parameters.screening_parameters);
        std::vector<std::unordered_map<std::string, std::string>> local_input_args_list(
                parameters.number_of_threads); // Vector to store intermediate results

#pragma omp parallel for schedule(dynamic)
        for (int sim = parameters.screen_start_idx - 1; sim < value_combinations.size(); ++sim) {
            int thread_id = omp_get_thread_num();
            std::unordered_map<std::string, std::string> local_input_args; // Create a local copy of input_args
            std::stringstream sim_para{};

            for (int i = 0; i < parameter_names.size(); ++i) {
                sim_para << parameter_names[i] << value_combinations[sim][i] << "_";
                local_input_args[parameter_names[i]] = value_combinations[sim][i]; // Update the local copy
            }

            local_input_args_list[thread_id] = std::move(
                    local_input_args); // Store the local copy in the thread-specific vector

            SYSTEM_STDOUT("Start " << parameters.runs << " runs of simulation " << sim + 1 << "/"
                                   << value_combinations.size());
            auto simulator = std::make_unique<Simulator>(parameters.config_path,
                                                         local_input_args_list[thread_id]); // Use the local copy
            simulator->setConfigPath(parameters.config_path);
            simulator->setCmdInputArgs(local_input_args_list[thread_id]);
            simulator->setOutputPath(parameters.output_dir);
            simulator->executeRuns(parameters.runs, parameters.system_seed, parameters.output_dir, parameters.input_dir,
                                   sim, sim_para.str());
        }
    }
    return 0;
}