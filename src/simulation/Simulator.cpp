//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "simulation/Simulator.h"
#include <chrono>
#include <omp.h>
#include <string>

#include "io/output_handler.h"
#include "simulation/CellFactory.h"
#include "utils/macros.h"
#include "utils/time_util.h"

#include "simulation/Site.h"
#include "simulation/site/CuboidSite.h"

#include "analyser/Analyser.h"



Simulator::Simulator(std::string config_path, std::unordered_map<std::string, std::string> cmd_input_args)
        : config_path_(config_path), cmd_input_args_(cmd_input_args) {
}

void Simulator::executeRuns(int runs, int seed, const std::string& output_dir, const std::string& input_dir, int sim, const std::string& parameter_string) const {
    // Initializes the seed for the current parameter configuration
    // Take system seed from command line IF provided
    seed = cmd_input_args_.find("seed") != cmd_input_args_.end() ? std::stoi(cmd_input_args_.find("seed")->second) : seed;
    SYSTEM_STDOUT("System Seed: " << seed);
    int const current_sim_seed = seed + runs * sim;
    SYSTEM_STDOUT("Current Simulation Seed: " << current_sim_seed);
    // Initializes project name and output directory
    std::ostringstream project_name;
    project_name << abm::util::getCurrentLocalTimeAsString() << "_" << parameter_string << current_sim_seed;
    auto project_name_str = project_name.str();
    int const max_filename_size = 254;
    if (project_name_str.length() > max_filename_size) {
        project_name_str = project_name_str.substr(0, max_filename_size);
    }

    // initializes session key sid for different result folders, e.g. ./coreABM ../../config.json -sid abc123
    auto sid = cmd_input_args_.find("sid") != cmd_input_args_.end() ? cmd_input_args_.find("sid")->second : "";
    std::string const output_folder_name = "results" + sid;
    const auto project_dir = static_cast<boost::filesystem::path>(output_dir).append(output_folder_name).append(project_name_str).string();


    const auto output_handler = std::make_unique<const OutputHandler>(config_path_, project_dir);
    const auto analyser = std::make_unique<const Analyser>(config_path_, project_dir);

    using timer = std::chrono::steady_clock;
    const auto start = timer::now();
#pragma omp parallel for schedule(dynamic)
    for (int current_run = 1; current_run <= runs; ++current_run) {
        SYSTEM_STDOUT("Thread " << omp_get_thread_num() << ": Start Run " << current_run << "/" << runs);
        if (cmd_input_args_.size() > 0) {
            if (current_run == 1) {
                const auto current_screening_csv = static_cast<boost::filesystem::path>(project_dir).append("screening.csv").string();
                std::ofstream csv_file(current_screening_csv, std::ofstream::out | std::ofstream::app);
                csv_file << "seed,";
                for (const auto&[key, value] : cmd_input_args_) {
                    csv_file << key << ",";
                }

                csv_file << "\n" << seed << ",";
                for (const auto&[key, value] : cmd_input_args_) {
                    csv_file << value << ",";
                }
                csv_file.close();
            }
        }
        int run_seed = current_run + current_sim_seed;
        const auto random_generator = std::make_unique<Randomizer>(run_seed);
        const auto site = createSites(current_run, random_generator.get(), analyser.get());
        SimulationTime time{site->getTimeStepping(), site->getMaxTime()}; time.updateTimestep(0);
        for (time.updateTimestep(0); !time.endReached(); ++time) { //inner loop: t -> t + dt
            output_handler->outputCurrentConfiguration(*site, time, current_run);
            site->doAgentDynamics(random_generator.get(), time);

            if (time.checkForNumberOfExecutions(11, true)) {
                SYSTEM_STDOUT("Run: " << current_run << ", Time: " << time.getCurrentTime());
            }
            if (time.getCurrentTime() == time.getMaxTime()*0.1 || time.getCurrentTime() == time.getMaxTime()*0.25 ||time.getCurrentTime() == time.getMaxTime()*0.5 ||time.getCurrentTime() == time.getMaxTime()*0.75){ //save results at 10%, 25, 50 and 75% of the simulation
                analyser->outputAllMeasurements();
                auto grid = site->molecule_manager_->getGridSize();
                auto mol_conc = site->molecule_manager_->get_conc().at("AMP");
                for (auto conc: mol_conc){
                    if (conc == NAN || conc < 0 || conc > 1e100){
                        SYSTEM_STDOUT("Invalid molecule concentration - likely due to timestep too large");
                        break;
                    }
                }
            }
            if(site->get_steady_state()){
                SYSTEM_STDOUT("Stopping criteria met:");
                SYSTEM_STDOUT("\t---> Simulation stopped at t=" << time.getCurrentTime());
                break;
            }
        }
        output_handler->outputCurrentConfiguration(*site, time, current_run, run_seed, true, cmd_input_args_);
        site->getMeasurments()->post_simulations_measurements(time);
    }
    const auto end = timer::now();
    analyser->outputAllMeasurements();
    output_handler->concludeSimulation(runs,
                                       seed,
                                       static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()));
}

std::unique_ptr<Site> Simulator::createSites(int run, Randomizer* random_generator, const Analyser* analyser) const {
    return std::make_unique<CuboidSite>(random_generator,
                                            analyser->generateMeasurement(std::to_string(run)),
                                            config_path_,
                                            cmd_input_args_,
                                            output_dir_);
}

