//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "output_handler.h"

#include <iostream>
#include <utility>
#include <boost/range/iterator_range.hpp>

#include "utils/time_util.h"
#include "utils/io_util.h"
#include "utils/macros.h"
#include "simulation/Site.h"
#include "utils/misc_util.h"

using boost::filesystem::path;
using boost::filesystem::exists;

OutputHandler::OutputHandler(const std::string &config_path, std::string project_dir)
    : current_project_string_(static_cast<path>(project_dir).filename().string()), current_project_dir_(std::move(project_dir)) {
  parameters_ = abm::util::getOutputParameters(static_cast<path>(config_path).append("output-config.json").string());
  config_param_ = abm::util::getMainConfigParameters(static_cast<path>(config_path).append("config.json").string());
  try {
      const auto copyOptions = boost::filesystem::copy_option::overwrite_if_exists;
    boost::filesystem::create_directories(current_project_dir_);
    //copy configuration files ending with "-config.xml"
    for (const auto &entry : boost::make_iterator_range(boost::filesystem::directory_iterator(config_path), {})) {
      if (const auto file_name = entry.path().filename().string(); file_name.find("-config.json") != std::string::npos) {
        boost::filesystem::copy_file(entry, static_cast<path>(current_project_dir_).append(file_name),copyOptions);
      }
    }

//    boost::filesystem::directory_entry("../../config.json")
    boost::filesystem::path config_file{"config.json"};
    boost::filesystem::copy_file(config_file, static_cast<path>(current_project_dir_).append("config.json"), copyOptions);

  } catch (const std::exception &e) {
    throw;
  }

}

void OutputHandler::concludeSimulation(int runs, int seed, double runtime) const {

  std::ifstream buffStream{static_cast<path>(current_project_dir_).append("runs.csv").string()};
  std::vector<std::string> hashes;
  for (std::string line; std::getline(buffStream, line);) {
    hashes.emplace_back(line.substr(line.find_last_of(';') + 1, line.length()));
  }
  std::sort(hashes.begin(), hashes.end());
  std::ostringstream buff_string;
  std::copy(hashes.begin(), hashes.end(), std::ostream_iterator<std::string>(buff_string, " "));

  const auto hash_of_runs = std::hash<std::string>{}(buff_string.str());
  DEBUG_STDOUT("Hash: " << hash_of_runs);

  if (parameters_.output_graphs){
      if (config_param_.conda_dir == ""){
          DEBUG_STDOUT("Please fill a value for conda_dir in output-config.json");
      } else {
          // Generate HTML files containing plots
          std::string command = config_param_.conda_dir; // conda folder from json?
          command.append(" run -n abmenv python ");
          command.append("../../python_scripts/output_graphs.py ");
          command.append(current_project_dir_);
          command.append("/measurements");
          DEBUG_STDOUT(command);
          system(command.c_str());
      }
  }

  if (parameters_.output_ovito){
      std::string command = config_param_.conda_dir;
      command.append(" run -n abmenv python ");
      command.append("../../python_scripts/ovito_preprocessing.py ");
      command.append(current_project_dir_);
      command.append("/measurements");
      system(command.c_str());
  }
}

void OutputHandler::concludeSimulationRun(int current_run, int seed, const std::string &hash, std::unordered_map<std::string, std::string> cmd_input_args) const {
  SYSTEM_STDOUT("Run: " << current_run << ", Hash: " << hash);
  const auto current_run_string = current_project_string_ + "-configs-" + std::to_string(current_run);
  const auto current_run_tar = static_cast<path>(current_project_dir_).append(current_run_string + ".tar.gz");
  const auto current_runs_csv = static_cast<path>(current_project_dir_).append("runs.csv").string();
  const auto current_run_dir = static_cast<path>(current_project_dir_).append(current_run_string);
  std::ofstream csv_file(current_runs_csv, std::ofstream::out | std::ofstream::app);
  csv_file << current_run_string << ";" << seed << ';' << hash << '\n';
  csv_file.close();
  if (exists(current_run_dir)) {
    std::ostringstream cmd_compression;
    try {
      cmd_compression << "tar czf " << current_run_tar << " -C " << current_project_dir_ << " " << current_run_string
                      << " --remove-files &";
      abm::util::executeShellCommand(cmd_compression.str());
    } catch (const std::exception &e) {
      throw;
    }
  }
}
void OutputHandler::setupOutputForRun(int run) const {
  if (parameters_.activated) {
    std::ostringstream run_directory;
    run_directory << current_project_string_ << "-configs-" << run;
  }
}

void OutputHandler::outputCurrentConfiguration(const Site &site,
                                               const SimulationTime &time,
                                               int run,
                                               int seed,
                                               bool simulation_end,
                                               std::unordered_map<std::string, std::string> cmd_input_args) const {

  if (simulation_end && !current_project_dir_.empty()) {
    concludeSimulationRun(run, seed, abm::util::generateHashFromMolecules(time.getCurrentTime(), site.molecule_manager_->get_conc()), cmd_input_args);
  }
}