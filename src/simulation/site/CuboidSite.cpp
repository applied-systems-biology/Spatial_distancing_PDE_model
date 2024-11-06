//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#include "CuboidSite.h"

#include "simulation/AgentManager.h"
#include "simulation/CellFactory.h"

CuboidSite::CuboidSite(Randomizer* random_generator,
                       std::shared_ptr<InSituMeasurements> measurements,
                       std::string config_path,
                       std::unordered_map<std::string, std::string> cmd_input_args,
                       std::string output_path)
        : Site(random_generator, std::move(measurements), config_path,  cmd_input_args, output_path) {


    // Read in simulator parameters
    auto sid = cmd_input_args.find("sid") != cmd_input_args.end() ? cmd_input_args.find("sid")->second : "";
    const auto simulator_config = static_cast<boost::filesystem::path>(config_path).append("simulator-config.json");
    parameters_ = abm::util::getSimulationParameters(simulator_config.string());


    // Read in cmd inputs and screening parameters
    parameters_.cmd_input_args = cmd_input_args;
    handleCmdInputArgs(cmd_input_args);

    cell_factory_ = std::make_unique<CellFactory>(parameters_.site_parameters);
    agent_manager_ = std::make_unique<AgentManager>(parameters_.time_stepping, this);

    auto* cuboid_parameters = static_cast<abm::util::SimulationParameters::CuboidSiteParameters*>(parameters_.site_parameters.get());
    identifier_ = cuboid_parameters->identifier;
    upper_bound_ = cuboid_parameters->upper_bound;
    lower_bound_ = cuboid_parameters->lower_bound;
    molecules_grid_size_ = cuboid_parameters->molecules_grid_size;
    molecular_layer_ = cuboid_parameters->molecular_layer_;

    molecule_manager_ = std::make_unique<MoleculeManager>(MoleculeManager(cuboid_parameters->molecule_manager_parameters, molecules_grid_size_, upper_bound_, lower_bound_));

    initializeAgents(cuboid_parameters->agent_manager_parameters, 0, parameters_.time_stepping);

    for (auto agent: agent_manager_->getAllAgents()) {
        molecule_manager_->initial_remove_mol_inside_cell(agent->get_inside());
    }

    steady_state_ = false;
    stopping_criteria_ = cuboid_parameters->stopping_criteria;

    //Automated time stepping computation (for uniform grid)
    if (parameters_.dt_auto) {
        auto D = 1.0;

        for (auto &name: molecule_manager_->get_molecules_names()) {
            if (D < std::get<0>(molecule_manager_->getParameters(name))) {
                D = std::get<0>(molecule_manager_->getParameters(name));
            }
        }
        auto cfl = std::pow(std::get<0>(molecule_manager_->getDistances()), 2) / (6 * D);
        auto max_rate = 1.0;
        for (auto &rate: parameters_.site_parameters->molecule_manager_parameters.molecule_molecule_interactions) {
            if (max_rate < rate.second.binding) {
                max_rate = rate.second.binding;
            }
            if (max_rate < rate.second.unbinding) {
                max_rate = rate.second.unbinding;
            }
        }
        for (auto &ag: parameters_.site_parameters->agent_manager_parameters.agents) {
            for (auto &uptake: ag->molecule_interactions) {
                if (max_rate < uptake.second.uptake) {
                    max_rate = uptake.second.uptake;
                }
            }
        }
        parameters_.time_stepping = (1 / ceil(1 / std::min(cfl, 1 / max_rate)))/10;
    }
}

bool CuboidSite::containsPosition(Coordinate3D position) {
    bool contains = (position.x >= upper_bound_.x) && (position.x <= lower_bound_.x) &&
                   (position.y >= upper_bound_.y) && (position.y <= lower_bound_.y) &&
                   (position.z >= upper_bound_.z) && (position.z <= lower_bound_.z);

    return contains;
}

Coordinate3D CuboidSite::getRandomPosition(double radius) {

    double   x = random_generator_->generateDouble(upper_bound_.x + radius, lower_bound_.x - radius);
    double   y = random_generator_->generateDouble(upper_bound_.y + radius, lower_bound_.y - radius);
    double   z = random_generator_->generateDouble(upper_bound_.z + radius, lower_bound_.z - radius);
    return Coordinate3D{x, y, z};
}

Coordinate3D CuboidSite::getCenterPosition() {

    double   x = (upper_bound_.x + lower_bound_.x) / 2.0;
    double   y = (upper_bound_.y + lower_bound_.y) / 2.0;
    double   z = (upper_bound_.z + lower_bound_.z) / 2.0;

    return Coordinate3D{x, y, z};
}

Coordinate3D CuboidSite::getLowerLimits() {
    return upper_bound_;
}

Coordinate3D CuboidSite::getUpperLimits() {
    return lower_bound_;
}

std::pair<std::vector<std::tuple<int, int, int>>, std::vector<std::tuple<int, int, int>>> CuboidSite::cell_grid_info(double radius, Coordinate3D position){

    const auto&[n_x, n_y, n_z] = molecule_manager_->getGridSize();
    const auto h_x = std::get<0>(molecule_manager_->getDistances());
    const auto h_y = std::get<1>(molecule_manager_->getDistances());
    const auto h_z = std::get<2>(molecule_manager_->getDistances());

    const auto normalized_position = position - upper_bound_; // was lower_bound_-position; // upper_bound_ = (-500, -500, -500)

    int x_min = std::max(0, static_cast<int>(std::floor((normalized_position.x-radius)/h_x)));
    int y_min = std::max(0, static_cast<int>(std::floor((normalized_position.y-radius)/h_y)));
    int z_min = std::max(0, static_cast<int>(std::floor((normalized_position.z-radius)/h_z)));
    int x_max = std::min(n_x-1, static_cast<int>(std::ceil((normalized_position.x+radius)/h_x)));
    int y_max = std::min(n_y-1, static_cast<int>(std::ceil((normalized_position.y+radius)/h_y)));
    int z_max = std::min(n_z-1, static_cast<int>(std::ceil((normalized_position.z+radius)/h_z)));

    std::vector<std::tuple<int, int, int>> covered_points;
    std::vector<std::tuple<int, int, int>> surface_points;

    for (int i = z_min; i <= z_max; ++i) {
        for (int j = y_min; j <= y_max; ++j) {
            for (int k = x_min; k <= x_max; ++k) {
                auto distance = position.calculateEuclidianDistance(molecule_manager_->get_loc(i, j, k));
                //SYSTEM_STDOUT("distance = " << distance);
                if (distance < radius-h_x/2) {
                    covered_points.emplace_back(i, j, k);
                }
                else if ((distance <= radius + h_x/2) && (distance >= radius - h_x/2)) {
                    surface_points.emplace_back(i, j, k);
                }
            }
        }
    }
    return std::make_pair(covered_points, surface_points);
}

std::vector<std::tuple<int, int, int>> CuboidSite::covered_points(double radius, Coordinate3D position){

    const auto&[n_x, n_y, n_z] = molecule_manager_->getGridSize();
    const auto h_x = std::get<0>(molecule_manager_->getDistances());
    const auto h_y = std::get<1>(molecule_manager_->getDistances());
    const auto h_z = std::get<2>(molecule_manager_->getDistances());

    const auto normalized_position = position - upper_bound_; // was lower_bound_-position; // upper_bound_ = (-500, -500, -500)

    int x_min = static_cast<int>(std::floor((normalized_position.x-radius)/h_x));
    int y_min = static_cast<int>(std::floor((normalized_position.y-radius)/h_y));
    int z_min = static_cast<int>(std::floor((normalized_position.z-radius)/h_z));
    int x_max = static_cast<int>(std::ceil((normalized_position.x+radius)/h_x));
    int y_max = static_cast<int>(std::ceil((normalized_position.y+radius)/h_y));
    int z_max = static_cast<int>(std::ceil((normalized_position.z+radius)/h_z));

    std::vector<std::tuple<int, int, int>> covered_points;

    for (int i = x_min; i <= x_max; ++i) {
        for (int j = y_min; j <= y_max; ++j) {
            for (int k = z_min; k <= z_max; ++k) {
                auto distance = position.calculateEuclidianDistance(molecule_manager_->get_loc(i, j, k));
                if (distance < radius-h_x/2) {
                    covered_points.emplace_back(i, j, k);
                }
            }
        }
    };
    return covered_points;
}

std::vector<std::tuple<int, int, int>> CuboidSite::surface_points(double radius, Coordinate3D position){
    const auto&[n_x, n_y, n_z] = molecule_manager_->getGridSize();
    const auto h_x = std::get<0>(molecule_manager_->getDistances());
    const auto h_y = std::get<1>(molecule_manager_->getDistances());
    const auto h_z = std::get<2>(molecule_manager_->getDistances());

    const auto normalized_position = position - upper_bound_; // was lower_bound_-position; // upper_bound_ = (-500, -500, -500)

    int x_min = static_cast<int>(std::floor((normalized_position.x-radius)/h_x));
    int y_min = static_cast<int>(std::floor((normalized_position.y-radius)/h_y));
    int z_min = static_cast<int>(std::floor((normalized_position.z-radius)/h_z));
    int x_max = static_cast<int>(std::ceil((normalized_position.x+radius)/h_x));
    int y_max = static_cast<int>(std::ceil((normalized_position.y+radius)/h_y));
    int z_max = static_cast<int>(std::ceil((normalized_position.z+radius)/h_z));

    std::vector<std::tuple<int, int, int>> surface_points;

    for (int i = x_min; i <= x_max; ++i) {
        for (int j = y_min; j <= y_max; ++j) {
            for (int k = z_min; k <= z_max; ++k) {
                auto distance = position.calculateEuclidianDistance(molecule_manager_->get_loc(i, j, k));
                if ((distance < radius + h_x/2) && (distance > radius - h_x/2)) {
                    surface_points.emplace_back(i, j, k);
                }
            }
        }
    }
    return surface_points;
}

float CuboidSite::average_conc_in_agent(const std::string& name, double radius, Coordinate3D position) {
    std::vector<std::tuple<int, int, int>> points = covered_points(radius, position);
    const auto&[n_x, n_y, n_z] = molecule_manager_->getGridSize();
    double mean_conc = 0;
    int count = 0;

    if(auto &pair = *molecule_manager_->find(name); pair != *molecule_manager_->end()){
        auto &data = pair.second;

        for (auto point: points){
            int i = std::get<2>(point);
            int j = std::get<1>(point);
            int k = std::get<0>(point);
            mean_conc += data[k+ n_x *j+ n_y * n_x *i];
            ++count;
        }
    }
    mean_conc /= count;
    return (mean_conc);
}



Coordinate3D CuboidSite::getCloseToCenterPosition() {
    double x, y, z;
    switch (dimensions) {
    case 2:x = random_generator_->generateDouble(upper_bound_.x*0.75, lower_bound_.x*0.75);
        y = random_generator_->generateDouble(upper_bound_.y*0.75, lower_bound_.y*0.75);
        z = 0.0;
        break;


    case 3:x = random_generator_->generateDouble(upper_bound_.x*0.75, lower_bound_.x*0.75);
        y = random_generator_->generateDouble(upper_bound_.y*0.75, lower_bound_.y*0.75);
        z = random_generator_->generateDouble(upper_bound_.z*0.75, lower_bound_.z*0.75);
        break;

    default:x = random_generator_->generateDouble(upper_bound_.x/2, lower_bound_.x/2);
        y = random_generator_->generateDouble(upper_bound_.y/2, lower_bound_.y/2);
        z = random_generator_->generateDouble(upper_bound_.z/2, lower_bound_.z/2);
        break;
    }
    return Coordinate3D{x, y, z};
}

Coordinate3D CuboidSite::getZEqualsZeroPosition(double radius) {

    double   x = random_generator_->generateDouble(upper_bound_.x + radius, lower_bound_.x - radius);
    double   y = random_generator_->generateDouble(upper_bound_.y + radius, lower_bound_.y - radius);

    return Coordinate3D{x, y, 0.0};
}

Coordinate3D CuboidSite::getZandYequalsZeroPosition(double radius) {

    double   x = random_generator_->generateDouble(upper_bound_.x + radius, lower_bound_.x - radius);

    return Coordinate3D{x, 0.0, 0.0};
}

void CuboidSite::initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                                  double current_time,
                                  double time_delta){
    for (const auto &agent_parameters: parameters.agents) {
        auto init_distribution = agent_parameters->initial_distribution;
        double pos_dist = agent_parameters->pos_dist;
        const auto name = agent_parameters->type;
        const auto number_of_agents = agent_parameters->number;
        for (auto i = 0; i < number_of_agents; ++i) {
            auto count_cell = 0;
            auto radius = agent_parameters->morphology_parameters.radius;
            Coordinate3D initial_position;
            if (i >= 1 && init_distribution == 1){
                SYSTEM_STDOUT("Initial distribution not possible: No more than 1 cell can be in the center of the environment.\n Switching to random position for remaining cells.");
                init_distribution = 0;
            }
            do {
                switch (init_distribution) {
                    case 0:
                        initial_position = getRandomPosition(radius);
                        break;
                    case 1:
                        initial_position = getCenterPosition();
                        break;
                    case 2:
                        initial_position = getCloseToCenterPosition();
                        break;
                    case 3:
                        initial_position = getZEqualsZeroPosition(radius);
                        break;
                    case 4:
                        initial_position = getZandYequalsZeroPosition(radius);
                        break;
                    case 5:
                        ++count_cell;
                        initial_position = two_cells_layout(radius, count_cell, pos_dist);
                        break;
                    default:
                        initial_position = getRandomPosition(radius);
                        break;
                }
            } while (check_collision(initial_position));
            auto agent = agent_manager_->emplace_back(cell_factory_->createCell(name,
                                                                                std::make_unique<Coordinate3D>(initial_position),
                                                                                agent_manager_->generateNewID(),
                                                                                this,
                                                                                time_delta,
                                                                                current_time));
        }
    }
}

bool CuboidSite::check_collision(Coordinate3D pos_to_check) {
    bool check = false;
    for (const auto &ag: agent_manager_->getAllAgents()){
        if (ag->getPosition().calculateEuclidianDistance(pos_to_check) <= ag->getAgentProperties()->getMorphology()->getBasicSphereOfThis()->getRadius()*2) {
            check = true;
            break;
        }
    }
    return check;
}

Coordinate3D CuboidSite::two_cells_layout(double radius, int cell, double pos_dist) {
    double x;
    double l = lower_bound_.x - upper_bound_.x;
    double m = std::min(pos_dist, l-4*radius);
    if (m< 0){
        m = (l-4*radius)/3;
    }
    auto scale = (l-4*radius-m)/2;
    if (cell == 1){
        x = upper_bound_.x + radius +scale;
    }
    else{
        x = lower_bound_.x - scale - radius;
    }
    return Coordinate3D{x, 0.0, 0.0};
}

