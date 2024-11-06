//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.


#include "MoleculeManager.h"

#include <cmath>
#include <sstream>

#include "basic/Coordinate3D.h"

MoleculeManager::MoleculeManager(abm::util::SimulationParameters::MoleculeManagerParameters molecules_manager_parameter, std::tuple<int,int,int> grid, Coordinate3D upper, Coordinate3D lower) : grid_(std::move(grid)){

    auto distance_x = (std::abs(upper.x) + std::abs(lower.x)) / std::get<0>(grid);
    auto distance_y = (std::abs(upper.y) + std::abs(lower.y)) / std::get<1>(grid);
    auto distance_z = (std::abs(upper.z) + std::abs(lower.z)) / std::get<2>(grid);
    if (std::get<0>(grid) == 0) {
        distance_x = 0;
    }
    if (std::get<1>(grid) == 0) {
        distance_y = 0;
    }
    if (std::get<2>(grid) == 0){
       distance_z = 0;
    }
    size_ = std::make_tuple(distance_x,distance_y,distance_z);
    align_factors = std::make_tuple(upper.x + distance_x/2, upper.y + distance_y/2, upper.z + distance_z/2);
    //align_factors = std::make_tuple(upper.x, upper.y, upper.z);
    Coordinate3D tmp;
    for(int i = 0; i < std::get<2>(grid); ++i){
        tmp.z=i * distance_z + std::get<2>(align_factors);
        for(int j = 0; j < std::get<1>(grid); ++j){
            tmp.y = j * distance_y + std::get<1>(align_factors);
            for(int k = 0; k < std::get<0>(grid); ++k){
                tmp.x = k* distance_x + std::get<0>(align_factors);
                location_.push_back(tmp);
            }
        }
    }

    //interactions
    for(auto &[name, data]:molecules_manager_parameter.molecule_molecule_interactions) {
        molecule_molecule_interaction.emplace(std::pair<std::string, abm::util::SimulationParameters::MoleculeMoleculeInteractionParameters>(name, data));
        concentrations_[name] =  std::vector<double>(std::get<0>(grid)*std::get<1>(grid)*std::get<2>(grid),0.0);
        concentrations_prev_[name] = std::vector<double>(std::get<0>(grid)*std::get<1>(grid)*std::get<2>(grid),0.0);
    }

    //Set Parameters
    for(const auto molecule: molecules_manager_parameter.molecules){
        std::string identifier = molecule->name;
        molecule_names.push_back(identifier);
        concentrations_[identifier] =  std::vector<double>(std::get<0>(grid)*std::get<1>(grid)*std::get<2>(grid),molecule->initial_concentration);
        concentrations_prev_[identifier] =  std::vector<double>(std::get<0>(grid)*std::get<1>(grid)*std::get<2>(grid),molecule->initial_concentration);
        decay[identifier] = molecule->decay;
        diffusion_coefficients[identifier] = molecule->diffusion_coefficient;
        initial_concentration[identifier] = molecule->initial_concentration;
        flow[identifier] = molecule->boundary_flow;
        spontaneous[identifier] = molecule->spontaneous_rate;

        setInitialDistribution(initial_concentration);
    }
    boundary_condition = molecules_manager_parameter.boundary_condition;
}

std::tuple<double,double> MoleculeManager::getParameters(const std::string &identifier){
    return std::make_tuple(diffusion_coefficients[identifier],decay[identifier]);

}

MoleculeManager::iterator MoleculeManager::find(const std::string & name){
    return concentrations_.find(name);
}
MoleculeManager::iterator MoleculeManager::end(const std::string & name){
    return concentrations_.end();
}

std::tuple<int, int, int> MoleculeManager::getGridSize() {return grid_;}
std::tuple<double, double, double> MoleculeManager::getDistances() {return size_;}

void MoleculeManager::setInitialDistribution(std::map<std::string , double> value) {
    const auto&[n_x,n_y,n_z] = grid_;
    auto volume = (std::get<0>(size_)*n_x) * (std::get<1>(size_)*n_y) * (std::get<2>(size_)*n_z);
    for(auto &[name,vec]: value){
        for (int i = 0; i < n_x * n_y * n_z; ++i){
            auto val = value[name];//*volume/(n_x * n_y * n_z);
            concentrations_.at(name)[i] = val;
            //std::cout << "Init conc per grid point = " << val << std::endl;
        }
    }
}

void MoleculeManager::initial_remove_mol_inside_cell(std::vector<std::tuple<int, int, int>> covered_points) {
    const auto&[n_x,n_y,n_z] = grid_;
    for (auto name: molecule_names) {
        for (auto points: covered_points){
            int i = std::get<0>(points);
            int j = std::get<1>(points);
            int k = std::get<2>(points);
            concentrations_.at(name)[k + n_x * j + n_y * n_x * i] = 0.0;
        }
    }
}

void MoleculeManager::do_reactions(double time_delta, double radius){
    const auto&[n_x,n_y,n_z] = grid_;
    for (const auto &[name, interaction] : molecule_molecule_interaction){
        auto binding_rate = interaction.binding;
        auto unbinding_rate = interaction.unbinding;
        auto first_mol = interaction.first;
        auto second_mol = interaction.second;

        for (int i = 0; i < n_z; ++i) {
            for (int j = 0; j < n_y; ++j) {
                for (int k = 0; k < n_x; ++k) {

                    // COMPLEX FORMATION
                    double complex_formation = time_delta * binding_rate * concentrations_.at(first_mol)[k + n_x * j + n_y * n_x * i] * concentrations_.at(second_mol)[k + n_x * j + n_y * n_x * i];
                    double complex_unbinding = time_delta * unbinding_rate * concentrations_.at(name)[k + n_x * j + n_y * n_x * i];

                    concentrations_.at(name)[k + n_x * j + n_y * n_x * i] += complex_formation - complex_unbinding;
                    concentrations_.at(first_mol)[k + n_x * j + n_y * n_x * i] += complex_unbinding - complex_formation;
                    concentrations_.at(second_mol)[k + n_x * j + n_y * n_x * i] += complex_unbinding - complex_formation;
                }
            }
        }
    }
}

void MoleculeManager::doDiffusion(double time_delta, double current_time, std::vector<std::tuple<int, int, int>> covered_points, std::vector<std::tuple<int, int, int>> surface_points, double radius) {
    const auto& [n_x, n_y, n_z] = getGridSize();
    std::vector<double> tmp(n_x * n_y * n_z, 0);

    const auto& [h_x, h_y, h_z] = getDistances();
    auto i_min = 0, j_min = 0, k_min = 0, i_max = 0, j_max = 0, k_max = 0;
    // Cube around the cells
    if (!covered_points.empty() && !surface_points.empty()){
        i_min = std::get<0>(*std::min_element(std::begin(surface_points), std::end(surface_points), [](auto lhs, auto rhs) {
            return (std::get<0>(lhs) < std::get<0>(rhs));}));
        j_min = std::get<1>(*std::min_element(std::begin(surface_points), std::end(surface_points), [](auto lhs, auto rhs) {
            return (std::get<1>(lhs) < std::get<1>(rhs));}));
        k_min = std::get<2>(*std::min_element(std::begin(surface_points), std::end(surface_points), [](auto lhs, auto rhs) {
            return (std::get<2>(lhs) < std::get<2>(rhs));}));
        i_max = std::get<0>(*std::min_element(std::begin(surface_points), std::end(surface_points), [](auto lhs, auto rhs) {
            return (std::get<0>(lhs) > std::get<0>(rhs));}));
        j_max = std::get<1>(*std::min_element(std::begin(surface_points), std::end(surface_points), [](auto lhs, auto rhs) {
            return (std::get<1>(lhs) > std::get<1>(rhs));}));
        k_max = std::get<2>(*std::min_element(std::begin(surface_points), std::end(surface_points), [](auto lhs, auto rhs) {
            return (std::get<2>(lhs) > std::get<2>(rhs));}));
    }
    tmp.reserve(n_x * n_y * n_z);

    for (auto& [name, data] : *this) {

        const auto& [diffusion_coefficient, decay] = getParameters(name);

        const double scalar = diffusion_coefficient * time_delta;

        const double spontaneous_val = spontaneous[name] *time_delta;

        //ALL GRID POINTS AROUND CUBE OF INTERESTS
        for (int i = 1; i < i_min; ++i) {
            for (int j = 1; j < n_y-1; ++j) {
                for (int k = 1; k < n_x-1; ++k) {
                    // Diffusion 3D
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Spontaneous creation
                    tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                }
            }
        }
        for (int i = i_min; i < n_z-1; ++i) {
            for (int j = 1; j < j_min; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    // Diffusion 3D
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Spontaneous creation
                    tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                }
            }
        }

        for (int i = i_min; i < n_z - 1; ++i) {
            for (int j = j_max + 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    // Diffusion 3D
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Spontaneous creation
                    tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                }
            }
        }
        for (int i = i_max + 1; i < n_z - 1; ++i) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    // Diffusion 3D
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Spontaneous creation
                    tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                }
            }
        }
        for (int i = i_min; i <= i_max; ++i) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int k = 1; k < k_min; ++k) {
                    // Diffusion 3D
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Spontaneous creation
                    tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                }
            }
        }
        for (int i = i_min; i <= i_max; ++i) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int k = k_max + 1; k < n_x-1 ; ++k) {
                    // Diffusion 3D
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Spontaneous creation
                    tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                }
            }
        }
        // CUBE around CELL
        for (int i = i_min; i <= i_max; ++i) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int k = k_min; k <= k_max; ++k) {
                    if (std::find(surface_points.begin(), surface_points.end(), std::tuple<int, int, int>(i, j, k)) != surface_points.end()){
                        // CELL BOUNDARY CONDITIONS: REFLECTING
                        double i_count = 2.0, j_count = 2.0, k_count = 2.0;

                        tmp[k + n_x * j + n_y * n_x * i] = 0.0;

                        if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i - 1, j, k)) != covered_points.end()) {
                            --i_count;
                        }
                        else {
                            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                        }
                        if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i + 1, j, k)) != covered_points.end()) {
                            --i_count;
                        }
                        else {
                            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                        }
                        if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i, j - 1, k)) != covered_points.end()) {
                            --j_count;
                        }
                        else {
                            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j-1) + n_y * n_x * i] / (h_z * h_z);
                        }
                        if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i, j + 1, k)) != covered_points.end()) {
                            --j_count;
                        }
                        else {
                            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j+1) + n_y * n_x * i] / (h_z * h_z);
                        }
                        if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i, j, k - 1)) != covered_points.end()) {
                            --k_count;
                        }
                        else {
                            tmp[k + n_x * j + n_y * n_x * i] += data[k-1 + n_x * j + n_y * n_x * i] / (h_z * h_z);
                        }
                        if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i, j, k + 1)) != covered_points.end()) {
                            --k_count;
                        }
                        else {
                            tmp[k + n_x * j + n_y * n_x * i] += data[k+1 + n_x * j + n_y * n_x * i] / (h_z * h_z);
                        }

                        tmp[k + n_x * j + n_y * n_x * i] += -k_count * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - j_count * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - i_count * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);
                        tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                        tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point

                        // Decay
                        tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    }
                    else if (std::find(covered_points.begin(), covered_points.end(), std::tuple<int, int, int>(i, j, k)) != covered_points.end()){
                        continue ;
                    }
                    else {
                        // Diffusion 3D
                        tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_x * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_x * j + n_y * n_x * i] / (h_z * h_z);

                        tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                        tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);

                        tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                        tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);

                        tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                        tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);

                        tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                        tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i]; // data from previous time at the same point
                        // Decay
                        tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                        // Spontaneous creation
                        tmp[k + n_x * j + n_y * n_x * i] += spontaneous_val;
                    }
                }
            }
        }

        // EXTERNAL BOUNDARY CONDITION
        auto volume_cell = (4.0/3.0) * M_PI * pow(radius, 3);
        auto volume = (std::get<0>(size_)*n_x) * (std::get<1>(size_)*n_y) * (std::get<2>(size_)*n_z) - volume_cell;
        auto surface = 2 * (std::get<0>(size_)*n_x * n_y) + 2* (std::get<1>(size_)*n_z*n_y) +2 * (std::get<0>(size_)*n_x * n_z);
        double flow_val = flow[name].amplitude; //* time_delta ;//* surface * (cos(flow[name].frequence * current_time)+1) /(2*n_x*n_y + 2*n_y*n_z + 2*n_z*n_x);
        if (flow_val != 0){
            if (flow[name].target_conc != 0){
                auto time_limit_flow = (flow[name].target_conc * volume)/(flow_val*surface); //To double-check
                if (current_time > time_limit_flow){
                    flow_val = 0;
                    flow[name].amplitude = 0;
                }
            }
        }
        //std::cout << "flow = " << flow_val << std::endl;
        if (boundary_condition == "periodic") {
            // 6 external faces (without edges)
            int i = 0;
            for (int j = 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }
            i = n_z - 1;
            for (int j = 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * 0] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            int j = 0;
            for (i = 1; i < n_z - 1; ++i) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }
            j = n_y - 1;
            for (i = 1; i < n_z - 1; ++i) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            int k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                for (j = 1; j < n_y - 1; ++j) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                for (j = 1; j < n_y - 1; ++j) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            // 12 edges
            i = 0;
            j = 0;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = n_z - 1;
            j = 0;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = 0;
            j = n_y - 1;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = n_z - 1;
            j = n_y - 1;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = 0;
            k = 0;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            i = 0;
            k = n_x - 1;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            i = n_z - 1;
            k = 0;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            i = n_z - 1;
            k = n_x - 1;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = 0;
            k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = 0;
            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = n_y - 1;
            k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = n_y - 1;
            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            // 8 corner points
            i = 0;
            j = 0;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = 0;
            j = n_y - 1;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = 0;
            j = n_y - 1;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = 0;
            j = 0;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (n_z - 1)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = 0;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = 0;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (n_y - 1) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = n_y - 1;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[n_x - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = n_y - 1;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[0 + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (0) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (0)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
        }
        else if (boundary_condition == "reflecting") {
            //6 external faces (without edges)
            int i = 0;
            for (int j = 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }
            i = n_z - 1;
            for (int j = 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            int j = 0;
            for (i = 1; i < n_z - 1; ++i) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }
            j = n_y - 1;
            for (i = 1; i < n_z - 1; ++i) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            int k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                for (j = 1; j < n_y - 1; ++j) {
                    tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                for (j = 1; j < n_y - 1; ++j) {
                    tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            // 12 edges
            i = 0;
            j = 0;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = n_z - 1;
            j = 0;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = 0;
            j = n_y - 1;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = n_z - 1;
            j = n_y - 1;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }
            i = 0;
            k = 0;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            i = 0;
            k = n_x - 1;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            i = n_z - 1;
            k = 0;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            i = n_z - 1;
            k = n_x - 1;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = 0;
            k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = 0;
            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = n_y - 1;
            k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            j = n_y - 1;
            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
            }

            // 8 corner points
            i = 0;
            j = 0;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = 0;
            j = n_y - 1;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = 0;
            j = n_y - 1;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = 0;
            j = 0;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = 0;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = 0;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = n_y - 1;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }

            i = n_z - 1;
            j = n_y - 1;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -1.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 1.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];

            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
        }
        else { // ABSORBING
            //6 external faces (without edges)
            int i = 0;
            for (int j = 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }
            i = n_z - 1;
            for (int j = 1; j < n_y - 1; ++j) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            int j = 0;
            for (i = 1; i < n_z - 1; ++i) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }
            j = n_y - 1;
            for (i = 1; i < n_z - 1; ++i) {
                for (int k = 1; k < n_x - 1; ++k) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            int k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                for (j = 1; j < n_y - 1; ++j) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                for (j = 1; j < n_y - 1; ++j) {
                    tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k - 1 + n_x * j + n_y * n_x * i] / (h_x * h_x);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (j - 1) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i - 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (i + 1)] / (h_z * h_z);
                    tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                    tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                    // Decay
                    tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                    // Flow
                    if (name == "AMP" && flow_val != 0){
                        tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                    }
                }
            }

            // 12 edges
            i = 0;
            j = 0;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }
            i = n_z - 1;
            j = 0;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }
            i = 0;
            j = n_y - 1;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }
            i = n_z - 1;
            j = n_y - 1;
            for (k = 1; k < n_x - 1; ++k) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }
            i = 0;
            k = 0;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            i = 0;
            k = n_x - 1;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            i = n_z - 1;
            k = 0;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            i = n_z - 1;
            k = n_x - 1;
            for (j = 1; j < n_y - 1; ++j) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            j = 0;
            k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            j = 0;
            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            j = n_y - 1;
            k = 0;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            j = n_y - 1;
            k = n_x - 1;
            for (i = 1; i < n_z - 1; ++i) {
                tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
                tmp[k + n_x * j + n_y * n_x * i] *= scalar;
                tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
                // Decay
                tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
                // Flow
                if (name == "AMP" && flow_val != 0){
                    tmp[k + n_x * j + n_y * n_x * i] = flow_val;
                }
            }

            // 8 corner points
            i = 0;
            j = 0;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = 0;
            j = n_y - 1;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = 0;
            j = n_y - 1;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = 0;
            j = 0;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = n_z - 1;
            j = 0;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = n_z - 1;
            j = 0;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = n_z - 1;
            j = n_y - 1;
            k = 0;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }

            i = n_z - 1;
            j = n_y - 1;
            k = n_x - 1;
            tmp[k + n_x * j + n_y * n_x * i] = -2.0 * data[k + n_y * j + n_y * n_x * i] / (h_x * h_x) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_y * h_y) - 2.0 * data[k + n_y * j + n_y * n_x * i] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] += data[-1 + k + n_x * j + n_y * n_x * i] / (h_x * h_x);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * (-1 + j) + n_y * n_x * i] / (h_y * h_y);
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * (-1 + i)] / (h_z * h_z);
            tmp[k + n_x * j + n_y * n_x * i] *= scalar;
            tmp[k + n_x * j + n_y * n_x * i] += data[k + n_x * j + n_y * n_x * i];
            // Decay
            tmp[k + n_x * j + n_y * n_x * i] -= decay * data[k + n_x * j + n_y * n_x * i] * time_delta;
            // Flow
            if (name == "AMP" && flow_val != 0){
                tmp[k + n_x * j + n_y * n_x * i] = flow_val;
            }
        }
        tmp.swap(data);
    }
}

double MoleculeManager::uptakeMolecules(std::string name,  double rate, std::vector<std::tuple<int, int, int>> surface_points, double timestep){
    const auto&[n_x, n_y, n_z] = getGridSize();

    double uptake = 0.0;
    double total_uptake = 0.0;
    for (auto& [name2, data] : *this) {
        if (name == name2) {
            for (auto points: surface_points){
                int i = std::get<0>(points);
                int j = std::get<1>(points);
                int k = std::get<2>(points);
                uptake = timestep * rate * data[k + n_x * j + n_y * n_x * i];
                total_uptake += uptake;
                data[k + n_x * j + n_y * n_x * i] -= uptake;
            }
        }
    }
    return total_uptake;
}

void MoleculeManager::secreteMolecules(const std::string& name,  double rate, std::vector<std::tuple<int, int, int>> surface_points, double timestep, double current_time, double radius){
    const auto&[n_x, n_y, n_z] = getGridSize();
    double secretion = timestep * rate/surface_points.size();
    for (auto& [name2, data] : *this) {
        if (name == name2) {
            for (auto points: surface_points){
                int i = std::get<0>(points);
                int j = std::get<1>(points);
                int k = std::get<2>(points);
                data[k + n_x * j + n_y * n_x * i] += secretion;
            }
        }
    }
}

Coordinate3D MoleculeManager::get_loc(int i, int j, int k) const {
    float z = i * std::get<2>(size_) + std::get<2>(align_factors);
    float y = j * std::get<1>(size_) + std::get<1>(align_factors);
    float x = k * std::get<0>(size_) + std::get<0>(align_factors);
    return Coordinate3D{x, y, z};
}

bool MoleculeManager::steady_state_reached(abm::util::SimulationParameters::StoppingCriteria& stopping_criteria, double delta_t){
    if (stopping_criteria.molecule == ""){
        return false;
    }
    else{
        auto name = stopping_criteria.molecule;
        double sum = 0.0;
        for (const auto &value: concentrations_.at(name)) {
            sum += value;
        }
        if (sum == 0.0) {
            return false;
        }
        else {
            double conc_now = sum / concentrations_.at(name).size();
            double sum2 = 0.0;
            for (const auto &val: concentrations_prev_.at(name)) {
                sum2 += val;
            }
            double conc_prev = sum2 / concentrations_prev_.at(name).size();
            auto diff = std::abs(conc_now - conc_prev);
            if (diff < stopping_criteria.threshold){
                stopping_criteria.time -= delta_t;
            }
            else{
                set_previous_conc();
                stopping_criteria.time = stopping_criteria.init_time;
                //std::cout << "TIME RESET" << std::endl;
            }
        }
    }
    //std::cout << "TIME BEFORE STEADY STATE REACHED" << stopping_criteria.time << std::endl;
    return (stopping_criteria.time <= 0);
}

void MoleculeManager::set_previous_conc() {
    concentrations_prev_ = concentrations_;
}
