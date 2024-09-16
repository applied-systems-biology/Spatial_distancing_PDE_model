//   Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
//   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//   https://www.leibniz-hki.de/en/applied-systems-biology.html
//   HKI-Center for Systems Biology of Infection
//   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//   This code is licensed under BSD 2-Clause
//   See the LICENSE file provided with this code for the full license.

#ifndef AGENT_BASED_FRAMEWORK2_MOLECULEMANAGER_H
#define AGENT_BASED_FRAMEWORK2_MOLECULEMANAGER_H

#include <tuple>
#include <utility>
#include <vector>
#include <map>
#include <memory>
#include "utils/io_util.h"

class Site;
class Coordinate3D;


class MoleculeManager {
  public:
    // Class to handle all actions and parameters related to the molecules simulation.
    MoleculeManager(abm::util::SimulationParameters::MoleculeManagerParameters molecules_manager_parameter, std::tuple<int,int,int> grid, Coordinate3D upper, Coordinate3D lower);

    // map iterators for range loops in site (or in general)
    using iterator = std::map<std::string,std::vector<double>>::iterator;
    using const_iterator = std::map<std::string,std::vector<double>>::const_iterator;
    iterator begin() { return concentrations_.begin(); }
    [[nodiscard]] const_iterator begin() const { return concentrations_.begin(); }
    iterator end() { return concentrations_.end(); }
    [[nodiscard]] const_iterator end() const { return concentrations_.end(); }

    iterator find(const std::string & name);
    iterator end(const std::string & name);

    /*!
        * Set the initial condition of the molecule
        * @param value map with name of the molecule and initial concentration
        */
    void setInitialDistribution(std::map<std::string , double> value);

    void initial_remove_mol_inside_cell(std::vector<std::tuple<int, int, int>> covered_points);

    std::tuple<int,int,int> getGridSize();
    std::tuple<double,double,double> getDistances();
    std::tuple<double,double> getParameters(const std::string &key);

    [[nodiscard]] std::vector<Coordinate3D> get_location() const { return location_; }

    /*!
        * Molecule - molecule interactions: do the reaction dynamics between molecules (binding, unbinding, ...)
        * @param value map with name of the molecule and initial concentration
        */
    void do_reactions(double time_delta, double radius);

    [[nodiscard]] std::vector<std::string> get_molecules_names() const {return molecule_names;};

    std::string get_boundary_condition() const {return boundary_condition;};

    /*!
        * Molecule diffusion
        * @param time_delta timestep of the simulation
        * @param current_time current time of the simulation
        * @param covered_points Points that are covered by agents (where molecules should not diffuse)
        * @param surface_points Points that are on the surface of agents (ie. boundary)
        */
    void doDiffusion(double time_delta, double current_time, std::vector<std::tuple<int, int, int>> covered_points, std::vector<std::tuple<int, int, int>> surface_points, double radius);

    /*!
        * Molecule uptake by agent
        * @param name name of the molecule that is handled
        * @param rate uptake rate
        * @param surface_points Points that are on the surface of agents, where uptake is happening
        * @param timestep timestep of the simulation
        */
    double uptakeMolecules(std::string name,  double rate, std::vector<std::tuple<int, int, int>> surface_points, double timestep);

    /*!
        * Molecule secretion by agent
        * @param name name of the molecule that is handled
        * @param rate secretion rate
        * @param surface_points Points that are on the surface of agents, where secretion is happening
        * @param timestep timestep of the simulation
        * @param current_time current time of the simulation
        * @param radius radius of the cell (agent) that is secreting
        */
    void secreteMolecules(const std::string& name,  double rate, std::vector<std::tuple<int, int, int>> surface_points, double timestep, double current_time, double radius);

    [[nodiscard]] Coordinate3D get_loc(int i, int j, int k) const;

    [[nodiscard]] std::map<std::string,std::vector<double>> get_conc() const {return concentrations_;};

    bool steady_state_reached(abm::util::SimulationParameters::StoppingCriteria& stopping_criteria, double delta_t);
    void set_previous_conc();

private:

    std::map<std::string,std::vector<double>> concentrations_;
    std::map<std::string,std::vector<double>> concentrations_prev_;
    std::map<std::string,double> diffusion_coefficients;
    std::map<std::string,double> decay;
    std::map<std::string,double> initial_concentration;
    std::map<std::string, abm::util::SimulationParameters::Flow> flow;
    std::map<std::string, double> spontaneous;

    std::vector<Coordinate3D> location_;

    std::tuple<double,double,double> size_;

    std::tuple<double,double,double> align_factors;

    std::tuple<int,int,int> grid_;

    std::vector<std::string> molecule_names;

    std::map<std::string, abm::util::SimulationParameters::MoleculeMoleculeInteractionParameters> molecule_molecule_interaction;

    std::string boundary_condition;

    double dist;

};

#endif //AGENT_BASED_FRAMEWORK2_MOLECULEMANAGER_H
