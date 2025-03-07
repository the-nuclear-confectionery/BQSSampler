#ifndef SAMPLER_H
#define SAMPLER_H

#include <random>
#include <iostream>
#include <chrono>
#include <random>


#include "surface.h"
#include "particle.h"
#include "particle_system.h"
#include "lrf.h"
#include "tools.h"
#include "tables.h"


class Sampler {
public:
    // Constructor that takes a Surface object reference
    Sampler();


    // Method to sample 
    void sample(ParticleSystem& particle_system, const Surface& surface, const NumericalIntegrator& integrator);

    bool check_weight_region(double mbar, double omega);
    bool check_weight_region_massive(double mbar, double omega);
 
    double get_max_w(const ThermalParams& params);
    double get_max_w_massive(const ThermalParams& params);

    void sample_momentum(const ThermalParams& params, double pLRF[4], std::default_random_engine& generator_momentum);


private:

    unsigned sampler_seed;  // Store the random number
    double poly_neg[5] = {-0.04357527, 0.03834501, -0.07525446, 0.8744597, -0.69273194};
    double poly_pos[5] = {-3.74011558e-08, 1.26120477e-06, 1.01637978e-03, 7.61645176e-05, -1.25334615e-04};

    double poly_pos_massive_2[2] = {1.0227939, -0.23578976};
    double poly_pos_massive_1[4] = {0.02973458, -0.34499695 , 2.36223265,-2.00694215};
    double poly_neg_massive[4] ={9.19848542, -20.06276524, 17.01588659, -6.03130873};

    int tries;
    int accepted;

    int nabove;
    int nabove_massive;

    Table2D max_w_table;
    Table2D max_w_table_massive;

    std::vector<int> sampled_pid;
    std::vector<double> sampled_px;
    std::vector<double> sampled_py;
    std::vector<double> sampled_pz;
    std::vector<double> sampled_E;
    std::vector<double> sampled_x;
    std::vector<double> sampled_y;
    std::vector<double> sampled_z;
    std::vector<double> sampled_t;
    std::vector<double> sampled_mass;


};



#endif // SAMPLER_H 