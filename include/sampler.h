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
    //use ThermalParames    
    double get_max_w(const ThermalParams& params);
    void sample_momentum(const ThermalParams& params, double pLRF[4]);


private:

    unsigned sampler_seed;  // Store the random number
    double poly_neg[5] = {-0.04357527, 0.03834501, -0.07525446, 0.8744597, -0.69273194};
    double poly_pos[5] = {-3.74011558e-08, 1.26120477e-06, 1.01637978e-03, 7.61645176e-05, -1.25334615e-04};

    int tries;
    int accepted;

    Table2D max_w_table;


};



#endif // SAMPLER_H 