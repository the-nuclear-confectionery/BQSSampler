#ifndef SAMPLER_H
#define SAMPLER_H

#include <random>
#include <iostream>
#include <chrono>
#include <random>
#include <sstream>
#include <cstdio> 

#include "surface.h"
#include "fluid_cell.h"
#include "particle.h"
#include "particle_system.h"
#include "lrf.h"
#include "tools.h"
#include "tables.h"
#include "settings.h"
#include "deltaf_corrections.h"

class Sampler {
public:
    // Constructor that takes a Surface object reference
    Sampler(const Settings& settings);


    // Method to sample 
    void sample(ParticleSystem& particle_system, const Surface& surface, const NumericalIntegrator& integrator);
    void conserved_charge_sampling(ParticleSystem& particle_system,  Surface& surface, const NumericalIntegrator& integrator);
    void sample_unconstrained(ParticleSystem& all_particles, Surface& surface, const NumericalIntegrator& integrator);

    bool check_weight_region(double mbar, double omega);
    bool check_weight_region_massive(double mbar, double omega);
 
    double get_max_w(const ThermalParams& params);
    double get_max_w_massive(const ThermalParams& params);

    void sample_momentum(const ThermalParams& params, double pLRF[4], std::mt19937& generator_momentum);


    void save_particles(const std::string& filename) const;
    void save_particles(const std::string& filename,
                        const std::vector<std::vector<Particle>>& particles) const;
    void calculate_integrated_ntot(ParticleSystem& particle_system, Surface& surface, const NumericalIntegrator& integrator);

    void check_total_charge_average(double totalB, double totalS, double totalQ, int Nsamples);



    std::vector<Particle> sample_fixed_yield_from_surface(
    const ParticleSystem& group,
    const Surface& surface,
    const std::vector<double>& N_cell_vector,
    const std::vector<std::vector<double>>& N_species_cell,
    int required,
    const std::string& coord,
    std::vector<Particle>* neg_out = nullptr);

    // Sample particles from a single fluid cell. The mean particle number is
    // computed internally (thermal integral + Poisson draw), so callers only
    // need to supply the cell, the particle group, and the integrator.
    std::vector<Particle> sample_cell(
    const FluidCell& cell,
    ParticleSystem& group,
    const NumericalIntegrator& integrator);

    double net_baryon(const std::vector<Particle>& particles);
    double net_strangeness(const std::vector<Particle>& particles);
    double net_charge(const std::vector<Particle>& particles);

private:
    // Core single-cell sampling routine shared by sample_cell and the
    // surface-loop entry points. Takes a precomputed mean yield and species
    // weights, draws a Poisson number of hadrons, and appends accepted
    // particles to `out`. `max_to_add` caps how many particles this call may
    // append (-1 = the full Poisson draw); used by fixed-yield callers to stop
    // exactly on target while preserving the RNG stream.
    void sample_cell_core(
    const FluidCell& cell,
    const ParticleSystem& group,
    double N_tot_mean,
    const std::vector<double>& species_weights,
    const std::string& coord,
    std::vector<Particle>& out,
    int max_to_add = -1,
    std::vector<Particle>* neg_out = nullptr);

    //settings reference
    const Settings& settings;

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

    double y_max;
    int D;
    int Nsamples;

    Table2D max_w_table;
    Table2D max_w_table_massive;

    Table4D deltaf_table;

    //sampled particle vectors
    std::vector< std::vector <Particle> > sampled_particles;
    // negative Cooper-Frye contributions (p.dSigma < 0), one vector per event
    std::vector< std::vector <Particle> > negative_particles;

    std::mt19937 gen_poisson;
    std::mt19937 gen_type;
    std::mt19937 gen_mom;
    std::mt19937 gen_keep;
    std::mt19937 gen_y;
    std::mt19937 gen_trim;
    std::mt19937 gen_pos;
    std::mt19937 gen_pos_neg;   // position smearing for negative particles (keeps gen_pos stream unperturbed)

    double pos_smearing;
    bool sample_negative;       // opt-in: also collect/write p.dSigma < 0 contributions

};



#endif // SAMPLER_H 