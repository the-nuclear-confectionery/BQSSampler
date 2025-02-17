#include "sampler.h"


// Constructor implementation
Sampler::Sampler(){
    // Initialize the random number generator with a seed based on the current time

    sampler_seed =  static_cast<unsigned>(
            std::chrono::system_clock::now().time_since_epoch().count() );
}

void Sampler::sample(ParticleSystem& particle_system, const Surface& surface, const NumericalIntegrator& integrator) {
    double y_max = 5.0; // Maximum rapidity
    std::default_random_engine generator_poisson(sampler_seed);
    std::default_random_engine generator_type(sampler_seed + 10000);
    std::default_random_engine generator_momentum(sampler_seed + 20000);
    std::default_random_engine generator_keep(sampler_seed + 30000);
    std::default_random_engine generator_rapidity(sampler_seed + 40000);

    std::cout << "Sampling particles..." << std::endl;
    int progress_step = surface.npoints / 20; // Step size for 5% progress updates

    int Ntot = 0;
    double Ncell = 0.;
    for (int icell = 0; icell < surface.npoints; icell++) {
        LRF lrf(surface.ut[icell],
                surface.ux[icell],
                surface.uy[icell],
                surface.ueta[icell],
                surface.tau[icell]);

        if (icell % (surface.npoints / 20) == 0) {
            std::cout << "Progress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "%" << std::endl;
        }

        double tau = surface.tau[icell];
        double tau_squared = tau * tau;
        double u_dot_dsigma = surface.u_dot_dsigma[icell];
        double T = surface.T[icell];
        double muB = surface.muB[icell];
        double muS = surface.muS[icell];
        double muQ = surface.muQ[icell];

        double N_tot_cell = particle_system.calculate_particle_number(T, muB, muS, muQ, integrator);

        lrf.dsigma_t = surface.dsigma_t[icell];
        lrf.dsigma_x = surface.dsigma_x[icell];
        lrf.dsigma_y = surface.dsigma_y[icell];
        lrf.dsigma_n = surface.dsigma_eta[icell];

        lrf.boost_dsigma_to_lrf(tau_squared);
        if(icell == 20424){
            std::cout << "surface dsigma_t: " << surface.dsigma_t[icell] << std::endl;
            std::cout << "surface dsigma_x: " << surface.dsigma_x[icell] << std::endl;
            std::cout << "surface dsigma_y: " << surface.dsigma_y[icell] << std::endl;
            std::cout << "surface dsigma_n: " << surface.dsigma_eta[icell] << std::endl;
            std::cout << "dsigma_t: " << lrf.dsigma_t << std::endl;
            std::cout << "dsigma_x: " << lrf.dsigma_x << std::endl;
            std::cout << "dsigma_y: " << lrf.dsigma_y << std::endl;
            std::cout << "dsigma_n: " << lrf.dsigma_n << std::endl;
        }
        lrf.compute_dsigma_magnitude();
        if(icell == 20424){
            std::cout << "dsigma_magnitude: " << lrf.dsigma_magnitude << std::endl;
            std::cout << "N_tot_cell: " << N_tot_cell << std::endl;
            std::cout << "N_tot_cell * 2.0 * y_max: " << N_tot_cell * 2.0 * y_max << std::endl;
        }


        N_tot_cell *= 2.0 * y_max * lrf.dsigma_magnitude;
        Ncell += N_tot_cell;
        if(N_tot_cell <=0.) continue;

        std::poisson_distribution<int> poisson_hadrons(N_tot_cell);
        int N_hadrons = poisson_hadrons(generator_poisson);
        Ntot += N_hadrons;
        if (icell % progress_step == 0) {
            std::cout << "Progress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "% completed" << std::endl;
        }
    }

    std::cout << "Total number of particles sampled: " << Ntot << std::endl;
    std::cout << "Total number of particles from cells: " << Ncell << std::endl;
}