#include "sampler.h"
#include <chrono>


// Constructor implementation
Sampler::Sampler(){
    // Initialize the random number generator with a seed based on the current time

    sampler_seed =  static_cast<unsigned>(
            std::chrono::system_clock::now().time_since_epoch().count() );
            Tools::read_2d_table("../utils/max_w_table.dat", mbar_vals_wtable, omega_vals_wtable, max_w_values);
    
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
    int npions = 0;
    int nabove = 0;
    int nabove2 = 0;
    int nabovepp = 0;
    int nabovep0 = 0;
    int nabovepm = 0;
    int nabovekp = 0;
    int nabovek0 = 0;
    int nabovekm = 0;
    int naboveantik0 = 0;
    int nabovebosons = 0;
    int nabovefermions = 0;
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
        lrf.compute_dsigma_magnitude();
        //if(icell == 20424){
        //    std::cout << "dsigma_magnitude: " << lrf.dsigma_magnitude << std::endl;
        //    std::cout << "N_tot_cell: " << N_tot_cell << std::endl;
        //    std::cout << "N_tot_cell * 2.0 * y_max: " << N_tot_cell * 2.0 * y_max << std::endl;
        //}


        N_tot_cell *= 2.0 * y_max * lrf.dsigma_magnitude;
        Ncell += N_tot_cell;
        if(N_tot_cell <=0.) continue;

        std::poisson_distribution<int> poisson_hadrons(N_tot_cell);
        std::discrete_distribution<int> particle_type_distribution(
                                                particle_system.particle_species_number.begin(), 
                                                particle_system.particle_species_number.end());
        int N_hadrons = poisson_hadrons(generator_poisson);
        Ntot += N_hadrons;
        for (int i = 0; i < N_hadrons; i++) {
            // Sample particle type
            int sampled_index = particle_type_distribution(generator_type);
            int pid = particle_system.pid[sampled_index];
            double mass = particle_system.mass[sampled_index];
            double theta = particle_system.theta[sampled_index];
            double baryon = particle_system.baryon[sampled_index];
            double strange = particle_system.strange[sampled_index];
            double charge = particle_system.charge[sampled_index];
            if(pid == 211 || pid == -211 || pid == 111) npions++;
            //check weight region
            double mbar = mass / T;
            double omega = (baryon * muB + strange * muS + charge * muQ)/T;
            double max_w =1.0;
            if(theta == -1.0){
                if(check_weight_region(mbar, omega)){
                    nabove++;
                    // interpolate max_w
                    max_w = get_max_w(mbar, omega);
                    // check for positive, neutral, negative pions
                    //if(pid == 211){
                    //    nabovepp++;
                    //    //print mbar, omega and max_w
                    //    std::cout << "mbar: " << mbar << ", omega: " << omega << ", max_w: " << max_w << std::endl;
                    //}
                    //if(pid == 111) nabovep0++;
                    //if(pid == -211) nabovepm++;
                    ////check for kaons
                    //if(pid == 321) nabovekp++;
                    //if(pid == 311) nabovek0++;
                    //if(pid == -321) nabovekm++;
                    //if(pid == -311) naboveantik0++;
                    //std::cout << "max_w: " << max_w << std::endl;
                }
            }

            bool rejected = false;
            double pbar, Ebar, kbar, phi_over_2pi, costheta, weight;
            //pions 
            ///@todo: check which value is correct
            if(mbar <= 1.008){
                while(rejected)
                {
                
                    double r1= 1. - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);
                    double r2= 1. - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);
                    double r3= 1. - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);

                    double l1 = log(r1);
                    double l2 = log(r2);
                    double l3 = log(r3);

                    double l1_plus_l2 = l1 + l2;

                    pbar = - (l1 + l2 + l3);
                    Ebar = sqrt(pbar * pbar  +  mbar*mbar);

                    phi_over_2pi = l1_plus_l2 * l1_plus_l2 / (pbar * pbar);
                    costheta = (l1 - l2) / l1_plus_l2;

                    weight = 1.0 / (exp(Ebar-omega) + theta) / max_w / (r1 * r2 * r3);


                    // check pLRF acceptance
                    if(std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum) < weight) break;

                } // rejection loop
            }
            else{
                //heavy particles
                std::uniform_real_distribution<double> costheta_distribution(-1.0, nextafter(1.0, std::numeric_limits<double>::max()));

                // integrated weights
                std::vector<double> K_weight;
                K_weight.resize(3);

                K_weight[0] = mbar*mbar; // heavy (sampled_distribution = 0)
                K_weight[1] = 2.0 * mbar;   // moderate (sampled_distribution = 1)
                K_weight[2] = 2.0;          // light (sampled_distribution = 2)

                std::discrete_distribution<int> K_distribution(K_weight.begin(), K_weight.end());

                double kbar;  // kinetic energy / T
                while(rejected)
                {
                    int sampled_distribution = K_distribution(generator_momentum);
                    // select distribution to sample from based on integrated weights
                    if(sampled_distribution == 0)
                    {
                      // draw k from exp(-k/T).dk by sampling r1
                      // sample direction uniformly
                      kbar = - log(1.0 - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum));
                      phi_over_2pi = std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);
                      costheta = costheta_distribution(generator_momentum);
                    } // distribution 1 (very heavy)
                    else if(sampled_distribution == 1)
                    {
                      // draw (k,phi) from k.exp(-k/T).dk.dphi by sampling (r1,r2)
                      // sample costheta uniformly
                      double l1 = log(1.0 - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum));
                      double l2 = log(1.0 - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum));
                      kbar = - (l1 + l2);
                      phi_over_2pi = - l1 / kbar;
                      costheta = costheta_distribution(generator_momentum);
                    } // distribution 2 (moderately heavy)
                    else if(sampled_distribution == 2)
                    {
                      // draw (k,phi,costheta) from k^2.exp(-k/T).dk.dphi.dcostheta by sampling (r1,r2,r3)
                      double l1 = log(1.0 - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum));
                      double l2 = log(1.0 - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum));
                      double l3 = log(1.0 - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum));
                      double l1_plus_l2 = l1 + l2;
                      kbar = - (l1 + l2 + l3);
                      phi_over_2pi = l1_plus_l2 * l1_plus_l2 / (kbar * kbar);
                      costheta = (l1 - l2) / l1_plus_l2;
                    } // distribution 3 (light)
                    else
                    {
                      printf("Error: did not sample any K distribution\n");
                      exit(-1);
                    }
                    Ebar = kbar + mbar;                        // energy / T
                    pbar = sqrt(Ebar * Ebar  -  mbar*mbar); // momentum magnitude / T
                    double exponent = exp(Ebar - omega);
                    double weight = pbar/Ebar * exponent / (exponent + theta) / max_w;
                    //if(fabs(weight - 0.5) > 0.5) printf("Sample momemtum error: weight = %f out of bounds\n", weight);
                    // check pLRF acceptance
                    if(std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum) < weight) break;
                } // rejection loop
            }
            double E = Ebar * T;
            double p = pbar * T;
            double phi = phi_over_2pi * 2.0 * M_PI;
            double sintheta = sqrt(1.0  -  costheta * costheta);    // sin(theta)

            double px_lrf = p * sintheta * cos(phi);
            double py_lrf = p * sintheta * sin(phi);
            double pz_lrf = p * costheta;


        

  
        }



        if (icell % progress_step == 0) {
            std::cout << "Progress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "% completed" << std::endl;
        }
    }

    std::cout << "Total number of particles sampled: " << Ntot << std::endl;
    std::cout << "Total number of particles from cells: " << Ncell << std::endl;
    std::cout << "Total number of pions: " << npions << std::endl;
    std::cout << "Total number of particles above boundary: " << nabove << std::endl;
    std::cout << "Total number of particles above boundary (pions): " << nabovepp << std::endl;
    std::cout << "Total number of particles above boundary (neutral): " << nabovep0 << std::endl;
    std::cout << "Total number of particles above boundary (negative): " << nabovepm << std::endl;
    std::cout << "Total number of particles above boundary (positive kaons): " << nabovekp << std::endl;
    std::cout << "Total number of particles above boundary (neutral kaons): " << nabovek0 << std::endl;
    std::cout << "Total number of particles above boundary (anti kaons): " << naboveantik0 << std::endl;
    std::cout << "Total number of particles above boundary (negative kaons): " << nabovekm << std::endl;
    std::cout << "Total number of particles above boundary (bosons): " << nabovebosons << std::endl;
    std::cout << "Total number of particles above boundary (fermions): " << nabovefermions << std::endl;
    
    std::cout << "Total number of particles above boundary2: " << nabove2 << std::endl;
}


bool Sampler::check_weight_region(double mbar, double omega) {
    double boundary_value;
    if (omega < 0) {
        boundary_value = poly_neg[0] *pow(mbar, 4) + poly_neg[1] * pow(mbar, 3) + poly_neg[2] * pow(mbar, 2) + poly_neg[3] * mbar + poly_neg[4];
    } else {
        boundary_value = poly_pos[0] *pow(mbar, 4) + poly_pos[1] * pow(mbar, 3) + poly_pos[2] * pow(mbar, 2) + poly_pos[3] * mbar + poly_pos[4];
    }
    return omega>boundary_value;
}


double Sampler::get_max_w(double mbar, double omega) {
    return Tools::bilinear_interpolate(mbar_vals_wtable, omega_vals_wtable, max_w_values, mbar, omega);
}
