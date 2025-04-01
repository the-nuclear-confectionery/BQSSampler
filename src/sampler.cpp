#include "sampler.h"



// Constructor implementation
Sampler::Sampler(){
    // Initialize the random number generator with a seed based on the current time

    sampler_seed =  static_cast<unsigned>(
            std::chrono::system_clock::now().time_since_epoch().count() );
            max_w_table.load_from_file("../utils/max_w_table.dat");
            max_w_table_massive.load_from_file("../utils/max_w_table_massive.dat");
    sampler_seed = 123123;
    accepted = 0;
    tries = 0;
    nabove = 0;
    nabove_massive = 0;
    
}

void Sampler::sample(ParticleSystem& particle_system, const Surface& surface, const NumericalIntegrator& integrator) {
    double y_max = 5.0; // Maximum rapidity
    int D = 3; // Number of dimensions
    int Nsamples = 200;
    y_max = 0.5;

    //resize the sampled particles vector
    sampled_particles.resize(Nsamples);

    std::default_random_engine generator_poisson(sampler_seed);
    std::default_random_engine generator_type(sampler_seed + 10000);
    std::default_random_engine generator_momentum(sampler_seed + 20000);
    std::default_random_engine generator_keep(sampler_seed + 30000);
    std::default_random_engine generator_rapidity(sampler_seed + 40000);

    std::cout << "Sampling particles..." << std::endl;
    int progress_step = surface.npoints / 20; // Step size for 5% progress updates

    int Ntot = 0;
    double Ncell = 0.;
    
    //variables for time acumulation
    double t_integrator = 0.;
    double t_momentum = 0.;
    double t_species = 0.;
    double t_total = 0.;

    auto start_total = std::chrono::high_resolution_clock::now();

    for (int icell = 0; icell < surface.npoints; icell++) {
        LRF lrf(surface.ut[icell],
                surface.ux[icell],
                surface.uy[icell],
                surface.ueta[icell],
                surface.dsigma_t[icell],
                surface.dsigma_x[icell],
                surface.dsigma_y[icell],
                surface.dsigma_eta[icell],
                surface.tau[icell]);
        double sinheta = sinh(surface.eta[icell]);
        double cosheta = sqrt(1.0  +  sinheta * sinheta);

        if (icell % (surface.npoints / 20) == 0) {
            std::cout << "Progress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "%" << std::endl;
        }

        //calculate udsigma 
        double udsigma = surface.ut[icell] * surface.dsigma_t[icell] 
                       + surface.ux[icell] * surface.dsigma_x[icell] 
                       + surface.uy[icell] * surface.dsigma_y[icell] 
                       + surface.ueta[icell] * surface.dsigma_eta[icell];
        if (udsigma <= 0.0) continue;
        double tau = surface.tau[icell];
        double tau_squared = tau * tau;
        double T = surface.T[icell];
        double muB = surface.muB[icell];
        double muS = surface.muS[icell];
        double muQ = surface.muQ[icell];

        auto start_integrator = std::chrono::high_resolution_clock::now();
        double N_tot_cell = particle_system.calculate_particle_number(T, muB, muS, muQ, integrator);
        auto end_integrator = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_integrator = end_integrator - start_integrator;
        t_integrator += elapsed_integrator.count();

        lrf.boost_dsigma_to_lrf(tau_squared);
        lrf.compute_dsigma_magnitude();
        //std::cout << "dsigma magnitude: " << lrf.dsigma_magnitude << std::endl;
        N_tot_cell *= 2.0 * y_max * lrf.dsigma_magnitude;
        //std::cout << "N_tot_cell: " << N_tot_cell << std::endl;
        if(N_tot_cell <=0.) continue;
        Ncell += N_tot_cell;
        //print all informations from surface
        //std::cout << "icell: " << icell << " tau: " << tau << " T: " << T << " muB: " << muB << " muS: " << muS << " muQ: " << muQ << std::endl;
        //std::cout << "dsigma_t: " << surface.dsigma_t[icell] << " dsigma_x: " << surface.dsigma_x[icell] << " dsigma_y: " << surface.dsigma_y[icell] << " dsigma_eta: " << surface.dsigma_eta[icell] << std::endl;
        //std::cout << "ut: " << surface.ut[icell] << " ux: " << surface.ux[icell] << " uy: " << surface.uy[icell] << " ueta: " << surface.ueta[icell] << std::endl;
        //std::cout << "dsigma_t_lrf: " << lrf.dsigma_t_lrf << " dsigma_x_lrf: " << lrf.dsigma_x_lrf << " dsigma_y_lrf: " << lrf.dsigma_y_lrf << " dsigma_z_lrf: " << lrf.dsigma_z_lrf << std::endl;
        //std::cout << "dsigma_magnitude: " << lrf.dsigma_magnitude << std::endl;
        //std::cout << "N_tot_cell: " << N_tot_cell << std::endl;

        std::poisson_distribution<int> poisson_hadrons(N_tot_cell);
        std::discrete_distribution<int> particle_type_distribution(
                                                particle_system.particle_species_number.begin(), 
                                                particle_system.particle_species_number.end());

        
        ThermalParams sampled_params;
        sampled_params.T = T;
        sampled_params.alphaB = muB / T;
        sampled_params.alphaQ = muQ / T;
        sampled_params.alphaS = muS / T;
        auto momentum_sampling_start = std::chrono::high_resolution_clock::now();
        for (int isample = 0; isample < Nsamples; isample++) {

            std::vector<Particle> single_sample_particles;

            int N_hadrons = poisson_hadrons(generator_poisson);
            Ntot += N_hadrons;

            //Create thermal parms struct for the sampled particles


            for (int i = 0; i < N_hadrons; i++) {
                // Sample particle type
                int sampled_index = particle_type_distribution(generator_type);

                int pid = particle_system.pid[sampled_index];

                double mass = particle_system.mass[sampled_index];
                //name
                std::string name = particle_system.name[sampled_index];
                //double theta = particle_system.theta[sampled_index];
                //double baryon = particle_system.baryon[sampled_index];
                //double strange = particle_system.strange[sampled_index];
                //double charge = particle_system.charge[sampled_index];

                sampled_params.mbar     = particle_system.mass[sampled_index] / T;
                sampled_params.baryon   = particle_system.baryon[sampled_index];
                sampled_params.strange  = particle_system.strange[sampled_index];
                sampled_params.charge   = particle_system.charge[sampled_index];
                sampled_params.sign     = particle_system.theta[sampled_index];

                double sampled_pLRF[4] = {0.0, 0.0, 0.0, 0.0};
                sample_momentum(sampled_params, sampled_pLRF, generator_momentum);


                double weight_flux = std::max(0.0, lrf.dsigma_t_lrf * sampled_pLRF[0] 
                                                 - lrf.dsigma_x_lrf * sampled_pLRF[1] 
                                                 - lrf.dsigma_y_lrf * sampled_pLRF[2]
                                                 - lrf.dsigma_z_lrf * sampled_pLRF[3])
                                                 /(lrf.dsigma_magnitude * sampled_pLRF[0]);

                double feq = 0.0;
                double delta_f = 0.0;
                ///\todo add feq and delta_f calculation
                //double weight_visc = 0.5*(1.+ delta_f/feq);
                double weight_visc = 0.5;
                //boost to lab frame
                lrf.boost_momentum_to_lab(tau_squared, sampled_pLRF);
                bool add_particle = (std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_keep) < (weight_flux*weight_visc));
                double E, pz, yp;
                if(D == 2 && add_particle)
                {
                  double random_number = std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_rapidity);
                  yp = y_max*(2.0*random_number - 1.0);

                  double sinhy = sinh(yp);
                  double coshy = sqrt(1.0 + sinhy * sinhy);

                  double ptau = lrf.pLab_tau;
                  double tau_pn = tau * lrf.pLab_eta;
                  double mT = sqrt(mass*mass  + lrf.pLab_x * lrf.pLab_x + lrf.pLab_y * lrf.pLab_y);

                  sinheta = (ptau*sinhy - tau_pn*coshy) / mT;
                  double eta = asinh(sinheta);
                  cosheta = sqrt(1.0 + sinheta * sinheta);

                  pz = mT * sinhy;
                  E = mT * coshy;
                  //if pz>5.0, print inf
                    //if(abs(pz) > 200.0){
                    //    std::cout << "Sampled id: " << pid << " name: " << name << std::endl;
                    //    std::cout << "pz: " << pz << " E: " << E << " yp: " << yp << " mT: " << mT << " ptau: " << ptau << " tau_pn: " << tau_pn << "px: " << lrf.pLab_x << " py: " << lrf.pLab_y << std::endl;
                    //    std::cout << "sinhy: " << sinhy << " coshy: " << coshy << " sinheta: " << sinheta << " cosheta: " << cosheta << std::endl;
                    //    std::cout << "mass: " << mass << " eta: " << eta << "peta: " << lrf.pLab_eta << std::endl;
                    //}
                  //cout << pz << "\t" << tau_pn * cosheta  +  ptau * sinheta << endl;
                }
                else
                {
                  pz = lrf.pLab_eta;
                  E = sqrt(mass * mass + lrf.pLab_x * lrf.pLab_x + lrf.pLab_y * lrf.pLab_y + pz * pz);
                  yp = 0.5 * log((E + pz) / (E - pz));
                  double eta =  surface.eta[icell];
                }

                double t = tau;
                double z = surface.eta[icell];
                if(add_particle){
                    double x = surface.x[icell];
                    double y = surface.y[icell];
                    Particle sampled_particle(pid, mass, E, lrf.pLab_x, lrf.pLab_y, pz, t, x, y, z);
                    sampled_particles[isample].push_back(sampled_particle);
                }


            }
        }
        auto momentum_sampling_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_momentum_sampling = momentum_sampling_end - momentum_sampling_start;
        t_momentum += elapsed_momentum_sampling.count();



        if (icell % progress_step == 0) {
            std::cout << "Progress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "% completed" << std::endl;
        }
    }
    auto end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_total = end_total - start_total;

    t_total += elapsed_total.count();

    std::cout << "Total tries: " << tries << std::endl;
    std::cout << "Total accepted: " << accepted << std::endl;
    std::cout << "Fraction accepted: " << (double)accepted/(double)tries << std::endl;

    std::cout << "Total sampling time: " << t_total << " seconds" << std::endl;
    std::cout << "Integrator computation time: " << t_integrator << " seconds" << std::endl;
    std::cout << "Momentum sampling time: " << t_momentum << " seconds" << std::endl;
    std::cout << "Species sampling time: " << t_species << " seconds" << std::endl;

    std::cout << "Total number of particles sampled: " << Ntot << std::endl;
    std::cout << "Total number of particles from cells: " << Ncell << std::endl;

    std::cout << "Massless above: " << nabove << std::endl;
    std::cout << "Massive above: " << nabove_massive << std::endl;

    std::string filename = "sampled_particles.dat";
    std::cout << "Saving sampled particles to " << filename << std::endl;
    save_particles(filename);


}



void Sampler::sample_momentum(const ThermalParams& params, double pLRF[4], std::default_random_engine& generator_momentum) {
    
    double max_w;   

    double mbar = params.mbar;
    double T = params.T;
    double omega = (params.baryon * params.alphaB + params.charge * params.alphaQ + params.strange * params.alphaS);
    double theta = params.sign;

    double MCUT =  1.0067; 
    ///\todo Determine the correct parameter

    bool rejected = true;
    double mbar_squared = mbar * mbar;
    double phi_over_2pi;
    double pbar, Ebar, phi, costheta, weight;
    //pion sampling (considered massless)
    if(mbar <= MCUT){
        max_w = get_max_w(params);
        while(rejected)
        {
            tries++;
            

            //generate canonical generate in [0,1), we then subtract from 1 to get (0,1] to avoid log(0)
            double r1= 1. - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);
            double r2= 1. - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);
            double r3= 1. - std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum);

            double l1 = log(r1);
            double l2 = log(r2);
            double l3 = log(r3);


            double l1_plus_l2 = l1 + l2;

            pbar = - (l1 + l2 + l3);
            Ebar = sqrt(pbar * pbar  +  mbar_squared);

            phi_over_2pi = l1_plus_l2 * l1_plus_l2 / (pbar * pbar);
            costheta = (l1 - l2) / l1_plus_l2;
            //using exp(pbar) = 1/ (r1 * r2 * r3)
            weight = 1.0 / (exp(Ebar-omega) + theta) / max_w / (r1 * r2 * r3);
            // check pLRF acceptance
            if(std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum) < weight) break;
        } // rejection loop
        accepted++;
    }
    else{
        //heavy particles
        std::uniform_real_distribution<double> costheta_distribution(-1.0, nextafter(1.0, std::numeric_limits<double>::max()));

        max_w = get_max_w_massive(params);

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
            tries++;
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
            weight = pbar/Ebar * exponent / (exponent + theta) / max_w;


            //if(fabs(weight - 0.5) > 0.5) printf("Sample momemtum error: weight = %f out of bounds\n", weight);
            // check pLRF acceptance
            if(std::generate_canonical<double, std::numeric_limits<double>::digits>(generator_momentum) < weight) break;
        } // rejection loop
        accepted++;
    }


    double E = Ebar * T;
    double p = pbar * T;
    phi = phi_over_2pi * 2.0 * M_PI;          // azimuthal angle
    double sintheta = sqrt(1.0  -  costheta * costheta);    // sin(theta)
    double px_lrf = p * sintheta * cos(phi);
    double py_lrf = p * sintheta * sin(phi);
    double pz_lrf = p * costheta;

    pLRF[0] = E;
    pLRF[1] = px_lrf;
    pLRF[2] = py_lrf;
    pLRF[3] = pz_lrf;

}





bool Sampler::check_weight_region(double mbar, double omega) {
    double boundary_value;
    double mbar_threshold = 0.8551;
    if (mbar < mbar_threshold) {
        boundary_value = poly_neg[0] *pow(mbar, 4) + poly_neg[1] * pow(mbar, 3) + poly_neg[2] * pow(mbar, 2) + poly_neg[3] * mbar + poly_neg[4];
    } else {
        boundary_value = poly_pos[0] *pow(mbar, 4) + poly_pos[1] * pow(mbar, 3) + poly_pos[2] * pow(mbar, 2) + poly_pos[3] * mbar + poly_pos[4];
    }
    return omega>boundary_value;
}


bool Sampler::check_weight_region_massive(double mbar, double omega) {
    double boundary_value;
    double MCUT =  1.0067;
    ///\todo make a read in parameter maybe?
    if (mbar < MCUT) {
        boundary_value =  poly_neg_massive[0] * pow(mbar, 3) + poly_neg_massive[1] * pow(mbar, 2) + poly_neg_massive[2] * mbar + poly_neg_massive[3];
    } else if(mbar >= MCUT && mbar < 5.1191) {
        boundary_value = poly_pos_massive_1[0] * pow(mbar, 3) + poly_pos_massive_1[1] * pow(mbar, 2) + poly_pos_massive_1[2] * mbar + poly_pos_massive_1[3];
    }
    else{
        boundary_value = poly_pos_massive_2[0] * mbar + poly_pos_massive_2[1];
    }

    return omega>boundary_value;
}




double Sampler::get_max_w(const ThermalParams& params) {
    double max_w = 1.0;
    double mbar = params.mbar;
    double omega = (params.baryon * params.alphaB + params.charge * params.alphaQ + params.strange * params.alphaS);
    if(params.sign == -1.0){
        if(check_weight_region(mbar, omega)){
            nabove++;
            max_w = max_w_table.bilinear_interpolate(mbar, omega);
            //max_w = 1.;

        }
    }
    return max_w;
}

double Sampler::get_max_w_massive(const ThermalParams& params) {
    double max_w = 1.0;
    double mbar = params.mbar;
    double omega = (params.baryon * params.alphaB + params.charge * params.alphaQ + params.strange * params.alphaS);
    if(params.sign == -1.0){
        if(check_weight_region_massive(mbar, omega)){
            nabove_massive++;
            max_w = max_w_table_massive.bilinear_interpolate(mbar, omega);
            //max_w = 1.;
        }
    }
    return max_w;
}


void Sampler::save_particles(const std::string& filename) const {
    std::ofstream output_file(filename);
    for (size_t i = 0; i < sampled_particles.size(); i++) {
        output_file << "# event " << i << "\n";
        for (size_t j = 0; j < sampled_particles[i].size(); j++) {
            output_file << sampled_particles[i][j].pid << " "
                        << sampled_particles[i][j].t << " "
                        << sampled_particles[i][j].x << " "
                        << sampled_particles[i][j].y << " "
                        << sampled_particles[i][j].z << " "
                        << sampled_particles[i][j].mass << " "
                        << sampled_particles[i][j].E << " "
                        << sampled_particles[i][j].px << " "
                        << sampled_particles[i][j].py << " "
                        << sampled_particles[i][j].pz << " " << std::endl;
        }
        output_file << "# event " << i << " end" << "\n";
    }
    output_file.close();
}