#include "sampler.h"



// Constructor implementation
Sampler::Sampler(const Settings& settings): settings(settings) 
{
    // Initialize the random number generator with a seed based on the current time



    int seed_from_config = settings.get_int("seed");
    if (seed_from_config > 0) {
        sampler_seed = seed_from_config;
        std::cout << "Using user-defined seed: " << sampler_seed << std::endl;
    }
    else {
        sampler_seed =  static_cast<unsigned>(
        std::chrono::system_clock::now().time_since_epoch().count() );
        std::cout << "Using random seed: " << sampler_seed << std::endl;
    }
    std::string tables_path = settings.get_string("tables_path");
    max_w_table.load_from_file(tables_path + "/max_w_table.dat");
    max_w_table_massive.load_from_file(tables_path + "/max_w_table_massive.dat");
    deltaf_table.load_from_file(tables_path + "/grad_coeffs.dat", 14); 
    accepted = 0;
    tries = 0;
    nabove = 0;
    nabove_massive = 0;
    this->gen_poisson = std::mt19937(sampler_seed);
    this->gen_type    = std::mt19937(sampler_seed + 10000);
    this->gen_mom     = std::mt19937(sampler_seed + 20000);
    this->gen_keep    = std::mt19937(sampler_seed + 30000);
    this->gen_y       = std::mt19937(sampler_seed + 40000);
    this->gen_trim    = std::mt19937(sampler_seed + 50000);
    this->gen_pos     = std::mt19937(sampler_seed + 60000);
    pos_smearing = settings.has_key("position_smearing") ? settings.get_double("position_smearing") : 0.0;


    Nsamples = settings.get_int("samples");
    D = settings.get_int("dimension");
    y_max = settings.get_double("y_max");
    if (D==3) y_max = 0.5; // we use 2.*ymax, canceling this parameter for 3D (where there is a y size)
    
}

void Sampler::sample_unconstrained(ParticleSystem& all_particles,
                                   Surface& surface,
                                   const NumericalIntegrator& integrator)
{
    std::string coordinate_system = settings.get_string("coordinate_system");
    sampled_particles.resize(Nsamples);

    if (surface.npoints <= 0) throw std::runtime_error("surface.npoints is zero or negative");


    // --- Step 1: Compute cell integrals
    std::vector<double>& N_cell = surface.N_all_particles_cell;
    N_cell.resize(surface.npoints, 0.0);
    double N_skip = 0.0;
    std::cout << "Calculating surface cell integrals ..." << std::endl;
    for (int icell = 0; icell < surface.npoints; ++icell) {

        //print percentage of progress
        { int step = std::max(1, surface.npoints / 20);
          if (icell % step == 0) {
            std::cout << "Integration progress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "% completed" << std::endl;
          }
        }

        double udsigma = surface.ut[icell] * surface.dsigma_t[icell]
                       + surface.ux[icell] * surface.dsigma_x[icell]
                       + surface.uy[icell] * surface.dsigma_y[icell]
                       + surface.ueta[icell] * surface.dsigma_eta[icell];

        if (udsigma <= 0.0){
            N_skip += 1.;
            continue;
        }

        double T   = surface.T[icell];
        double muB = surface.muB[icell];
        double muS = surface.muS[icell];
        double muQ = surface.muQ[icell];

        double N = all_particles.calculate_particle_number(T, muB, muQ, muS, integrator);
        N_cell[icell] = N;
    }
    std::cout << "Finished calculating surface cell integrals." << std::endl;
    // --- Step 2: Sample events ---
    double total_B = 0, total_S = 0, total_Q = 0;
    std::cout << "Sampling events with Nsamples = " << Nsamples << std::endl;
    for (int isample = 0; isample < Nsamples; ++isample) {
        //print percentage of samples sampled 
        int progress_step = std::max(1, Nsamples / 20);
        if (isample % progress_step == 0) {
            std::cout << "Sampling progress: " << (static_cast<double>(isample) / Nsamples) * 100 << "% completed" << std::endl;
        }
        std::vector<Particle> event = sample_fixed_yield_from_surface(
            all_particles, surface, N_cell, -1, coordinate_system);

        // Accumulate total BSQ
        for (const auto& p : event) {
            total_B += p.baryon;
            total_S += p.strange;
            total_Q += p.charge;
        }

        sampled_particles[isample] = std::move(event);
    }

    // --- Step 3: Charge check and save ---
    std::cout << "Finished sampling events." << std::endl;
    check_total_charge_average(total_B, total_S, total_Q, Nsamples);

    std::string filename = settings.get_string("output_file");
    std::cout << "Saving sampled particles to " << filename << std::endl;
    //print the percentage of skipped cells
    std::cout << "Skipped " << N_skip << " cells out of " << surface.npoints << " cells." << std::endl;
    std::cout << "Percentage of skipped cells: " 
              << (N_skip / surface.npoints) * 100 << "%" << std::endl;
    save_particles(filename);
}



void Sampler::sample_momentum(const ThermalParams& params, double pLRF[4], std::mt19937& generator_momentum) {
    
    double max_w;   

    double mbar = params.mbar;
    double T = params.T;
    double omega = (params.baryon * params.alphaB + params.charge * params.alphaQ + params.strange * params.alphaS);
    double theta = params.sign;

    double MCUT =  1.67; 
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
    FILE* f = fopen(filename.c_str(), "w");
    if (!f) {
        perror("Error opening file for writing");
        return;
    }

    for (size_t i = 0; i < sampled_particles.size(); ++i) {
        fprintf(f, "# event %zu\n", i);

        for (const auto& p : sampled_particles[i]) {
            // write each particle as a single line
            fprintf(f, "%d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
                    p.pid, p.t, p.x, p.y, p.z, p.mass, p.E, p.px, p.py, p.pz);
        }

        fprintf(f, " end %zu \n", i);
    }

    fclose(f);
}


void Sampler::conserved_charge_sampling(ParticleSystem& particle_system, Surface& surface, const NumericalIntegrator& integrator) {
    std::string coordinate_system = settings.get_string("coordinate_system");
    sampled_particles.resize(Nsamples);

    int netB = surface.Btotal;
    int netS = surface.Stotal;
    int netQ = surface.Qtotal;

    if (surface.npoints <= 0) throw std::runtime_error("surface.npoints is zero or negative");



    // --- Initialize BSQ groups ---
    ParticleSystem baryons, antibaryons, strange_mesons_sminus, strange_mesons_splus;
    ParticleSystem charged_mesons_qplus, charged_mesons_qminus, neutral_mesons;

    for (int i = 0; i < particle_system.nparticles; ++i) {
        double B = particle_system.baryon[i];
        double S = particle_system.strange[i];
        double Q = particle_system.charge[i];
        int pid = particle_system.pid[i];

        if (B == 1) particle_system.copy_particle(i, baryons);
        else if (B == -1) particle_system.copy_particle(i, antibaryons);
        else if (B == 0 && S == -1) particle_system.copy_particle(i, strange_mesons_sminus);
        else if (B == 0 && S == +1) particle_system.copy_particle(i, strange_mesons_splus);
        else if (B == 0 && S == 0 && Q == +1) particle_system.copy_particle(i, charged_mesons_qplus);
        else if (B == 0 && S == 0 && Q == -1) particle_system.copy_particle(i, charged_mesons_qminus);
        else if (B == 0 && S == 0 && Q == 0 && pid != 22) particle_system.copy_particle(i, neutral_mesons);
    }

    std::cout << " Calculating surface cell integrals ..." << std::endl;

    surface.N_baryons_cell.resize(surface.npoints, 0.0);
    surface.N_antibaryons_cell.resize(surface.npoints, 0.0);
    surface.N_strange_mesons_sminus_cell.resize(surface.npoints, 0.0);
    surface.N_strange_mesons_splus_cell.resize(surface.npoints, 0.0);
    surface.N_charged_mesons_qplus_cell.resize(surface.npoints, 0.0);
    surface.N_charged_mesons_qminus_cell.resize(surface.npoints, 0.0);
    surface.N_neutral_mesons_cell.resize(surface.npoints, 0.0);

    for (int icell = 0; icell < surface.npoints; ++icell) {
        double T = surface.T[icell];
        double muB = surface.muB[icell];
        double muS = surface.muS[icell];
        double muQ = surface.muQ[icell];
        { int step = std::max(1, surface.npoints / 20);
          if (icell % step == 0) {
            std::cout << "Integration rogress: " << (static_cast<double>(icell) / surface.npoints) * 100 << "% completed" << std::endl;
          }
        }
        double udsigma = surface.ut[icell] * surface.dsigma_t[icell]
                       + surface.ux[icell] * surface.dsigma_x[icell]
                       + surface.uy[icell] * surface.dsigma_y[icell]
                       + surface.ueta[icell] * surface.dsigma_eta[icell];
        if (udsigma <= 0.0) std::cout << "Skipping cell " << icell << " due to non-positive u dot dsigma: " << udsigma << std::endl;
        if (udsigma <= 0.0) continue;

        surface.N_baryons_cell[icell]               = baryons.calculate_particle_number(T, muB, muQ, muS, integrator);
        surface.N_antibaryons_cell[icell]           = antibaryons.calculate_particle_number(T, muB, muQ, muS, integrator);
        surface.N_strange_mesons_sminus_cell[icell] = strange_mesons_sminus.calculate_particle_number(T, muB, muQ, muS, integrator);
        surface.N_strange_mesons_splus_cell[icell]  = strange_mesons_splus.calculate_particle_number(T, muB, muQ, muS, integrator);
        surface.N_charged_mesons_qplus_cell[icell]  = charged_mesons_qplus.calculate_particle_number(T, muB, muQ, muS, integrator);
        surface.N_charged_mesons_qminus_cell[icell] = charged_mesons_qminus.calculate_particle_number(T, muB, muQ, muS, integrator);
        surface.N_neutral_mesons_cell[icell]        = neutral_mesons.calculate_particle_number(T, muB, muQ, muS, integrator);
    }

    std::cout << "Finished calculating surface cell integrals." << std::endl;
    std::cout << "Sampling events with Nsamples = " << Nsamples << std::endl;
    for (int isample = 0; isample < Nsamples; ++isample) {
        std::vector<Particle> event;
        //print percentage of samples sampled 
        int progress_step = std::max(1, Nsamples / 20);
        if (isample % progress_step == 0) {
            std::cout << "Sampling progress: " << (static_cast<double>(isample) / Nsamples) * 100 << "% completed" << std::endl;
        }

        // === Step 1: Baryons ===
        std::vector<Particle> baryon_sample;
        // set N_baryons to numeric limits, so we always sample baryons even if netB < 0
        int N_baryons = -9999999;
        while (N_baryons < netB) {
            std::cout << "Sampling baryons..." << std::endl;
            baryon_sample = sample_fixed_yield_from_surface(baryons, surface, surface.N_baryons_cell, -1, coordinate_system);
            N_baryons = net_baryon(baryon_sample);
        }
        std::cout << "Sampled " << N_baryons << " baryons." << std::endl;

        int N_antibaryons = N_baryons - netB;
        std::vector<Particle> antibaryon_sample = sample_fixed_yield_from_surface(antibaryons, surface, surface.N_antibaryons_cell, N_antibaryons, coordinate_system);
        std::cout << "Sampled " << antibaryon_sample.size() << " antibaryons." << std::endl;

        std::vector<Particle> s_mesons = sample_fixed_yield_from_surface(strange_mesons_sminus, surface, surface.N_strange_mesons_sminus_cell, -1, coordinate_system);
        std::cout << "Sampled " << s_mesons.size() << " strange mesons (S = -1)." << std::endl;

        int NS_baryons = net_strangeness(baryon_sample) + net_strangeness(antibaryon_sample);
        int NS_strange = net_strangeness(s_mesons);
        int NS_needed = NS_baryons + NS_strange - netS;
        //take the absolute value
        NS_needed = std::abs(NS_needed);
        //std::cout << "NS_needed = " << NS_needed << std::endl;
        //std::cout << "NS_baryons = " << NS_baryons << ", NS_strange = " << NS_strange << std::endl;
        //std::cout << "netS = " << netS << std::endl;
        //std::cout << "Baryon strangeness: " << net_strangeness(baryon_sample) << std::endl;
        //std::cout << "Antibaryon strangeness: " << net_strangeness(antibaryon_sample) << std::endl;
        std::vector<Particle> anti_strange_mesons = sample_fixed_yield_from_surface(strange_mesons_splus, surface, surface.N_strange_mesons_splus_cell, NS_needed, coordinate_system);
        std::cout << "Sampled " << anti_strange_mesons.size() << " strange mesons (S = +1)." << std::endl;

        std::vector<Particle> pos_mesons;
        int Q_baryons = net_charge(baryon_sample) + net_charge(antibaryon_sample);
        int Q_strange = net_charge(s_mesons) + net_charge(anti_strange_mesons);
        int Q_remain = netQ - Q_baryons - Q_strange;
        
        //set to numeric limits, so we always sample electric charge even if Q_remain < 0
        int N_qplus = -9999999;
        while (N_qplus < Q_remain) {
            pos_mesons = sample_fixed_yield_from_surface(charged_mesons_qplus, surface, surface.N_charged_mesons_qplus_cell, -1, coordinate_system);
            N_qplus = net_charge(pos_mesons);
        }
        std::cout << "Sampled " << pos_mesons.size() << " charged mesons (Q = +1)." << std::endl;

        int Q_BQ = net_charge(baryon_sample) + net_charge(antibaryon_sample);
        int Q_SM = net_charge(s_mesons) + net_charge(anti_strange_mesons);
        int Q_pos = net_charge(pos_mesons);
        int N_neg = Q_BQ + Q_SM + Q_pos - netQ;

        std::vector<Particle> neg_mesons = sample_fixed_yield_from_surface(charged_mesons_qminus, surface, surface.N_charged_mesons_qminus_cell, N_neg, coordinate_system);
        std::cout << "Sampled " << neg_mesons.size() << " charged mesons (Q = -1)." << std::endl;

        std::vector<Particle> neutrals = sample_fixed_yield_from_surface(neutral_mesons, surface, surface.N_neutral_mesons_cell, -1, coordinate_system);
        std::cout << "Sampled " << neutrals.size() << " neutral mesons (Q = 0)." << std::endl;

        // Final event particle list
        event.insert(event.end(), baryon_sample.begin(), baryon_sample.end());
        event.insert(event.end(), antibaryon_sample.begin(), antibaryon_sample.end());
        event.insert(event.end(), s_mesons.begin(), s_mesons.end());
        event.insert(event.end(), anti_strange_mesons.begin(), anti_strange_mesons.end());
        event.insert(event.end(), pos_mesons.begin(), pos_mesons.end());
        event.insert(event.end(), neg_mesons.begin(), neg_mesons.end());
        event.insert(event.end(), neutrals.begin(), neutrals.end());

        // === New: Print totals ===
        int B_total = net_baryon(event);
        int S_total = net_strangeness(event);
        int Q_total = net_charge(event);
        std::cout << "[Event Totals] B = " << B_total << ", S = " << S_total << ", Q = " << Q_total << std::endl;

        std::cout << "[Event Subgroup Info] Strangeness in B/B̄ = " << NS_baryons << std::endl;
        std::cout << "[Event Subgroup Info] Charge in B/B̄ = " << Q_BQ << ", in S/S̄ = " << Q_SM << std::endl;

        sampled_particles[isample] = std::move(event);
    }
    std::cout << "Finished sampling events." << std::endl;
    std::string filename = settings.get_string("output_file");
    std::cout << "Saving sampled particles to " << filename << std::endl;
    save_particles(filename);
}


double Sampler::net_baryon(const std::vector<Particle>& plist) {
    double netB = 0;
    for (const auto& p : plist) netB += p.baryon;
    return netB;
}

double Sampler::net_strangeness(const std::vector<Particle>& plist) {
    double netS = 0;
    for (const auto& p : plist) netS += p.strange;
    return netS;
}

double Sampler::net_charge(const std::vector<Particle>& plist) {
    double netQ = 0;
    for (const auto& p : plist) netQ += p.charge;
    return netQ;
}

void Sampler::check_total_charge_average(double totalB, double totalS, double totalQ, int Nsamples) {
    double avgB = totalB / static_cast<double>(Nsamples);
    double avgS = totalS / static_cast<double>(Nsamples);
    double avgQ = totalQ / static_cast<double>(Nsamples);

    std::cout << "[check_total_charge_average] Avg B: " << avgB
              << ", Avg S: " << avgS
              << ", Avg Q: " << avgQ << std::endl;
}



std::vector<Particle> Sampler::sample_fixed_yield_from_surface(
    const ParticleSystem& group,
    const Surface& surface,
    const std::vector<double>& N_cell_vector,
    int required,
    const std::string& coord)
{
    std::vector<Particle> result;

    std::discrete_distribution<int> type_dist(group.particle_species_number.begin(), group.particle_species_number.end());
    bool fixed_yield = (required >= 0);
    int target_yield = fixed_yield ? required : 0;

    int max_attempts = fixed_yield ? 1e4 : 1;  // repeat sampling in fixed-yield mode only
    int attempt = 0;

    while ((fixed_yield && static_cast<int>(result.size()) < target_yield && attempt < max_attempts) ||
           (!fixed_yield && attempt < 1))
    {
        ++attempt;

        for (int icell = 0; icell < surface.npoints; ++icell) {
            if (fixed_yield && static_cast<int>(result.size()) >= target_yield)
                break;

            double T = surface.T[icell];
            double muB = surface.muB[icell];
            double muS = surface.muS[icell];
            double muQ = surface.muQ[icell];
            double tau = surface.tau[icell];
            double s0_cell = surface.s[icell];
            double tau_squared = tau * tau;

            double udsigma = surface.ut[icell] * surface.dsigma_t[icell] 
                           + surface.ux[icell] * surface.dsigma_x[icell] 
                           + surface.uy[icell] * surface.dsigma_y[icell] 
                           + surface.ueta[icell] * surface.dsigma_eta[icell];
            if (udsigma <= 0.0) continue;

            double N_tot_cell = N_cell_vector[icell];

            LRF lrf(coord,
                    surface.ut[icell], surface.ux[icell], surface.uy[icell], surface.ueta[icell],
                    surface.dsigma_t[icell], surface.dsigma_x[icell],
                    surface.dsigma_y[icell], surface.dsigma_eta[icell], tau);

            lrf.boost_dsigma_to_lrf(tau_squared);
            lrf.compute_dsigma_magnitude();
            N_tot_cell *= 2.0 * y_max * lrf.dsigma_magnitude;

            if (N_tot_cell <= 0.0) continue;

            std::poisson_distribution<int> poisson_hadrons(N_tot_cell);
            int N_hadrons = poisson_hadrons(gen_poisson);

            DissipativeParams diss_params;
            diss_params.shv_tt = surface.shv_tt[icell];
            diss_params.shv_tx = surface.shv_tx[icell];
            diss_params.shv_ty = surface.shv_ty[icell];
            diss_params.shv_teta = surface.shv_teta[icell];
            diss_params.shv_xx = surface.shv_xx[icell];
            diss_params.shv_xy = surface.shv_xy[icell];
            diss_params.shv_xeta = surface.shv_xeta[icell];
            diss_params.shv_yy = surface.shv_yy[icell]; 
            diss_params.shv_yeta = surface.shv_yeta[icell];
            diss_params.shv_etaeta = surface.shv_etaeta[icell];
            diss_params.bulk = surface.bulk[icell]; 
            diss_params.q_B0 = surface.diff_B0[icell];
            diss_params.q_Q0 = surface.diff_Q0[icell];
            diss_params.q_S0 = surface.diff_S0[icell];
            diss_params.q_Bx = surface.diff_Bx[icell];
            diss_params.q_Qx = surface.diff_Qx[icell];
            diss_params.q_Sx = surface.diff_Sx[icell];
            diss_params.q_By = surface.diff_By[icell];
            diss_params.q_Qy = surface.diff_Qy[icell];
            diss_params.q_Sy = surface.diff_Sy[icell];
            diss_params.q_Beta = surface.diff_Beta[icell];
            diss_params.q_Qeta = surface.diff_Qeta[icell];
            diss_params.q_Seta = surface.diff_Seta[icell];

            for (int i = 0; i < N_hadrons; ++i) {
                if (fixed_yield && static_cast<int>(result.size()) >= target_yield)
                    break;

                int sampled_index = type_dist(gen_type);
                int pid = group.pid[sampled_index];
                double mass = group.mass[sampled_index];
                double B = group.baryon[sampled_index];
                double S = group.strange[sampled_index];
                double Q = group.charge[sampled_index];

                ThermalParams sampled_params;
                sampled_params.T = T;
                sampled_params.alphaB = muB / T;
                sampled_params.alphaQ = muQ / T;
                sampled_params.alphaS = muS / T;
                sampled_params.mbar = mass / T;
                sampled_params.baryon = B;
                sampled_params.strange = S;
                sampled_params.charge = Q;
                sampled_params.sign = group.theta[sampled_index];
                
                double E_cell = surface.E[icell];
                double P_cell = surface.P[icell];
                double sampled_pLRF[4] = {0.0};
                sample_momentum(sampled_params, sampled_pLRF, gen_mom);

                double flux = std::max(0.0, lrf.dsigma_t_lrf * sampled_pLRF[0]
                                             - lrf.dsigma_x_lrf * sampled_pLRF[1]
                                             - lrf.dsigma_y_lrf * sampled_pLRF[2]
                                             - lrf.dsigma_z_lrf * sampled_pLRF[3])
                              / (lrf.dsigma_magnitude * sampled_pLRF[0]);

                double feq = group.spin_degeneracy[sampled_index]/(std::exp((sampled_pLRF[0]
                                    - B * muB
                                    - Q * muQ
                                    - S * muS) / T) + group.theta[sampled_index]);
                double delta_f = feq*(1.-group.theta[sampled_index]*feq/group.spin_degeneracy[sampled_index])
                                    *df_corrections(settings,deltaf_table ,lrf, tau_squared, sampled_pLRF, T, E_cell,s0_cell, P_cell, mass,sampled_params ,diss_params);
                         
                //normalize delta_f
                if (settings.get_bool("normalize_deltaf")) {
                    if (delta_f >  feq) delta_f =  feq;
                    if (delta_f < -feq) delta_f = -feq;
                }
                
                double weight_visc = (1. + delta_f/feq);

                const double u_keep =
                    std::generate_canonical<double, std::numeric_limits<double>::digits>(gen_keep);
                const double u_y =
                    std::generate_canonical<double, std::numeric_limits<double>::digits>(gen_y);  

                bool add_particle = u_keep < (flux * weight_visc);
                if (!add_particle) continue;

                lrf.boost_momentum_to_lab(tau_squared, sampled_pLRF);

                double sinheta = sinh(surface.eta[icell]);
                double cosheta = sqrt(1.0 + sinheta * sinheta);

                double px = lrf.pLab_x;
                double py = lrf.pLab_y;
                double pz = 0.0;
                double E = 0.0;
                double yp = 0.0;
                double t,z;
                if (D == 2) {
                    yp = y_max * (2.0 * u_y - 1.0);
                    double sinhy = sinh(yp);
                    double coshy = sqrt(1.0 + sinhy * sinhy);

                    double mT = sqrt(mass * mass + px * px + py * py);
                    double ptau = lrf.pLab_tau;
                    double tau_pn = tau * lrf.pLab_eta;

                    sinheta = (ptau * sinhy - tau_pn * coshy) / mT;
                    cosheta = sqrt(1.0 + sinheta * sinheta);
                    pz = mT * sinhy;
                    E  = mT * coshy;

                } else if (D == 3) {
                    if (coord == "cartesian") {
                        pz = lrf.pLab_eta;
                    } else {
                        pz = tau * lrf.pLab_eta * cosheta + lrf.pLab_tau * sinheta;
                    }
                    E = sqrt(mass * mass + px * px + py * py + pz * pz);
                    yp = 0.5 * log((E + pz) / (E - pz));
                }

                double x = surface.x[icell];
                double y = surface.y[icell];

                if(coord == "cartesian"){
                    t = tau;
                    z = surface.eta[icell];
                }
                else{
                    t = tau * cosheta;
                    z = tau * sinheta;
                }

                if (pos_smearing > 0.0) {
                    std::uniform_real_distribution<double> smear(-pos_smearing, pos_smearing);
                    x += smear(gen_pos);
                    y += smear(gen_pos);
                    z += smear(gen_pos);
                }

                Particle p(pid, mass, E, px, py, pz, t, x, y, z,
                           static_cast<int>(B), static_cast<int>(S), static_cast<int>(Q));
                result.push_back(p);
            }
        }
    }

    // Trim if too many
    if (fixed_yield && static_cast<int>(result.size()) > target_yield) {
        std::shuffle(result.begin(), result.end(), gen_trim);
        result.resize(target_yield);
    }

    return result;
}
