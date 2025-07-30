#include "particle_system.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


/// @brief Default constructor for ParticleSystem
ParticleSystem::ParticleSystem() : nparticles(0) {}

/// @brief Read the particle list from a file and populate the ParticleSystem
/// @details This function reads particle data from a file and populates 
/// the ParticleSystem with Particle objects. It skips lines for decay products, 
/// ensuring that only the main particles are stored in the system.
/// @param filename The path to the file containing particle data
void ParticleSystem::read_particle_list(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open particle list file " << filename << std::endl;
        return;
    }

    // Temporary variables to hold the data read from the file
    int pid, decays;
    std::string name;
    double mass, width, spin_degeneracy, baryon, strange, charm, bottom, isospin_degeneracy, charge;
    int baryons = 0;
    int mesons = 0;
    int antibaryons = 0;
    while (true) {
        std::string line;
        if (!std::getline(file, line)) {
            break; // Exit if we reach the end of the file
        }

        // Skip empty lines
        if (line.empty()) continue;

        // Read particle data from the current line
        std::istringstream iss(line);
        if (!(iss >> pid >> name >> mass >> width >> spin_degeneracy >> baryon >> strange >> charm
                 >> bottom >> isospin_degeneracy >> charge >> decays)) {
            std::cerr << "Error reading particle data from line: " << line << std::endl;
            continue; // Skip this line and move to the next
        }
        double theta;
        if (static_cast<int>(baryon)  % 2 == 0) {  // Check if baryon number is even
            theta = -1.;
        } else {
           theta = 1.;
        }

        // Construct the Particle object and add to the system
        this->pid.push_back(pid);
        this->name.push_back(name);
        this->mass.push_back(mass);
        this->width.push_back(width);
        this->spin_degeneracy.push_back(spin_degeneracy);
        this->baryon.push_back(baryon);
        this->strange.push_back(strange);
        this->charm.push_back(charm);
        this->bottom.push_back(bottom);
        this->isospin_degeneracy.push_back(isospin_degeneracy);
        this->charge.push_back(charge);
        this->decays.push_back(decays);
        this->theta.push_back(theta);
        this->particle_species_number.push_back(0.0);
        this->equilibrium_density.push_back(0.0);
        nparticles++;

        if (baryon == 0) {
            mesons++;
        } else if (baryon > 0) {
            baryons++;

            // add antibaryons to vectors
            this->pid.push_back(-pid);
            this->name.push_back("anti-" + name);
            this->mass.push_back(mass);
            this->width.push_back(width);
            this->spin_degeneracy.push_back(spin_degeneracy);
            this->baryon.push_back(-baryon);
            this->strange.push_back(-strange);
            this->charm.push_back(-charm);
            this->bottom.push_back(-bottom);
            this->isospin_degeneracy.push_back(isospin_degeneracy);
            this->charge.push_back(-charge);
            this->decays.push_back(decays);
            this->theta.push_back(theta);
            this->particle_species_number.push_back(0.0);
            this->equilibrium_density.push_back(0.0);
            nparticles++;
            antibaryons++;
        }


        // Skip decay product lines
        for (int i = 0; i < decays; ++i) {
            std::getline(file, line); // Skip the decay product line
        }
    }

    file.close();
    if(baryons != antibaryons) printf("Error: (anti)baryons not paired correctly\n");

    printf("\nThere are %d mesons, %d baryons and %d antibaryons\n\n", mesons, baryons, antibaryons);
    
}


///@brief Display all particles in the system
///@details This function iterates over all particles in the system and
/// calls the display method of each Particle object to print its details
void ParticleSystem::display_all_particles() const {
    for (int i = 0; i < nparticles; ++i) {
        std::cout << "Particle ID: " << pid[i] << ", Name: " << name[i] << ", Mass: " << mass[i]
          << ", Width: " << width[i] << ", Spin Degeneracy: " << spin_degeneracy[i] 
          << ", Baryon: " << baryon[i] << ", Strange: " << strange[i] 
          << ", Charm: " << charm[i] << ", Bottom: " << bottom[i] 
          << ", Isospin Degeneracy: " << isospin_degeneracy[i] 
          << ", Charge: " << charge[i] << ", Decays: " << decays[i] 
          << ", Equilibrium Density: " << equilibrium_density[i] << std::endl;
    }
}


/// @brief Copy a particle from this ParticleSystem to another ParticleSystem
/// @details This function copies the properties of a particle at a given index
/// from this ParticleSystem to the destination ParticleSystem.
/// @param index The index of the particle to copy
/// @param dest The destination ParticleSystem to which the particle will be copied
void ParticleSystem::copy_particle(int index, ParticleSystem& dest) const {
    if (index < 0 || index >= static_cast<int>(pid.size())) {
        throw std::out_of_range("ParticleSystem::copy_particle index out of bounds");
    }
    dest.pid.push_back(pid[index]);
    dest.name.push_back(name[index]);
    dest.mass.push_back(mass[index]);
    dest.width.push_back(width[index]);
    dest.spin_degeneracy.push_back(spin_degeneracy[index]);
    dest.baryon.push_back(baryon[index]);
    dest.strange.push_back(strange[index]);
    dest.charm.push_back(charm[index]);
    dest.bottom.push_back(bottom[index]);
    dest.isospin_degeneracy.push_back(isospin_degeneracy[index]);
    dest.charge.push_back(charge[index]);
    dest.decays.push_back(decays[index]);
    dest.theta.push_back(theta[index]);
    dest.particle_species_number.push_back(0.0);
    dest.equilibrium_density.push_back(0.0);
    dest.nparticles++;
}


// @brief  Calculate the equilibrium density of all particles
/// @details Calculate the equilibrium density of the particle species
/// in a given freeze out cell. The equilibrium density integral is 
/// int_{0}^{\infty} \frac{pLRF^2_bar}{e^{(E_bar - X \mu/T)} + theta} dpLRF_bar, which we write
/// as int_{0}^{\infty} pLRF_bar e^{pLRF_bar} integrand(pLRF_bar) dpLRF_bar, to be able to use the
/// generalized Gauss-Laguerre quadrature, where the integrand absorved the
/// an exponential and a factor of pLRF.
/// @param surface The surface object containing the thermodynamic data
/// @param integrator The numerical integrator object to perform the integration
/// @return The total yield of all particles
double ParticleSystem::calculate_particle_number(double T, double muB, double muQ, double muS, const NumericalIntegrator& integrator) {
    double particle_density = 0.0;
    for (int i = 0; i < nparticles; i++) {
        ThermalParams params;
        params.mbar = mass[i]/T; 
        params.alphaB = muB/T;
        params.alphaQ = muQ/T;
        params.alphaS = muS/T;
        params.sign = theta[i];
        params.T = T;
        params.baryon = baryon[i];
        params.charge = charge[i];
        params.strange = strange[i];
        //integrate over pTok, g
        double integral = integrator.gauss_quadrature(
            get_equilibrium_density,
            params,
            integrator.roots_gauss_laguerre[1],
            integrator.weights_gauss_laguerre[1]
        );
        double particle_equilibrium_density = spin_degeneracy[i] * pow(T,3)* integral / (2.* pow(M_PI, 2.0) * pow(HBARC, 3.0));

        if(mass[i] == 0.0) particle_equilibrium_density = 0.0; //dont sample photons
        
        equilibrium_density[i] = particle_equilibrium_density;
        particle_density += particle_equilibrium_density;
        particle_species_number[i] = 2. * particle_equilibrium_density;
        //std::cout << "Particle ID: " << pid[i] << ", Equilibrium Density: " << particle_equilibrium_density << ", Particle Species Number: " << particle_species_number[i] << std::endl;
        //check for positive pion 



    }
    /// @todo check this 2. factor
    return particle_density*2.;
}

