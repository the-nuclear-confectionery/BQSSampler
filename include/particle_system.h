#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <vector>
#include <fstream>
#include "particle.h"



class ParticleSystem {


public:
    int nparticles;

    std::vector<int> pid;                 // Particle ID
    std::vector<std::string> name;        // Name of the particle (e.g., π⁻, π⁺, π⁰)
    std::vector<double> mass;             // Mass of the particle
    std::vector<double> width;            // Width (if applicable)
    std::vector<double> spin_degeneracy;  // Spin degeneracy
    std::vector<double> baryon;           // Baryon number
    std::vector<double> strange;          // Strange quark number
    std::vector<double> charm;            // Charm quark number
    std::vector<double> bottom;           // Bottom quark number
    std::vector<double> isospin_degeneracy; // Isospin degeneracy
    std::vector<double> charge;           // Electric charge
    std::vector<double> equilibrium_density; // Equilibrium density
    std::vector<int> decays;              // Number of decay channels
    std::vector<double> theta;            // Statistics: -1 for boson, 1 for fermion
    std::vector<double> particle_species_number;   // particle number for a given surface point   


    // Constructor
    ParticleSystem();

    // Method to read particles from a file, skipping decay channels
    void read_particle_list(const std::string& filename);

    // Method to display all particles
    void display_all_particles() const;

    double calculate_particle_number(double T, double muB, double muQ, double muS, const NumericalIntegrator& integrator);


    void copy_particle(int index, ParticleSystem& dest) const;
};


#endif // PARTICLE_SYSTEM_H