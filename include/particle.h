#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <iostream>
#include <cmath>

#include "numerical_integrator.h"
#include "surface.h"


class Particle {
public:
    int pid;                 // Particle ID
    std::string name;        // Name of the particle (e.g., π⁻, π⁺, π⁰)
    double mass;             // Mass of the particle
    double width;            // Width (if applicable)
    double spin_degeneracy;  // Spin degeneracy
    double baryon;           // Baryon number
    double strange;          // Strange quark number
    double charm;            // Charm quark number
    double bottom;           // Bottom quark number
    double isospin_degeneracy; // Isospin degeneracy
    double charge;           // Electric charge
    double equilibrium_density; // Equilibrium density
    int decays;              // Number of decay channels
    double theta;            // Statistics: -1 for boson, 1 for fermion
    double particle_number;   // particle number for a given surface point

    // Constructor with default values
    Particle(int pid = 0, std::string name = "", double mass = 0.0, double width = 0.0, 
             double spin_degeneracy = 1.0, double baryon = 0.0, double strange = 0.0, double charm = 0.0,
             double bottom = 0.0, double isospin_degeneracy = 1.0, double charge = 0.0, int decays = 0, double theta = 1.0);

    // Method to display the particle's data
    void display() const;
    


};

#endif // PARTICLE_H  