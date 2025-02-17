#include "particle.h"

/// @brief Construct a Particle object with default values
Particle::Particle(int pid, std::string name, double mass, double width, 
                   double spin_degeneracy, double baryon, double strange, 
                   double charm, double bottom, double isospin_degeneracy, 
                   double charge, int decays, double theta)
    : pid(pid), name(name), mass(mass), width(width), spin_degeneracy(spin_degeneracy),
      baryon(baryon), strange(strange), charm(charm), bottom(bottom),
      isospin_degeneracy(isospin_degeneracy), charge(charge), 
      equilibrium_density(0.0), decays(decays), theta(theta) {}


/// @brief Display the particle's data
void Particle::display() const {
    std::cout << "Particle ID: " << pid << ", Name: " << name << ", Mass: " << mass
              << ", Width: " << width << ", Spin Degeneracy: " << spin_degeneracy 
              << ", Baryon: " << baryon << ", Strange: " << strange 
              << ", Charm: " << charm << ", Bottom: " << bottom 
              << ", Isospin Degeneracy: " << isospin_degeneracy 
              << ", Charge: " << charge << ", Decays: " << decays 
              << ", Equilibrium Density: " << equilibrium_density << std::endl;
}