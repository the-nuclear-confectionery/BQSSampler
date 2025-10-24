#include "surface.h"
#include "particle.h"
#include "numerical_integrator.h"
#include "sampler.h"
#include "particle_system.h"
#include "globals.h"
#include "settings.h"
#include <iostream>

int main(int argc, char** argv) {
    // Require user to provide config path
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>\n";
        return 1;
    }

    std::string config_file = argv[1];
    Settings config(config_file);
    // Create a Surface object
    Surface surface(config.get_string("input_file"), config);
    
    // check if using cartesian and 2D at same time
    if (config.get_string("coordinate_system") == "cartesian" && config.get_int("dimension") == 2) {
        std::cerr << "Error: Cannot use cartesian coordinates in 2D.\n";
        return 1;
    }

    // Read data from the surface file
    surface.read_data();
    std::cout << "Surface data read successfully." << std::endl;
    // print data
    std::cout << "Number of points: " << surface.npoints << std::endl;

    // Create a ParticleSystem object
    ParticleSystem particle_system;
    // Read the particle list from a file
    std::string tables_path = config.get_string("tables_path");
    //particle_system.read_particle_list("../particles/pdg.dat");
    particle_system.read_particle_list(tables_path + "/pdg.dat");
    std::cout << "Particle list read successfully." << std::endl;
    // print data
    std::cout << "Number of particles: " << particle_system.nparticles << std::endl;


    // Create a NumericalIntegrator object
    NumericalIntegrator integrator;
    // Load integration tables for Gauss-Legendre and Gauss-Laguerre
    //integrator.load_integration_tables("../integration_tables/gauss_legendre.dat", "legendre");
    //integrator.load_integration_tables("../integration_tables/gauss_laguerre_a0.dat", "laguerre");
    //integrator.load_integration_tables("../integration_tables/gauss_laguerre_a1.dat", "laguerre_alpha1");
    integrator.load_integration_tables(tables_path + "/gauss_legendre.dat", "legendre");
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a0.dat", "laguerre");
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a1.dat", "laguerre_alpha1");


    //test_Gauss_Legendre(integrator);
    //test_Gauss_Laguerre(integrator);
    //test_trapezoidal(integrator);
    //test_Gauss_Laguerre_alpha1(integrator);

    // Create a Sampler object
    Sampler sampler(config);
    // Sample particles from the particle system
    if (config.get_string("sampling_method") == "conserved_charge") {
        std::cout << "Using conserved charge sampling method." << std::endl;
        sampler.conserved_charge_sampling(particle_system, surface, integrator);
    } else {
        std::cout << "Using regular sampling method." << std::endl;
        sampler.sample_unconstrained(particle_system, surface, integrator);
    }


    return 0;
}