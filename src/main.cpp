#include "surface.h"
#include "particle.h"
#include "numerical_integrator.h"
#include "sampler.h"
#include "particle_system.h"
#include "globals.h"
#include <iostream>

/// @brief Test Gauss-Legendre integration using a known function
void test_Gauss_Legendre(const NumericalIntegrator& integrator) {
    // Function to integrate: f(x) = x^2 over [-1,1]
    auto f = [](double x) { return x * x; };
    
    // Analytical result: integral of x^2 from -1 to 1 is 2/3
    double result = integrator.Gauss_Legendre(f, -1.0, 1.0);
    double expected = 2.0 / 3.0;
    
    std::cout << "Gauss-Legendre Test: "
              << "Computed = " << result
              << ", Expected = " << expected
              << ", Error = " << fabs(result - expected) << std::endl;
}

/// @brief Test Gauss-Laguerre integration using a known function
void test_Gauss_Laguerre(const NumericalIntegrator& integrator) {
    // Function to integrate: f(x) = e^(-x) * x^2 over [0, infinity]
    auto f = [](double x, const ThermalParams&) { return x * x; };
    
    // Analytical result: integral of x^2 * e^(-x) from 0 to infinity is 2
    std::vector<double> roots = integrator.roots_gauss_laguerre[0];
    std::vector<double> weights = integrator.weights_gauss_laguerre[0];
    ThermalParams params{}; // Unused but required by function signature
    
    double result = integrator.gauss_quadrature(f, params, roots, weights);
    double expected = 2.0;
    
    std::cout << "Gauss-Laguerre Test: "
              << "Computed = " << result
              << ", Expected = " << expected
              << ", Error = " << fabs(result - expected) << std::endl;
}

/// @brief Test Trapezoidal integration using a known function
void test_trapezoidal(const NumericalIntegrator& integrator) {
    // Function to integrate: f(x) = x^2 over [0,1]
    auto f = [](double x) { return x * x; };
    
    // Analytical result: integral of x^2 from 0 to 1 is 1/3
    double result = integrator.trapezoidal(f, 0.0, 1.0, 1000);
    double expected = 1.0 / 3.0;
    
    std::cout << "Trapezoidal Test: "
              << "Computed = " << result
              << ", Expected = " << expected
              << ", Error = " << fabs(result - expected) << std::endl;
}





/// @brief Test Gauss-Laguerre integration with alpha=1 using a known function
void test_Gauss_Laguerre_alpha1(const NumericalIntegrator& integrator) {
    // Function to integrate: f(x) = e^(-x) * x^3 over [0, infinity]
    auto f = [](double x, const ThermalParams&) { return x * x ; };
    
    // Analytical result: integral of x^3 * e^(-x) from 0 to infinity is 6
    std::vector<double> roots = integrator.roots_gauss_laguerre[1];
    std::vector<double> weights = integrator.weights_gauss_laguerre[1];
    ThermalParams params{}; // Unused but required by function signature
    
    double result = integrator.gauss_quadrature(f, params, roots, weights);
    double expected = 6.0;
    
    std::cout << "Gauss-Laguerre Test (alpha=1): "
              << "Computed = " << result
              << ", Expected = " << expected
              << ", Error = " << fabs(result - expected) << std::endl;
}


int main() {

    // Create a Surface object
    Surface surface("surface.dat");

    // Read data from the surface file
    surface.read_data();
    std::cout << "Surface data read successfully." << std::endl;
    // print data
    std::cout << "Number of points: " << surface.npoints << std::endl;

    // Create a ParticleSystem object
    ParticleSystem particle_system;
    // Read the particle list from a file
    particle_system.read_particle_list("../particles/pdg.dat");
    std::cout << "Particle list read successfully." << std::endl;
    // print data
    std::cout << "Number of particles: " << particle_system.nparticles << std::endl;


    // Create a NumericalIntegrator object
    NumericalIntegrator integrator;
    // Load integration tables for Gauss-Legendre and Gauss-Laguerre
    integrator.load_integration_tables("../integration_tables/gauss_legendre.dat", "legendre");
    integrator.load_integration_tables("../integration_tables/gauss_laguerre_a0.dat", "laguerre");
    integrator.load_integration_tables("../integration_tables/gauss_laguerre_a1.dat", "laguerre_alpha1");

    //test_Gauss_Legendre(integrator);
    //test_Gauss_Laguerre(integrator);
    //test_trapezoidal(integrator);
    //test_Gauss_Laguerre_alpha1(integrator);

    // Create a Sampler object
    Sampler sampler;
    // Sample particles from the particle system
    sampler.sample(particle_system, surface, integrator);

    return 0;
}