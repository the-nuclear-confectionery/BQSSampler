#ifndef NUMERICAL_INTEGRATOR_H
#define NUMERICAL_INTEGRATOR_H

#include <vector>
#include <string>
#include <functional>
#include "integrands.h"

class NumericalIntegrator {
public:
    // Constructor: Initialize the method
    NumericalIntegrator();

    // Load roots and weights for a specified quadrature method (Gauss-Legendre or Gauss-Laguerre)
    void load_integration_tables(const std::string& filename, const std::string& method, int alpha = 0);

    // Integration functions
    double gauss_quadrature(double integrand(double root, const ThermalParams& params), const ThermalParams& params,
                            const std::vector<double>& roots, const std::vector<double>& weights) const;
    double Gauss_Legendre(const std::function<double(double)>& f, double a, double b) const;
    double trapezoidal(const std::function<double(double)>& f, double a, double b, int n) const;


    std::vector<std::vector<double>> roots_gauss_laguerre;   // One table for each alpha
    std::vector<std::vector<double>> weights_gauss_laguerre; // One table for each alpha
    std::vector<double> roots_gauss_legendre; // Roots for Gauss-Legendre
    std::vector<double> weights_gauss_legendre; // Weights for Gauss-Legendre
};




#endif 