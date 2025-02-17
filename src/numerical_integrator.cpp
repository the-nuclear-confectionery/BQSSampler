#include "numerical_integrator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>

/// @brief Default constructor for NumericalIntegrator
NumericalIntegrator::NumericalIntegrator() {}

/// @brief Load integration tables from a file
/// @details This function reads integration tables from a specified file 
/// and stores the roots and weights for different numerical integration methods.
/// @param filename The path to the file containing the integration tables
/// @param method The numerical integration method to load the tables for
void NumericalIntegrator::load_integration_tables(const std::string& filename, const std::string& method) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    int alpha = 0;  // Default to alpha = 0
    while (std::getline(file, line)) {
        double root, weight;
        std::istringstream ss(line);
        ss >> root >> weight;

        if (method == "legendre") {
            // Read Gauss-Legendre roots and weights
            roots_gauss_legendre.push_back(root);
            weights_gauss_legendre.push_back(weight);
        } else if (method == "laguerre") {
            // Read Gauss-Laguerre roots and weights for alpha = 0
            if (roots_gauss_laguerre.size() <= alpha) {
                roots_gauss_laguerre.resize(alpha + 1); // Resize if new alpha is encountered
                weights_gauss_laguerre.resize(alpha + 1);
            }
            roots_gauss_laguerre[alpha].push_back(root);
            weights_gauss_laguerre[alpha].push_back(weight);
        } else if (method == "laguerre_alpha1") {
            // Read Gauss-Laguerre roots and weights for alpha = 1
            alpha = 1;
            if (roots_gauss_laguerre.size() <= alpha) {
                roots_gauss_laguerre.resize(alpha + 1); // Resize if new alpha is encountered
                weights_gauss_laguerre.resize(alpha + 1);
            }
            roots_gauss_laguerre[alpha].push_back(root);
            weights_gauss_laguerre[alpha].push_back(weight);
        } else {
            std::cerr << "Error: Unknown method '" << method << "'!" << std::endl;
            file.close();
            return;
        }
    }

    file.close();
}


/// @brief Perform Gauss-Legendre integration
/// @details This function performs Gauss-Legendre integration using the
/// roots and weights stored for the method. The integral is computed as
/// int_{-1}^{1} f(x) dx, and the function is rescaled to the interval [a, b].
/// @param f The function to integrate, can be a lambda or std::function
/// @param a The lower limit of integration
/// @param b The upper limit of integration
double NumericalIntegrator::Gauss_Legendre(const std::function<double(double)>& f, double a, double b) const{
    double result = 0.0;

    // Rescale the points from [-1, 1] to [a, b] for Gauss-Legendre
    for (size_t i = 0; i < roots_gauss_legendre.size(); ++i) {
        double x = 0.5 * (b - a) * roots_gauss_legendre[i] + 0.5 * (a + b); // Rescaled point
        result += weights_gauss_legendre[i] * f(x);
    }

    return 0.5 * (b - a) * result;  // Final multiplication with the rescaling factor
}

/// @brief Perform trapezoidal integration
/// @details This function performs trapezoidal integration over the interval [a, b]
/// using n subintervals. 
/// @param f The function to integrate, can be a lambda or std::function
/// @param a The lower limit of integration
/// @param b The upper limit of integration
/// @param n The number of subintervals
double NumericalIntegrator::trapezoidal(const std::function<double(double)>& f, double a, double b, int n) const{
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }

    return h * sum;
}



double NumericalIntegrator::gauss_quadrature(double integrand(double root, const ThermalParams& params), const ThermalParams& params,
                            const std::vector<double>& roots, const std::vector<double>& weights) const {
    double integral = 0.0;
    for (size_t i = 0; i < roots.size(); ++i) {
        integral += weights[i] * integrand(roots[i], params);
    }
    return integral;
}

