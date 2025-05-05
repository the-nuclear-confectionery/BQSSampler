#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include "settings.h"
#include "globals.h"

// A single surface point's data
struct SurfacePoint {
    std::vector<double> position = std::vector<double>(4);  // tau, x, y, eta
    std::vector<double> normal = std::vector<double>(4);    // normal_t, normal_x, normal_y, normal_eta
    std::vector<double> u = std::vector<double>(4); // ut, ux, uy, un
    std::vector<double> thermodynamics = std::vector<double>(4); // E, T, P, s
    std::vector<double> chemical_potentials = std::vector<double>(3); // muB, muS, muQ
    std::vector<double> densities = std::vector<double>(3);  // nB, nS, nQ
    std::vector<double> stress_tensor = std::vector<double>(10); // shv_tt, shv_tx, shv_ty, shv_teta, shv_xx, shv_xy, shv_xeta, shv_yy, shv_yeta, shv_etaeta
    double bulk = 0.0; // Bulk pressure
    double u_dot_dsigma = 0.0; // u dot dsigma
};



class Surface {
public:
    // Constructor to initialize the path
    Surface(const std::string& path, const Settings& settings);

    //store settings reference
    const Settings& settings;


    // Method to read data from the file
    void read_data();

    // Class variables
    std::string path;       // Path to the surface data file
    int npoints;             // Number of surface points (calculated automatically)
    std::vector<SurfacePoint> points;  // Vector of surface points

    std::vector<double> tau, x, y, eta; // Spacetime coordinates
    std::vector<double> ut, ux, uy, ueta; // Fluid velocity components
    std::vector<double> dsigma_t, dsigma_x, dsigma_y, dsigma_eta; // Surface element components
    std::vector<double> E, T, P, s; // Thermodynamic quantities
    std::vector<double> muB, muS, muQ; // Chemical potentials
    std::vector<double> nB, nS, nQ; // Densities
    std::vector<double> u_dot_dsigma; // u dot dsigma
    std::vector<double> bulk; // Bulk pressure
    std::vector<double> shv_tt, shv_tx, shv_ty, shv_teta, 
                        shv_xx, shv_xy, shv_xeta, shv_yy, shv_yeta, shv_etaeta; // Shear stress tensor components

    double T_average = 0.0, E_average = 0.0, P_average = 0.0, nB_average = 0.0, total_surface_volume = 0.0;
    double  muB_average = 0.0, muQ_average = 0.0, muS_average = 0.0;

    
};

#endif // SURFACE_H 