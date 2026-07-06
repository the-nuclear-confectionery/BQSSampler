#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include "settings.h"
#include "globals.h"
#include "fluid_cell.h"

class Surface {
public:
    // Constructor to initialize the path
    Surface(const std::string& path, const Settings& settings);

    //store settings reference
    const Settings& settings;


    // Method to read data from the file
    void read_data();

    // Copy all per-cell fields of point icell into a self-contained FluidCell
    FluidCell get_cell(int icell) const;

    // Class variables
    std::string path;       // Path to the surface data file
    int npoints;             // Number of surface points (calculated automatically)
    double Btotal, Stotal, Qtotal; // Total charge in the initial condition

    std::vector<double> tau, x, y, eta; // Spacetime coordinates
    std::vector<double> ut, ux, uy, ueta; // Fluid velocity components
    std::vector<double> dsigma_t, dsigma_x, dsigma_y, dsigma_eta; // Surface element components
    std::vector<double> E, T, P, s; // Thermodynamic quantities
    std::vector<double> muB, muS, muQ; // Chemical potentials
    std::vector<double> u_dot_dsigma; // u dot dsigma
    std::vector<double> bulk; // Bulk pressure
    std::vector<double> shv_tt, shv_tx, shv_ty, shv_teta, 
                        shv_xx, shv_xy, shv_xeta, shv_yy, shv_yeta, shv_etaeta; // Shear stress tensor components

    //densities
    std::vector<double> rhoB, rhoS, rhoQ; // Densities of baryon, strangeness, and charge
    // diffusion components
    std::vector<double> diff_B0, diff_Bx, diff_By, diff_Beta;
    std::vector<double> diff_S0, diff_Sx, diff_Sy, diff_Seta;
    std::vector<double> diff_Q0, diff_Qx, diff_Qy, diff_Qeta;

    // EoS type per surface point (1 = table, 0 = HRG); populated only when eos_column=true
    std::vector<int> eos_type;

    // Aux quantities
    std::vector<double> N_baryons_cell;
    std::vector<double> N_antibaryons_cell;
    std::vector<double> N_strange_mesons_sminus_cell;
    std::vector<double> N_strange_mesons_splus_cell;
    std::vector<double> N_charged_mesons_qplus_cell;
    std::vector<double> N_charged_mesons_qminus_cell;
    std::vector<double> N_neutral_mesons_cell;
    std::vector<double> N_all_particles_cell;



    double T_average = 0.0, E_average = 0.0, P_average = 0.0, nB_average = 0.0, total_surface_volume = 0.0;
    double  muB_average = 0.0, muQ_average = 0.0, muS_average = 0.0;

    
};

#endif // SURFACE_H 