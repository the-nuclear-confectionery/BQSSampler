#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "surface.h"

/// @brief Construct a Surface object with the given path
Surface::Surface(const std::string& path) : path(path), npoints(0) {}

/// @brief Read data from the surface file
/// @details This function reads the surface data from a file, processes it, 
/// and stores it in the SurfacePoint objects. It also calculates the 
/// average thermodynamic quantities and total surface volume.
void Surface::read_data()
{
    std::cout << "Reading in freezeout surface in CCAKE boost invariant format" << std::endl;

    // Open the file to count lines
    // Prepare the file path using stringstream
    std::ostringstream surfdat_stream;
    surfdat_stream << path;
    std::ifstream surface_file(surfdat_stream.str().c_str());  // Open the file

    if (!surface_file.is_open()) {
        throw std::runtime_error("Could not open surface file.");
    }

    // Count the number of surface points (lines)
    std::string line;
    while (std::getline(surface_file, line)) {
        npoints++;  // Increment line counter for each valid line in the file
    }
    surface_file.close();  // Close the file after counting lines
    //npoints = 2;
    // Reopen the file to read data
    surface_file.open(surfdat_stream.str().c_str());
    if (!surface_file.is_open()) {
        throw std::runtime_error("Could not open surface file.");
    }

    double aux_tau, aux_x, aux_y, aux_eta;
    double aux_ut, aux_ux, aux_uy, aux_ueta;
    double aux_dsigma_t, aux_dsigma_x, aux_dsigma_y, aux_dsigma_eta;
    double aux_E, aux_T, aux_P, aux_s;
    double aux_muB, aux_muS, aux_muQ;
    double aux_nB, aux_nS, aux_nQ;
    double aux_u_dot_dsigma;
    double aux_bulk;
    double aux_shv_tt, aux_shv_tx, aux_shv_ty, aux_shv_teta, 
           aux_shv_xx, aux_shv_xy, aux_shv_xeta, aux_shv_yy, aux_shv_yeta, aux_shv_etaeta;

    int etacount = 0;

    for (int i = 0; i < npoints; i++)
    {

        // Read the data in the correct order and calculate directly
        surface_file >> aux_dsigma_t;
        surface_file >> aux_dsigma_x;
        surface_file >> aux_dsigma_y;
        surface_file >> aux_dsigma_eta;
        //std::cout << "dsigma_t: " << aux_dsigma_t << " dsigma_x: " << aux_dsigma_x << " dsigma_y: " << aux_dsigma_y << " dsigma_eta: " << aux_dsigma_eta << std::endl;
        // Read u and calculate gamma, ux, uy, un
        surface_file >> aux_ut;  // ut
        double gamma = aux_ut;
        surface_file >> aux_ux;  // ux
        surface_file >> aux_uy;  // uy
        surface_file >> aux_ueta;  // ueta
        //std::cout << "ut: " << aux_ut << " ux: " << aux_ux << " uy: " << aux_uy << " ueta: " << aux_ueta << std::endl;

        // Normalize normal vector components
        double m_over_sigma;
        surface_file >> m_over_sigma; // m/sigma from SPH
        double u_dot_n = (gamma * aux_dsigma_t + aux_ux *aux_dsigma_x + aux_uy * aux_dsigma_y + aux_ueta * aux_dsigma_eta) ;
        //std::cout << "u_dot_n: " << u_dot_n << std::endl;
       // std::cout << "m_over_sigma: " << m_over_sigma << std::endl;
       
        aux_dsigma_t = aux_dsigma_t* m_over_sigma / u_dot_n;
        aux_dsigma_x = aux_dsigma_x* m_over_sigma / u_dot_n;
        aux_dsigma_y = aux_dsigma_y* m_over_sigma / u_dot_n;
        aux_dsigma_eta = aux_dsigma_eta* m_over_sigma / u_dot_n;


        surface_file >> aux_bulk;

        // Read stress tensor components and apply conversion
        surface_file >> aux_shv_tt; // shv_tt
        aux_shv_tt *= HBARC;
        surface_file >> aux_shv_xx;  // shv_xx
        aux_shv_xx *= HBARC;
        surface_file >> aux_shv_yy;  // shv_yy
        aux_shv_yy *= HBARC;
        surface_file >> aux_shv_etaeta;  // shv_etaeta
        aux_shv_etaeta *= HBARC;
        surface_file >> aux_shv_xy;  // shv_xy
        aux_shv_xy *= HBARC;

        double dummy;
        surface_file >> dummy;  // shv_xeta
        surface_file >> dummy;  // shv_yeta

        aux_shv_xeta = 0.0;  // shv_xeta
        aux_shv_yeta = 0.0;  // shv_yeta
        aux_shv_tx = (aux_ux * aux_shv_xx + aux_uy * aux_shv_xy) / gamma;  // shv_tx
        aux_shv_ty = (aux_ux * aux_shv_xy + aux_uy * aux_shv_yy) / gamma;  // shv_ty
        aux_shv_teta = 0.0;  // shv_teta

        // Read spacetime position components
        surface_file >> aux_tau;  // tau
        surface_file >> aux_x;  // x
        surface_file >> aux_y;  // y
        surface_file >> aux_eta;  // eta

        double eta_ = 0.5 * log((aux_tau + aux_eta) / (aux_tau - aux_eta));
        if(aux_eta*aux_eta < aux_tau*aux_tau){
            if (eta_ < 0.5 && eta_ > -0.5){
                etacount++;
            }
        }

        // Read thermodynamic quantities
        surface_file >> aux_s; // s
        surface_file >> aux_E; // E
        aux_E *= HBARC;
        surface_file >> aux_T; // T


        // Read chemical potentials
        surface_file >> aux_muB;  // muB
        aux_muB *= HBARC;


        surface_file >> aux_muS;  // muS
        surface_file >> aux_muQ;  // muQ
        aux_muS *= HBARC;
        aux_muQ *= HBARC;

        
        std::string eos_name;
        //surface_file >> eos_name;


        double enthalpy;
        surface_file >> enthalpy;
        aux_P = enthalpy * HBARC - aux_E ;
        double cs2;
        surface_file >> cs2;

            

        // Calculate averages and total surface volume

        
        double ut_enforced = sqrt(1.0 + aux_ux * aux_ux + aux_uy * aux_uy + aux_ueta * aux_ueta);  // enforce normalization
        double dat = aux_dsigma_t;
        double dax = aux_dsigma_x;
        double day = aux_dsigma_y;
        double dan = aux_dsigma_eta;

        double udsigma = aux_ut * dat + aux_ux * dax + aux_uy * day + aux_ueta * dan;
        double dsigma_dsigma = dat * dat - dax * dax - day * day - dan * dan ;
        double dsigma_magnitude = fabs(udsigma) + sqrt(fabs(udsigma * udsigma - dsigma_dsigma));
        aux_u_dot_dsigma = udsigma;
        total_surface_volume += dsigma_magnitude;

        E_average += (aux_E * dsigma_magnitude);
        T_average += (aux_T * dsigma_magnitude);
        P_average += (aux_P * dsigma_magnitude);
        muB_average += (aux_muB * dsigma_magnitude);
        muQ_average += (aux_muQ * dsigma_magnitude);  // muQ is not computed here
        muS_average += (aux_muS * dsigma_magnitude);  // muS is not computed here
        nB_average += (0.0 * dsigma_magnitude);  // nB is not computed here
        
        // Store the data in the vectors
        tau.push_back(aux_tau);
        x.push_back(aux_x);
        y.push_back(aux_y);
        eta.push_back(aux_eta);
        ut.push_back(aux_ut);
        ux.push_back(aux_ux);
        uy.push_back(aux_uy);
        ueta.push_back(aux_ueta);
        dsigma_t.push_back(aux_dsigma_t);
        dsigma_x.push_back(aux_dsigma_x);
        dsigma_y.push_back(aux_dsigma_y);
        dsigma_eta.push_back(aux_dsigma_eta);
        E.push_back(aux_E);
        T.push_back(aux_T);
        P.push_back(aux_P);
        muB.push_back(aux_muB);
        muS.push_back(aux_muS);
        muQ.push_back(aux_muQ);
        u_dot_dsigma.push_back(aux_u_dot_dsigma);
        bulk.push_back(aux_bulk);
        shv_tt.push_back(aux_shv_tt);
        shv_tx.push_back(aux_shv_tx);
        shv_ty.push_back(aux_shv_ty);
        shv_teta.push_back(aux_shv_teta);
        shv_xx.push_back(aux_shv_xx);
        shv_xy.push_back(aux_shv_xy);
        shv_xeta.push_back(aux_shv_xeta);
        shv_yy.push_back(aux_shv_yy);
        shv_yeta.push_back(aux_shv_yeta);
        shv_etaeta.push_back(aux_shv_etaeta);
        //print all data
        //std::cout << "tau: " << aux_tau << " x: " << aux_x << " y: " << aux_y << " eta: " << aux_eta << std::endl;
        //std::cout << "ut: " << aux_ut << " ux: " << aux_ux << " uy: " << aux_uy << " ueta: " << aux_ueta << std::endl;
        //std::cout << "dsigma_t: " << aux_dsigma_t << " dsigma_x: " << aux_dsigma_x << " dsigma_y: " << aux_dsigma_y << " dsigma_eta: " << aux_dsigma_eta << std::endl;
        //std::cout << "E: " << aux_E << " T: " << aux_T << " P: " << aux_P << " s: " << aux_s << std::endl;
        //std::cout << "muB: " << aux_muB << " muS: " << aux_muS << " muQ: " << aux_muQ << std::endl;
        //std::cout << "bulk: " << aux_bulk << std::endl;
        //std::cout << "shv_tt: " << aux_shv_tt << " shv_tx: " << aux_shv_tx << " shv_ty: " << aux_shv_ty << " shv_teta: " << aux_shv_teta << std::endl;
        //std::cout << "shv_xx: " << aux_shv_xx << " shv_xy: " << aux_shv_xy << " shv_xeta: " << aux_shv_xeta << std::endl;
        //std::cout << "shv_yy: " << aux_shv_yy << " shv_yeta: " << aux_shv_yeta << " shv_etaeta: " << aux_shv_etaeta << std::endl;
        //std::cout << "u_dot_dsigma: " << aux_u_dot_dsigma << std::endl;
        

    }

    surface_file.close();
    std::cout << "Midrapidity cells: " << etacount << "out of" << npoints << std::endl;
    // Write averaged thermodynamic quantities to file
    T_average /= total_surface_volume;
    E_average /= total_surface_volume;
    P_average /= total_surface_volume;
    muB_average /= total_surface_volume;
    muS_average /= total_surface_volume;
    muQ_average /= total_surface_volume;
    nB_average /= total_surface_volume;

    std::cout << "Finished reading surface" << std::endl;
    std::ofstream thermal_average("average_thermodynamic_quantities.dat", std::ios_base::out);
    thermal_average << std::setprecision(15) << T_average << "\n" << E_average << "\n" << P_average << "\n" << muB_average << "\n" << nB_average;
    thermal_average.close();
}