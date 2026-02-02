#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "surface.h"


/// @brief Construct a Surface object with the given path
Surface::Surface(const std::string& path , const Settings& settings)
    : path(path), settings(settings), npoints(0) {}

/// @brief Read data from the surface file
/// @details This function reads the surface data from a file, processes it, 
/// and stores it in the SurfacePoint objects. It also calculates the 
/// average thermodynamic quantities and total surface volume.
void Surface::read_data()
{
    int D = settings.get_int("dimension");
    std::string mode = settings.get_string("mode");
    if (mode == "ccakev1" && D==3){
        throw std::runtime_error("ccakev1 mode is not supported in 3D");
    }

    std::string coordinate_system = settings.get_string("coordinate_system");
    if (D == 2) std::cout << "Reading surface data in 2D..." << std::endl;
    else if (D == 3) std::cout << "Reading surface data in 3D..." << std::endl;
    else throw std::runtime_error("D must be either 2 or 3.");



    
    //check coordinate system
    if (coordinate_system == "cartesian") std::cout << "Using cartesian coordinates..." << std::endl;
    else if (coordinate_system == "hyperbolic") std::cout << "Using hyperbolic coordinates..." << std::endl;
    else throw std::runtime_error("Coordinate system must be either cartesian or hyperbolic.");

    // Open the file to count lines
    std::ostringstream surfdat_stream;
    surfdat_stream << path;
    std::ifstream surface_file(surfdat_stream.str().c_str());
    if (!surface_file.is_open()) {
        throw std::runtime_error("Could not open surface file.");
    }

    // Count number of surface points (ignore header/comments)
    npoints = 0;
    std::string line;
    bool header_read = false;
    double Btmp, Stmp, Qtmp;
    while (std::getline(surface_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            // Try to parse first header "# B S Q"
            if (!header_read) {
                std::istringstream hs(line.substr(1));
                if (hs >> Btmp >> Stmp >> Qtmp) {
                    // Store if you have members for totals (optional):
                    // this->Btot = Btmp; this->Stot = Stmp; this->Qtot = Qtmp;
                    header_read = true;
                }
            }
            continue; // do not count comment lines
        }
        npoints++;
    }
    surface_file.close();
 
    //store the total charges
    Btotal = Btmp;
    Stotal = Stmp;
    Qtotal = Qtmp;
    // Reopen the file to read data
    surface_file.open(surfdat_stream.str().c_str());
    if (!surface_file.is_open()) {
        throw std::runtime_error("Could not open surface file.");
    }

    // Skip comment/header lines at the top (including "# B S Q")
    std::streampos data_start_pos;
    while (true) {
        data_start_pos = surface_file.tellg();
        if (!std::getline(surface_file, line)) break;
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        // first data line found: rewind to its beginning
        surface_file.clear();
        surface_file.seekg(data_start_pos);
        break;
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

        if ( D == 2) {
            aux_dsigma_eta = 0.0;
        } else {
            surface_file >> aux_dsigma_eta;  // dsigma_eta
        }
        // covariant dsigma components


        // Read u and calculate gamma, ux, uy, un
        surface_file >> aux_ut;  // ut
        double gamma = aux_ut;
        surface_file >> aux_ux;  // ux
        surface_file >> aux_uy;  // uy
        if ( D == 2) {
            aux_ueta = 0.0;  // ueta
        } else {
            surface_file >> aux_ueta;  // ueta
        }
        //contravariant four-velocity components

        // Normalize normal vector components
        double m_over_sigma;
        surface_file >> m_over_sigma; // m/sigma from SPH
        double u_dot_n = (gamma * aux_dsigma_t + aux_ux *aux_dsigma_x + aux_uy * aux_dsigma_y + aux_ueta * aux_dsigma_eta) ;

        aux_dsigma_t = aux_dsigma_t* m_over_sigma / u_dot_n;
        aux_dsigma_x = aux_dsigma_x* m_over_sigma / u_dot_n;
        aux_dsigma_y = aux_dsigma_y* m_over_sigma / u_dot_n;
        aux_dsigma_eta = aux_dsigma_eta* m_over_sigma / u_dot_n;


        surface_file >> aux_bulk;
        aux_bulk *= HBARC;

        // Read contravariant stress tensor components and apply conversion
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

        if(D == 2) {
            aux_shv_xeta = 0.0;  // shv_xeta
            aux_shv_yeta = 0.0;  // shv_yeta
        } else {
            surface_file >> aux_shv_xeta;  // shv_xeta
            surface_file >> aux_shv_yeta;  // shv_yeta
        }
        aux_shv_xeta *= HBARC;
        aux_shv_yeta *= HBARC;



        // Read contravariant spacetime position components
        surface_file >> aux_tau;  // tau
        surface_file >> aux_x;  // x
        surface_file >> aux_y;  // y

        if (D == 2) {
            aux_eta = 0.0;  // eta
        } else {
            surface_file >> aux_eta;  // eta
        }


        // Determine g33 based on coordinate system
        double g33;
        if (coordinate_system == "hyperbolic") {
            g33 = -aux_tau * aux_tau;  // g_ηη = -τ²
        } else {
            g33 = -1.0;                // g_zz = -1
        }


        // Lower the index of the four-velocity: u_mu = g_mu_nu * u^nu
        double ut_cov = aux_ut;           // g_ττ = 1 or g_tt = 1
        double ux_cov = -aux_ux;          // g_xx = -1
        double uy_cov = -aux_uy;          // g_yy = -1
        double un_cov = g33 * aux_ueta;     // g_ηη = -τ² or g_zz = -1

        // Compute temporal-spatial shear stress components π^{τi}
        aux_shv_tx = -(ux_cov * aux_shv_xx + uy_cov * aux_shv_xy + un_cov * aux_shv_xeta) / aux_ut;
        aux_shv_ty = -(ux_cov * aux_shv_xy + uy_cov * aux_shv_yy + un_cov * aux_shv_yeta) / aux_ut;
        aux_shv_teta = -(ux_cov * aux_shv_xeta + uy_cov * aux_shv_yeta + un_cov * aux_shv_etaeta) / aux_ut;


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
        if (mode == "ccakev1") {
            surface_file >> eos_name;  // eos_name
        };

        //surface_file >> eos_name;


        double enthalpy;
        surface_file >> enthalpy;
        aux_P = enthalpy * HBARC - aux_E ;
        double cs2;
        surface_file >> cs2;
        //print all the data read
        //std::cout << "Read surface point " << i+1 << "/" << npoints << ": "
        //          << "tau=" << aux_tau << ", x=" << aux_x << ", y=" << aux_y << ", eta=" << aux_eta
        //          << ", ut=" << aux_ut << ", ux=" << aux_ux << ", uy=" << aux_uy << ", ueta=" << aux_ueta
        //          << ", dsigma_t=" << aux_dsigma_t << ", dsigma_x=" << aux_dsigma_x << ", dsigma_y=" << aux_dsigma_y << ", dsigma_eta=" << aux_dsigma_eta
        //          << ", E=" << aux_E << ", T=" << aux_T << ", P=" << aux_P
        //          << ", muB=" << aux_muB << ", muS=" << aux_muS << ", muQ=" << aux_muQ
        //          << ", bulkPi=" << aux_bulk
        //          << ", shv_tt=" << aux_shv_tt << ", shv_tx=" << aux_shv_tx << ", shv_ty=" << aux_shv_ty << ", shv_teta=" << aux_shv_teta
        //          << ", shv_xx=" << aux_shv_xx << ", shv_xy=" << aux_shv_xy << ", shv_xeta=" << aux_shv_xeta
        //          << ", shv_yy=" <<  aux_shv_yy<< ", shv_yeta=" << aux_shv_yeta<< ", shv_etaeta=" << aux_shv_etaeta
        //          << std::endl;
        //for future
        double aux_nB = 0.0;
        double aux_nS = 0.0;
        double aux_nQ = 0.0;
        surface_file >> aux_nB;  // nB
        surface_file >> aux_nS;  // nS
        surface_file >> aux_nQ;  // nQ

        //diffusion components 
        double aux_diff_Bx = 0.0;
        double aux_diff_By = 0.0;
        double aux_diff_Beta = 0.0;
        double aux_diff_Sx = 0.0;
        double aux_diff_Sy = 0.0;
        double aux_diff_Seta = 0.0;
        double aux_diff_Qx = 0.0;
        double aux_diff_Qy = 0.0;
        double aux_diff_Qeta = 0.0;

        double dummy_diff_0; //will be discarded
        surface_file >> dummy_diff_0; //diff_B0
        surface_file >> aux_diff_Bx;
        surface_file >> aux_diff_By;
        surface_file >> aux_diff_Beta;
            
        surface_file >> dummy_diff_0; //diff_S0
        surface_file >> aux_diff_Sx;
        surface_file >> aux_diff_Sy;
        surface_file >> aux_diff_Seta;
        
        surface_file >> dummy_diff_0; //diff_Q0
        surface_file >> aux_diff_Qx;
        surface_file >> aux_diff_Qy;
        surface_file >> aux_diff_Qeta;
        //calculate q0 components
        double aux_diff_B0 = -(aux_diff_Bx * ux_cov + aux_diff_By * uy_cov + aux_diff_Beta * un_cov) / aux_ut;
        double aux_diff_S0 = -(aux_diff_Sx * ux_cov + aux_diff_Sy * uy_cov + aux_diff_Seta * un_cov) / aux_ut;
        double aux_diff_Q0 = -(aux_diff_Qx * ux_cov + aux_diff_Qy * uy_cov + aux_diff_Qeta * un_cov) / aux_ut;




        // Calculate averages and total surface volume

        
        double ut_enforced = sqrt(1.0 + aux_ux * aux_ux + aux_uy * aux_uy + aux_ueta * aux_ueta);  // enforce normalization
        double dat = aux_dsigma_t;
        double dax = aux_dsigma_x;
        double day = aux_dsigma_y;
        double dan = aux_dsigma_eta;

        double udsigma = aux_ut * dat + aux_ux * dax + aux_uy * day + aux_ueta * dan;

        aux_u_dot_dsigma = udsigma;

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

        if (settings.get_bool("use_mub")) {
            muB.push_back(aux_muB);
        } else {
            muB.push_back(0.0);
        }
        if (settings.get_bool("use_mus")) {
            muS.push_back(aux_muS);
        } else {
            muS.push_back(0.0);
        }
        if (settings.get_bool("use_muq")) {
            muQ.push_back(aux_muQ);
        } else {
            muQ.push_back(0.0);
        }
        rhoB.push_back(aux_nB);
        rhoS.push_back(aux_nS);
        rhoQ.push_back(aux_nQ);
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
        diff_B0.push_back(aux_diff_B0);
        diff_Bx.push_back(aux_diff_Bx);
        diff_By.push_back(aux_diff_By);
        diff_Beta.push_back(aux_diff_Beta);
        diff_S0.push_back(aux_diff_S0);
        diff_Sx.push_back(aux_diff_Sx);
        diff_Sy.push_back(aux_diff_Sy);
        diff_Seta.push_back(aux_diff_Seta);
        diff_Q0.push_back(aux_diff_Q0);
        diff_Qx.push_back(aux_diff_Qx);
        diff_Qy.push_back(aux_diff_Qy);
        diff_Qeta.push_back(aux_diff_Qeta);

        //aux quantities
        N_baryons_cell.push_back(0.0);
        N_antibaryons_cell.push_back(0.0);
        N_strange_mesons_sminus_cell.push_back(0.0);
        N_strange_mesons_splus_cell.push_back(0.0);
        N_charged_mesons_qplus_cell.push_back(0.0);
        N_charged_mesons_qminus_cell.push_back(0.0);
        N_neutral_mesons_cell.push_back(0.0);
        N_all_particles_cell.push_back(0.0);

        

    }

    surface_file.close();


    std::cout << "Finished reading surface" << std::endl;

}