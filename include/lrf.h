#ifndef LRF_H
#define LRF_H

#include <cmath>
#include <algorithm>
#include <iostream>

class LRF {
public:

    std::string coordinate_system;
    // Milne basis vectors: u^mu, x^mu, y^mu, z^mu (nonzero components)
    double ut, ux, uy, ueta;
    double xt, xx, xy, xeta;
    double yx, yy;
    double zt, zeta;

    // Shear stress tensor components
    double shv_tt, shv_tx, shv_ty, shv_tn;
    double shv_xx, shv_xy, shv_xn;
    double shv_yy, shv_yn;
    double shv_nn;

    // Shear stress tensor components in LRF
    double shv_xx_lrf, shv_xy_lrf, shv_xz_lrf;
    double shv_yy_lrf, shv_yz_lrf;
    double shv_zz_lrf;

    // Surface element vector components
    double dsigma_t, dsigma_x, dsigma_y, dsigma_n;
    double dsigma_t_lrf, dsigma_x_lrf, dsigma_y_lrf, dsigma_z_lrf;
    double dsigma_magnitude, dsigma_space;


    // Uperp and utperp
    double uperp, utperp;

    // momentums
    double pLRF_t, pLRF_x, pLRF_y, pLRF_z;
    // momentum in lab frame
    double pLab_tau, pLab_x, pLab_y, pLab_eta;

    LRF(std::string  coordinate_system,
        double ut_in, double ux_in, double uy_in, double ueta_in,
        double dsigma_t_in, double dsigma_x_in, 
        double dsigma_y_in, double dsigma_n_in, 
        double tau) {

        this->coordinate_system = coordinate_system;

        ut = ut_in;
        ux = ux_in;
        uy = uy_in;
        ueta = ueta_in;

        uperp =  sqrt(ux * ux  +  uy * uy);
        utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

        double sinh_eta = tau * ueta / utperp;
        double cosh_eta = ut / utperp;

        xt = 0;
        zt = 0.0;
        xeta =0.0;
        zeta = 1.0;
        if (coordinate_system == "hyperbolic") {
            xt = uperp * cosh_eta;
            zt = sinh_eta;
            xeta = uperp * sinh_eta / tau;
            zeta = cosh_eta / tau;
        }

        xx = 1.0;
        yx = 0.0;
        xy = 0.0;
        yy = 1.0;

        if (uperp > 1.e-5) {
            xx = utperp * ux / uperp;
            yx = -uy / uperp;
            xy = utperp * uy / uperp;
            yy = ux / uperp;
        }

        dsigma_t = dsigma_t_in;
        dsigma_x = dsigma_x_in;
        dsigma_y = dsigma_y_in;
        dsigma_n = dsigma_n_in;
    }



    void boost_dsigma_to_lrf(double tau_squared) {
        dsigma_t_lrf = dsigma_t * ut + dsigma_x * ux + dsigma_y * uy + dsigma_n * ueta;
        dsigma_x_lrf = -(dsigma_t * xt + dsigma_x * xx + dsigma_y * xy + dsigma_n * xeta);
        dsigma_y_lrf = -(dsigma_x * yx + dsigma_y * yy);
        dsigma_z_lrf = -(dsigma_t * zt + dsigma_n * zeta);
    }

    void compute_dsigma_magnitude() {
        dsigma_space = sqrt(dsigma_x_lrf * dsigma_x_lrf + dsigma_y_lrf * dsigma_y_lrf + dsigma_z_lrf * dsigma_z_lrf);
        dsigma_magnitude = fabs(dsigma_t_lrf) + dsigma_space;
    }

    void boost_momentum_to_lab(double tau_squared, double pLRF[4]) {
        pLRF_t = pLRF[0];
        pLRF_x = pLRF[1];
        pLRF_y = pLRF[2];
        pLRF_z = pLRF[3];

        pLab_tau = pLRF_t * ut + pLRF_x * xt + pLRF_z * zt;
        pLab_x   = pLRF_t * ux + pLRF_x * xx + pLRF_y * yx;
        pLab_y   = pLRF_t * uy + pLRF_x * xy + pLRF_y * yy;
        pLab_eta = pLRF_t * ueta + pLRF_x * xeta + pLRF_z * zeta; 
    }

    void boost_pimunu_to_lrf(double tau_squared) {
        double g33 = -1.0;
        if (coordinate_system == "hyperbolic") {
            g33 = -tau_squared;
        }

        shv_xx_lrf = shv_tt * xt * xt
                   + shv_xx * xx * xx
                   + shv_yy * xy * xy
                   + g33 * g33 * shv_nn * xeta * xeta
                   + 2.0 * ( -xt * (shv_tx * xx + shv_ty * xy)
                             + shv_xy * xx * xy
                             - g33 * xeta * (shv_xn * xx + shv_yn * xy - shv_tn * xt) );

        shv_xy_lrf = yx * ( -shv_tx * xt + shv_xx * xx + shv_xy * xy - g33 * shv_xn * xeta )
                   + yy * ( -shv_ty * xt + shv_xy * xx + shv_yy * xy - g33 * shv_yn * xeta );

        shv_xz_lrf = zt * ( shv_tt * xt - shv_tx * xx - shv_ty * xy + g33 * shv_tn * xeta )
                   + g33 * zeta * ( shv_tn * xt - shv_xn * xx - shv_yn * xy + g33 * shv_nn * xeta );

        shv_yy_lrf = shv_xx * yx * yx + 2.0 * shv_xy * yx * yy + shv_yy * yy * yy;

        shv_yz_lrf = -zt * ( shv_tx * yx + shv_ty * yy )
                     -g33 * zeta * ( shv_xn * yx + shv_yn * yy );

        shv_zz_lrf = - ( shv_xx_lrf + shv_yy_lrf );
    }


    void boost_Vmu_to_lrf(double tau_squared, double Vt, double Vx, double Vy, double Vn,
                          double& Vx_LRF, double& Vy_LRF, double& Vz_LRF) {
        double g33 = -1.0;
        if (coordinate_system == "hyperbolic") {
            g33 = -tau_squared;
        }

        // V^i_LRF = - X^i_mu * V^mu
        Vx_LRF = - Vt * xt + Vx * xx + Vy * xy - g33 * Vn * xeta;
        Vy_LRF = Vx * yx + Vy * yy;
        Vz_LRF = - Vt * zt - g33 * Vn * zeta;
    }


};

#endif

