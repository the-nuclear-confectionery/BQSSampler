#ifndef LRF_H
#define LRF_H

#include <cmath>
#include <algorithm>
#include <iostream>

class LRF {
public:
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

    // Momentums
    double pLRF_t, pLRF_x, pLRF_y, pLRF_z;
    double pLab_tau, pLab_x, pLab_y, pLab_eta;

    LRF(double ut_in, double ux_in, double uy_in, double ueta_in,
        double dsigma_t_in, double dsigma_x_in, 
        double dsigma_y_in, double dsigma_n_in, 
        double tau) {

        ut = ut_in;
        ux = ux_in;
        uy = uy_in;
        ueta = ueta_in;

        uperp =  sqrt(ux * ux  +  uy * uy);
        utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

        xt = uperp;
        zt = 1.;
        xeta = 0;
        zeta = 1.;

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

    /**
     * Perform a Lorentz boost of a four-momentum using SMASH formulation
     * 
     * @param beta_x, beta_y, beta_z: velocity components (in units of c)
     * @param p_in: input four-momentum (E, px, py, pz)
     * @param p_out: boosted four-momentum (output)
     */
    void lorentz_boost(double beta_x, double beta_y, double beta_z, 
                       double p_in[4], double p_out[4]) {

        double beta_sq = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;
        if (beta_sq >= 1.0) {
            std::cerr << "Warning: Attempted superluminal boost!" << std::endl;
            return;
        }

        // Compute gamma factor
        double gamma = 1.0 / sqrt(1.0 - beta_sq);

        // Compute p ⋅ β
        double p_dot_beta = p_in[1] * beta_x + p_in[2] * beta_y + p_in[3] * beta_z;

        // Compute transformed energy (time component)
        double p0_prime = gamma * (p_in[0] - p_dot_beta);

        // Compute space components
        double factor = (gamma - 1.0) * p_dot_beta / beta_sq - gamma * p_in[0];

        p_out[0] = p0_prime;
        p_out[1] = p_in[1] + beta_x * factor;
        p_out[2] = p_in[2] + beta_y * factor;
        p_out[3] = p_in[3] + beta_z * factor;
    }

    /**
     * Boost particle momentum from LRF to the Lab frame using the new Lorentz boost method.
     */
    void boost_momentum_to_lab(double tau_squared, double pLRF[4]) {
        pLRF_t = pLRF[0];
        pLRF_x = pLRF[1];
        pLRF_y = pLRF[2];
        pLRF_z = pLRF[3];

        double beta_x = ux / ut;
        double beta_y = uy / ut;
        double beta_z = ueta / ut;

        double pLab[4];
        lorentz_boost(beta_x, beta_y, beta_z, pLRF, pLab);

        pLab_tau = pLab[0];
        pLab_x = pLab[1];
        pLab_y = pLab[2];
        pLab_eta = pLab[3];
    }

    /**
     * Boost surface element from Lab frame to LRF using Lorentz transformation.
     */
    void boost_dsigma_to_lrf(double tau_squared) {
        double dsigma_lab[4] = {dsigma_t, dsigma_x, dsigma_y, dsigma_n};
        double dsigma_lrf[4];

        double beta_x = ux / ut;
        double beta_y = uy / ut;
        double beta_z = ueta / ut;

        lorentz_boost(-beta_x, -beta_y, -beta_z, dsigma_lab, dsigma_lrf);

        dsigma_t_lrf = dsigma_lrf[0];
        dsigma_x_lrf = dsigma_lrf[1];
        dsigma_y_lrf = dsigma_lrf[2];
        dsigma_z_lrf = dsigma_lrf[3];
    }

    /**
     * Compute dsigma magnitude in LRF.
     */
    void compute_dsigma_magnitude() {
        dsigma_space = sqrt(dsigma_x_lrf * dsigma_x_lrf +
                            dsigma_y_lrf * dsigma_y_lrf +
                            dsigma_z_lrf * dsigma_z_lrf);
        dsigma_magnitude = fabs(dsigma_t_lrf) + dsigma_space;
    }
};

#endif
