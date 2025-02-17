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

    LRF(double ut_in, double ux_in, double uy_in, double ueta_in,  double tau) {
        ut = ut_in;
        ux = ux_in;
        uy = uy_in;
        ueta = ueta_in;

        uperp =  sqrt(ux * ux  +  uy * uy);
        utperp = sqrt(1.0  +  ux * ux  +  uy * uy);

        double sinh_eta = tau * ueta / utperp;
        double cosh_eta = ut / utperp;

        xt = uperp * cosh_eta;
        zt = sinh_eta;
        xeta = uperp * sinh_eta / tau;
        zeta = cosh_eta / tau;

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
    }


    void test_orthonormality(double tau_squared) {
        // Test orthonormality of basis vectors
        double u_normal = fabs(ut * ut - ux * ux - uy * uy - tau_squared * ueta * ueta - 1.0);
        double x_normal = fabs(xt * xt - xx * xx - xy * xy - tau_squared * xeta * xeta + 1.0);
        double y_normal = fabs(-yx * yx - yy * yy + 1.0);
        double z_normal = fabs(zt * zt - tau_squared * zeta * zeta + 1.0);

        double u_dot_x = fabs(xt * ut - xx * ux - xy * uy - tau_squared * xeta * ueta);
        double u_dot_y = fabs(-yx * ux - yy * uy);
        double u_dot_z = fabs(zt * ut - tau_squared * zeta * ueta);

        double x_dot_y = fabs(-xx * yx - xy * yy);
        double x_dot_z = fabs(xt * zt - tau_squared * xeta * zeta);

        double u_orthogonal = std::max(u_dot_x, std::max(u_dot_y, u_dot_z));
        double x_orthogonal = std::max(x_dot_y, x_dot_z);

        double epsilon = 1.e-14;
        // Additional checks or logging can be added if needed
    }

    void boost_shv_to_lrf(double tau_squared) {
        // shv_ij_LRF = Xi.shv.Xj
        shv_xx_lrf = shv_tt * xt * xt + shv_xx * xx * xx + shv_yy * xy * xy + tau_squared * tau_squared * shv_nn * xeta * xeta
                     + 2.0 * (-xt * (shv_tx * xx + shv_ty * xy) + shv_xy * xx * xy + tau_squared * xeta * (shv_xn * xx + shv_yn * xy - shv_tn * xt));

        shv_xy_lrf = yx * (-shv_tx * xt + shv_xx * xx + shv_xy * xy + tau_squared * shv_xn * xeta)
                     + yy * (-shv_ty * xt + shv_xy * xx + shv_yy * xy + tau_squared * shv_yn * xeta);

        shv_xz_lrf = zt * (shv_tt * xt - shv_tx * xx - shv_ty * xy - tau_squared * shv_tn * xeta)
                     - tau_squared * zeta * (shv_tn * xt - shv_xn * xx - shv_yn * xy - tau_squared * shv_nn * xeta);

        shv_yy_lrf = shv_xx * yx * yx + 2.0 * shv_xy * yx * yy + shv_yy * yy * yy;
        shv_yz_lrf = -zt * (shv_tx * yx + shv_ty * yy) + tau_squared * zeta * (shv_xn * yx + shv_yn * yy);
        shv_zz_lrf = -(shv_xx_lrf + shv_yy_lrf);
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
};

#endif
