#ifndef DELTAF_CORRECTIONS_H
#define DELTAF_CORRECTIONS_H

#include <cmath>
#include "lrf.h"
#include "integrands.h"
#include "settings.h"
#include "tables.h"

// ======================================================
// Subfunctions: each one handles its own boosting
// ======================================================

// --- Shear correction ---
inline double df_shear(
    LRF& lrf,
    double tau_squared,
    const double pLRF[4],   // [E, px, py, pz]
    double T,
    double E,
    double P,
    double shv_tt, double shv_tx, double shv_ty, double shv_tn,
    double shv_xx, double shv_xy, double shv_xn,
    double shv_yy, double shv_yn,
    double shv_nn
) {
    // Load π^{μν} components in lab/Milne
    lrf.shv_tt = shv_tt; lrf.shv_tx = shv_tx; lrf.shv_ty = shv_ty; lrf.shv_tn = shv_tn;
    lrf.shv_xx = shv_xx; lrf.shv_xy = shv_xy; lrf.shv_xn = shv_xn;
    lrf.shv_yy = shv_yy; lrf.shv_yn = shv_yn;
    lrf.shv_nn = shv_nn;

    // Boost to LRF
    lrf.boost_pimunu_to_lrf(tau_squared);

    // Compute p_i p_j π^{ij} in LRF
    const double px = pLRF[1];
    const double py = pLRF[2];
    const double pz = pLRF[3];

    const double pimunu_pmu_pnu =
        px * px * lrf.shv_xx_lrf +
        py * py * lrf.shv_yy_lrf +
        pz * pz * lrf.shv_zz_lrf +
        2.0 * (px * py * lrf.shv_xy_lrf +
               px * pz * lrf.shv_xz_lrf +
               py * pz * lrf.shv_yz_lrf);

    return 0.5 * pimunu_pmu_pnu /(T * T * (E + P));
}

// --- Bulk correction (currently zero) ---
inline double df_bulk(
    Table4D& table,
    double T,
    double E,
    double mass,
    const ThermalParams& thermal_params,
    double Pi
) {

    //calcualte interpolated value from table
    double mub = thermal_params.alphaB * T;
    double muq = thermal_params.alphaQ * T;
    double mus = thermal_params.alphaS * T;
    //print mub,muq,mus, T
   // std::cout << "Calculating bulk delta-f coefficients for T=" << T << " GeV, muB=" << mub << " GeV, muQ=" << muq << " GeV, muS=" << mus << " GeV." << std::endl;

    //if abs(muB), abs(muQ), abs(muS) > 0.8, print warning
    if (std::abs(mub) >=  0.8 || std::abs(muq) >=  0.8 || std::abs(mus) >= 0.8) {
        std::cerr << "Warning: |muB|, |muQ|, or |muS| > 0.8 GeV. Diffusion delta-f coefficients may be inaccurate." << std::endl;
    }
    //check if t < 0.1 or t > 0.2 GeV
    if (T <= 0.1 || T >= 0.2) {
        std::cerr << "Warning: T out of bounds for diffusion delta-f coefficients table." << std::endl;
    }    
    double A_T = table.interpolate(mub,mus,muq,T, 0);
    double A_E = table.interpolate(mub,mus,muq,T, 1);
    double A_B_baryon = table.interpolate(mub,mus,muq,T, 2);
    double A_B_charge = table.interpolate(mub,mus,muq,T, 3);
    double A_B_strange = table.interpolate(mub,mus,muq,T, 4);

    double delta_f = A_T*Pi*mass*mass + A_E*Pi*E*E
                     + A_B_baryon*Pi*E*thermal_params.baryon
                     + A_B_charge*Pi*E*thermal_params.charge
                     + A_B_strange*Pi*E*thermal_params.strange;

    return delta_f;
}

// --- Diffusion correction 
// V^mu p_<mu> = - V^i p^j g_{ij} = V^i p^i in LRF
inline double df_diffusion(
    LRF& lrf,
    Table4D& deltaf_table,
    const double pLRF[4],
    double tau_squared,
    double E,
    double T,
    const ThermalParams& thermal_params,
    double q_B0, double q_Bx, double q_By, double q_Beta,
    double q_S0, double q_Sx, double q_Sy, double q_Seta,
    double q_Q0, double q_Qx, double q_Qy, double q_Qeta
) {

    double muB = thermal_params.alphaB * T;
    double muQ = thermal_params.alphaQ * T;
    double muS = thermal_params.alphaS * T;


    double Vb_x_lrf, Vb_y_lrf, Vb_z_lrf;
    double Vq_x_lrf, Vq_y_lrf, Vq_z_lrf;
    double Vs_x_lrf, Vs_y_lrf, Vs_z_lrf;
    lrf.boost_Vmu_to_lrf(
        tau_squared, q_B0, q_Bx, q_By, q_Beta,
        Vb_x_lrf, Vb_y_lrf, Vb_z_lrf
    );
    lrf.boost_Vmu_to_lrf(
        tau_squared, q_Q0, q_Qx, q_Qy, q_Qeta,
        Vq_x_lrf, Vq_y_lrf, Vq_z_lrf
    );
    lrf.boost_Vmu_to_lrf(
        tau_squared, q_S0, q_Sx, q_Sy, q_Seta,
        Vs_x_lrf, Vs_y_lrf, Vs_z_lrf
    );


    double A_V_bb = deltaf_table.interpolate(muB, muQ, muS, T, 5);
    double A_V_bq = deltaf_table.interpolate(muB, muQ, muS, T, 6);
    double A_V_bs = deltaf_table.interpolate(muB, muQ, muS, T, 7);
    double A_V_qq = deltaf_table.interpolate(muB, muQ, muS, T, 8);
    double A_V_qs = deltaf_table.interpolate(muB, muQ, muS, T, 9);
    double A_V_ss = deltaf_table.interpolate(muB, muQ, muS, T, 10);
    double A_q_b = deltaf_table.interpolate(muB, muQ, muS, T, 11);
    double A_q_q = deltaf_table.interpolate(muB, muQ, muS, T, 12);
    double A_q_s = deltaf_table.interpolate(muB, muQ, muS, T, 13);

    double C_V_bx = A_V_bb * Vb_x_lrf + A_V_bq * Vq_x_lrf + A_V_bs * Vs_x_lrf;
    double C_V_by = A_V_bb * Vb_y_lrf + A_V_bq * Vq_y_lrf + A_V_bs * Vs_y_lrf;
    double C_V_bz = A_V_bb * Vb_z_lrf + A_V_bq * Vq_z_lrf + A_V_bs * Vs_z_lrf;
    double C_V_qx = A_V_bq * Vb_x_lrf + A_V_qq * Vq_x_lrf + A_V_qs * Vs_x_lrf;
    double C_V_qy = A_V_bq * Vb_y_lrf + A_V_qq * Vq_y_lrf + A_V_qs * Vs_y_lrf;
    double C_V_qz = A_V_bq * Vb_z_lrf + A_V_qq * Vq_z_lrf + A_V_qs * Vs_z_lrf;  
    double C_V_sx = A_V_bs * Vb_x_lrf + A_V_qs * Vq_x_lrf + A_V_ss * Vs_x_lrf;
    double C_V_sy = A_V_bs * Vb_y_lrf + A_V_qs * Vq_y_lrf + A_V_ss * Vs_y_lrf;
    double C_V_sz = A_V_bs * Vb_z_lrf + A_V_qs * Vq_z_lrf + A_V_ss * Vs_z_lrf;

    double C_Q_x = A_q_b * Vb_x_lrf + A_q_q * Vq_x_lrf + A_q_s * Vs_x_lrf;
    double C_Q_y = A_q_b * Vb_y_lrf + A_q_q * Vq_y_lrf + A_q_s * Vs_y_lrf;
    double C_Q_z = A_q_b * Vb_z_lrf + A_q_q * Vq_z_lrf + A_q_s * Vs_z_lrf;

    // q_a C_V_a_mu p^mu
    double delta_f = thermal_params.baryon * (C_V_bx * pLRF[1] + C_V_by * pLRF[2] + C_V_bz * pLRF[3])
                     + thermal_params.charge * (C_V_qx * pLRF[1] + C_V_qy * pLRF[2] + C_V_qz * pLRF[3])
                     + thermal_params.strange * (C_V_sx * pLRF[1] + C_V_sy * pLRF[2] + C_V_sz * pLRF[3]);

    // Add E c_Q_mu p^mu term
    delta_f += E * (C_Q_x * pLRF[1] + C_Q_y * pLRF[2] + C_Q_z * pLRF[3]);

    return delta_f;
}

// ======================================================
// Main wrapper: sums all corrections
// ======================================================
inline double df_corrections(
    const Settings& settings,
    Table4D& coefficients_table,
    LRF& lrf,
    double tau_squared,
    const double pLRF[4],
    double T,
    double E,
    double P,
    double mass,
    const ThermalParams& thermal_params,
    const DissipativeParams& diss_params

) {





    //check for settings and add each df automatically
    double df_total = 0.;
    if (settings.get_bool("delta_f_shear")) {
        double df_s = df_shear(
            lrf, tau_squared, pLRF, T, E, P,
            diss_params.shv_tt, diss_params.shv_tx, diss_params.shv_ty, diss_params.shv_teta,
            diss_params.shv_xx, diss_params.shv_xy, diss_params.shv_xeta,
            diss_params.shv_yy, diss_params.shv_yeta, diss_params.shv_etaeta
        );
        df_total += df_s;
    }
    if (settings.get_bool("delta_f_bulk")) {
        double df_b = df_bulk(
            coefficients_table,
            T,
            E,
            mass,
            thermal_params,
            diss_params.bulk
        );
        df_total += df_b;
    }
    if (settings.get_bool("delta_f_diffusion")) {
        double df_d = df_diffusion(lrf, coefficients_table, pLRF, tau_squared, E, T,thermal_params, 
                                   diss_params.q_B0, diss_params.q_Bx, diss_params.q_By, diss_params.q_Beta,
                                   diss_params.q_S0, diss_params.q_Sx, diss_params.q_Sy, diss_params.q_Seta,
                                   diss_params.q_Q0, diss_params.q_Qx, diss_params.q_Qy, diss_params.q_Qeta);
        df_total += df_d;
    }
    return df_total;
}

#endif // DELTAF_CORRECTIONS_H
