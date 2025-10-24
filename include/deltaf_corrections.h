#ifndef DELTAF_CORRECTIONS_H
#define DELTAF_CORRECTIONS_H

#include <cmath>
#include "lrf.h"
#include "integrands.h"
#include "settings.h"

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
    double /*bulkPi*/
) {
    return 0.0;
}

// --- Diffusion correction (currently zero) ---
inline double df_diffusion(
    LRF& lrf,
    double q_B0, double q_Bx, double q_By, double q_Beta,
    double q_S0, double q_Sx, double q_Sy, double q_Seta,
    double q_Q0, double q_Qx, double q_Qy, double q_Qeta
) {
    return 0.0;
}

// ======================================================
// Main wrapper: sums all corrections
// ======================================================
inline double df_corrections(
    const Settings& settings,
    LRF& lrf,
    double tau_squared,
    const double pLRF[4],
    double T,
    double E,
    double P,
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
    if (settings.get_bool("delta_f_diffusion")) {
        double df_b = df_bulk(diss_params.bulk); 
        df_total += df_b;
    }
    if (settings.get_bool("delta_f_bulk")) {
        double df_d = df_diffusion(lrf, diss_params.q_B0, diss_params.q_Bx, diss_params.q_By, diss_params.q_Beta,
                                   diss_params.q_S0, diss_params.q_Sx, diss_params.q_Sy, diss_params.q_Seta,
                                   diss_params.q_Q0, diss_params.q_Qx, diss_params.q_Qy, diss_params.q_Qeta);
        df_total += df_d;
    }
    return df_total;
}

#endif // DELTAF_CORRECTIONS_H
