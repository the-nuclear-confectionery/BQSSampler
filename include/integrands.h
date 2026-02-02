
#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#include <cmath>

//thermal params
struct ThermalParams {
    double T;
    double alphaB;
    double alphaQ;
    double alphaS;
    double mbar;
    double baryon;
    double strange;
    double charge;
    double sign;
    double spin_degeneracy;
};

struct DissipativeParams {
    double shv_tt;
    double shv_tx;
    double shv_ty;
    double shv_teta;
    double shv_xx;
    double shv_xy;
    double shv_xeta;
    double shv_yy;
    double shv_yeta;
    double shv_etaeta;
    double bulk;
    double q_B0;
    double q_Bx;
    double q_By;
    double q_Beta;
    double q_S0;
    double q_Sx;
    double q_Sy;
    double q_Seta;
    double q_Q0;
    double q_Qx;
    double q_Qy;
    double q_Qeta;
};



inline constexpr double get_equilibrium_density(const double pLRF_bar, const ThermalParams& params) {
    const double Ebar = std::sqrt(pLRF_bar * pLRF_bar + params.mbar * params.mbar);
    const double chem = params.baryon * params.alphaB +
                        params.charge * params.alphaQ +
                        params.strange * params.alphaS;

    if (params.sign == -1 && chem >= Ebar) {
        return 0.0;
    }

    return pLRF_bar * std::exp(pLRF_bar) / (std::exp(Ebar - chem) + params.sign);
}



#endif // INTEGRANDS_H