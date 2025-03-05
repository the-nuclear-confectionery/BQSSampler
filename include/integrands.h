
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