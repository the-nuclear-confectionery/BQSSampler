
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



inline double get_equilibrium_density( double pLRF_bar, const ThermalParams& params) {
    double Ebar = sqrt(pLRF_bar * pLRF_bar + params.mbar * params.mbar);
    double chem = params.baryon  * params.alphaB +
                  params.charge  * params.alphaQ +
                  params.strange * params.alphaS;
    
    if (params.sign == -1 && chem >= Ebar) {
        return 0.0; 
    }
    
    return pLRF_bar * exp(pLRF_bar) / (exp(Ebar - chem) + params.sign);
}


#endif // INTEGRANDS_H