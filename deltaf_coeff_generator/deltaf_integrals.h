#include <cmath>


static inline double odd_double_factorial(int n) {
    if (n <= 1) return 1.0;
    double r = 1.0;
    for (int k = n; k >= 3; k -= 2) r *= double(k);
    return r;
}

// Use with Gauss-Laguerre table alpha = 2 + 2*q
// Returns F_{kq}(pbar) such that
// ∫ dpbar pbar^alpha e^{-pbar} J_{kq}(pbar)
// equals ∫ dpbar pbar^(2q+2) Ebar^(k-2q-1) ffbar.
template<int k, int q>
inline double get_Jkq_integrand(double pbar, const ThermalParams& p) {
    const double mbar = p.mbar;
    const double Ebar = std::sqrt(pbar*pbar + mbar*mbar);

    const double chem =
        p.baryon  * p.alphaB +
        p.charge  * p.alphaQ +
        p.strange * p.alphaS;

    // boson 
    if (p.sign == -1.0 && chem >= Ebar) return 0.0;

    const double z = std::exp(Ebar - chem);
    const double denom = z + p.sign;              
    if (denom == 0.0) return 0.0;

    // ffbar = f0(1-Theta f0) = z/(z+Theta)^2
    const double ffbar = z / (denom * denom);

    return std::exp(pbar) * std::pow(Ebar, k - 2*q - 1) * ffbar;
}

template<int k, int q>
static double J_kq(const NumericalIntegrator& integ, const ThermalParams& p) {
    constexpr int alpha = 2 + 2*q;
    const auto& roots = integ.roots_gauss_laguerre.at(alpha);
    const auto& wts   = integ.weights_gauss_laguerre.at(alpha);

    const double integral = integ.gauss_quadrature(get_Jkq_integrand<k,q>, p, roots, wts);

    double prefactor = p.spin_degeneracy / (2.0 * M_PI * M_PI * std::pow(HBARC,3.0));
    return prefactor * odd_double_factorial(k + 2) * std::pow(p.T, k + 2) * integral;
}

