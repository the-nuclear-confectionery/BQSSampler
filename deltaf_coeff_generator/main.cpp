// deltaf_coeff_generator/main.cpp
#include "numerical_integrator.h"
#include "particle_system.h"
#include "settings.h"
#include "integrands.h"
#include "globals.h"
#include "deltaf_integrals.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

#include <omp.h>

// ----------------- grid helper -----------------
static std::vector<double> linspace(double a, double b, int n) {
    if (n <= 0) throw std::runtime_error("linspace: n must be > 0");
    std::vector<double> v(n);
    if (n == 1) { v[0] = a; return v; }
    const double d = (b - a) / double(n - 1);
    for (int i = 0; i < n; ++i) v[i] = a + d * double(i);
    return v;
}

// ----------------- outputs -----------------
struct RowOut_GRAD {
    double T, muB, muQ, muS;

    // bulk sector 
    double A_T, A_E, A_b, A_q, A_s;

    // diffusion sector 
    double CV_bb, CV_bq, CV_bs, CV_qq, CV_qs, CV_ss;
    double CQ_b, CQ_q, CQ_s;
};

struct RowOut_CE {
    double T, muB, muQ, muS;

    double J30, J32;
    double N20_b, N20_q, N20_s;
    double M10_bb, M10_bq, M10_bs, M10_qq, M10_qs, M10_ss;
    double M11_bb, M11_bq, M11_bs, M11_qq, M11_qs, M11_ss;
};

// ----------------- generic 3x3 inverse (clear + works for symmetric too) -----------------
static inline bool invert_3x3(const double A[3][3], double Ainv[3][3], double& det, double eps = 1e-30) {
    const double c00 =  A[1][1]*A[2][2] - A[1][2]*A[2][1];
    const double c01 = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
    const double c02 =  A[1][0]*A[2][1] - A[1][1]*A[2][0];

    const double c10 = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]);
    const double c11 =  A[0][0]*A[2][2] - A[0][2]*A[2][0];
    const double c12 = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]);

    const double c20 =  A[0][1]*A[1][2] - A[0][2]*A[1][1];
    const double c21 = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]);
    const double c22 =  A[0][0]*A[1][1] - A[0][1]*A[1][0];

    det = A[0][0]*c00 + A[0][1]*c01 + A[0][2]*c02;

    if (std::abs(det) < eps) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Ainv[i][j] = 0.0;
        return false;
    }

    const double invdet = 1.0 / det;

    Ainv[0][0] = c00 * invdet;  Ainv[0][1] = c10 * invdet;  Ainv[0][2] = c20 * invdet;
    Ainv[1][0] = c01 * invdet;  Ainv[1][1] = c11 * invdet;  Ainv[1][2] = c21 * invdet;
    Ainv[2][0] = c02 * invdet;  Ainv[2][1] = c12 * invdet;  Ainv[2][2] = c22 * invdet;

    return true;
}

static inline double scalar_contraction(const double x[3], const double A[3][3], const double y[3]) {
    double sum = 0.0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            sum += x[a] * A[a][b] * y[b];
    return sum;
}


int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " config.dat\n";
        return 1;
    }

    Settings config(argv[1]);
    const std::string tables_path = config.get_string("tables_path");

    ParticleSystem particle_system;
    particle_system.read_particle_list(tables_path + "/pdg.dat");
    std::cout << "Particle list read. nparticles=" << particle_system.nparticles << "\n";

    NumericalIntegrator integrator;
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a0.dat", "laguerre", 0);
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a1.dat", "laguerre", 1);
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a2.dat", "laguerre", 2);
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a3.dat", "laguerre", 3);
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a4.dat", "laguerre", 4);
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a5.dat", "laguerre", 5);
    integrator.load_integration_tables(tables_path + "/gauss_laguerre_a6.dat", "laguerre", 6);

    // ----------------- 4D grid declared in code -----------------
    const double T_min   = 0.10,  T_max   = 0.20;  const int T_pts   = 41;
    const double muB_min = -0.80, muB_max = 0.80;  const int muB_pts = 21;
    const double muQ_min = -0.80, muQ_max = 0.80;  const int muQ_pts = 21;
    const double muS_min = -0.80, muS_max = 0.80;  const int muS_pts = 21;

    const auto T_grid   = linspace(T_min,   T_max,   T_pts);
    const auto muB_grid = linspace(muB_min, muB_max, muB_pts);
    const auto muQ_grid = linspace(muQ_min, muQ_max, muQ_pts);
    const auto muS_grid = linspace(muS_min, muS_max, muS_pts);

    const long long nT  = (long long)T_grid.size();
    const long long nB  = (long long)muB_grid.size();
    const long long nQ  = (long long)muQ_grid.size();
    const long long nS  = (long long)muS_grid.size();
    const long long total = nT * nB * nQ * nS;

    std::cout << "Total grid points = " << total << "\n";
    std::cout << "OpenMP max threads = " << omp_get_max_threads() << "\n";

    std::vector<RowOut_GRAD> rows_grad((size_t)total);
    std::vector<RowOut_CE>   rows_ce((size_t)total);

#pragma omp parallel for schedule(dynamic)
    for (long long idx = 0; idx < total; ++idx) {

        // unravel idx -> (iB,iQ,iS,iT)
        const long long iT = idx % nT;
        const long long tmp1 = idx / nT;
        const long long iS = tmp1 % nS;
        const long long tmp2 = tmp1 / nS;
        const long long iQ = tmp2 % nQ;
        const long long iB = tmp2 / nQ;

        const double T   = T_grid[(size_t)iT];
        const double muB = muB_grid[(size_t)iB];
        const double muQ = muQ_grid[(size_t)iQ];
        const double muS = muS_grid[(size_t)iS];

        // Charge ordering: 0=B, 1=Q, 2=S
        // We'll store vectors/matrices in that order.
        double J40 = 0.0, J41 = 0.0;
        double A20 = 0.0, A22 = 0.0;

        double N30[3] = {0.0, 0.0, 0.0};
        double N31[3] = {0.0, 0.0, 0.0};
        double B10[3] = {0.0, 0.0, 0.0};

        double M20[3][3] = {{0.0}};
        double M21[3][3] = {{0.0}};

        // CE
        double J30 = 0.0, J32 = 0.0;
        double N20[3] = {0.0, 0.0, 0.0};
        double M10[3][3] = {{0.0}};
        double M11[3][3] = {{0.0}};

        for (int ip = 0; ip < particle_system.nparticles; ++ip) {
            const double m = particle_system.mass[ip];
            if (m == 0.0) continue;

            ThermalParams p;
            p.T      = T;
            p.mbar   = m / T;
            p.alphaB = muB / T;
            p.alphaQ = muQ / T;
            p.alphaS = muS / T;
            p.baryon  = particle_system.baryon[ip];
            p.charge  = particle_system.charge[ip];
            p.strange = particle_system.strange[ip];
            p.sign    = particle_system.theta[ip];
            p.spin_degeneracy = particle_system.spin_degeneracy[ip];

            const double charge_vec[3] = { p.baryon, p.baryon, p.baryon };
            const double m2 = m * m;

            const double j40 = J_kq<4,0>(integrator, p);
            const double j41 = J_kq<4,1>(integrator, p);
            const double j30 = J_kq<3,0>(integrator, p);
            const double j31 = J_kq<3,1>(integrator, p);
            const double j32 = J_kq<3,2>(integrator, p);
            const double j20 = J_kq<2,0>(integrator, p);
            const double j21 = J_kq<2,1>(integrator, p);
            const double j10 = J_kq<1,0>(integrator, p);
            const double j11 = J_kq<1,1>(integrator, p);

            J40 += j40;
            J41 += j41;

            J30 += j30;
            J32 += j32;

            A20 += m2 * j20;
            A22 += m2 * j21;

            for (int a = 0; a < 3; ++a) {
                N30[a] += charge_vec[a] * j30;
                N31[a] += charge_vec[a] * j31;
                B10[a] += charge_vec[a] * (m2 * j10);
            }

            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    M20[a][b] += charge_vec[a] * charge_vec[b] * j20;
                    M21[a][b] += charge_vec[a] * charge_vec[b] * j21;
                }
            }

            for (int a = 0; a < 3; ++a) N20[a] += charge_vec[a] * j20;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    M10[a][b] += charge_vec[a] * charge_vec[b] * j10;
                    M11[a][b] += charge_vec[a] * charge_vec[b] * j11;
                }
            }
        }

        // ---------------- BULK coefficients (divided by Pi) ----------------
        double bulk_A_T = 0.0, bulk_A_E = 0.0;
        double bulk_A_charge[3] = {0.0, 0.0, 0.0};

        double M20_inv[3][3] = {{0.0}};
        double detM20 = 0.0;
        const bool okM20 = invert_3x3(M20, M20_inv, detM20);

        if (okM20) {
            const double N30_Minv_N30 = scalar_contraction (N30, M20_inv, N30);
            const double N30_Minv_B10 = scalar_contraction (N30, M20_inv, B10);
            const double N31_Minv_N30 = scalar_contraction (N31, M20_inv, N30);
            const double N31_Minv_B10 = scalar_contraction (N31, M20_inv, B10);

            const double X = (A20 - N30_Minv_B10);
            const double Y = (J41 - N31_Minv_N30);
            const double U = (J40 - N30_Minv_N30);
            const double W = (A22 - N31_Minv_B10);

            const double Delta = X*Y - U*W;

            if (std::abs(Delta) > 0.0) {
                bulk_A_T = -(U / Delta);
                bulk_A_E = +(X / Delta);

                double rhs_vec[3];
                for (int b = 0; b < 3; ++b)
                    rhs_vec[b] = B10[b] * bulk_A_T + N30[b] * bulk_A_E;

                for (int a = 0; a < 3; ++a) {
                    double s = 0.0;
                    for (int b = 0; b < 3; ++b) s += M20_inv[a][b] * rhs_vec[b];
                    bulk_A_charge[a] = -s;
                }
            }
        }

        // ---------------- DIFFUSION coefficients (divided by V_b^mu) ----------------
        // Build K = M21 - (1/J41) N31 N31^T, invert it, then:
        //   CV_ab = -(K^{-1})^{ab}
        //   CQ_b  = (1/J41) * sum_a N31[a] * (K^{-1})^{ab}
        double CV[3][3] = {{0.0}};     // will store -(K^{-1})
        double CQ[3] = {0.0, 0.0, 0.0};

        if (J41 != 0.0) {
            const double invJ41 = 1.0 / J41;

            double K[3][3];
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    K[a][b] = M21[a][b] - invJ41 * (N31[a] * N31[b]);
                }
            }

            double K_inv[3][3] = {{0.0}};
            double detK = 0.0;
            const bool okK = invert_3x3(K, K_inv, detK);

            if (okK) {
                // CV = -K_inv
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        CV[a][b] = -K_inv[a][b];

                // CQ[b] = (1/J41) * sum_a N31[a] * K_inv[a][b]
                for (int b = 0; b < 3; ++b) {
                    double s = 0.0;
                    for (int a = 0; a < 3; ++a) s += N31[a] * K_inv[a][b];
                    CQ[b] = invJ41 * s;
                }
            }
        }

        rows_grad[(size_t)idx] = RowOut_GRAD{
            T, muB, muQ, muS,
            bulk_A_T, bulk_A_E, bulk_A_charge[0], bulk_A_charge[1], bulk_A_charge[2],
            CV[0][0], CV[0][1], CV[0][2], CV[1][1], CV[1][2], CV[2][2],
            CQ[0], CQ[1], CQ[2]
        };

        rows_ce[(size_t)idx] = RowOut_CE{
            J30, J32,
            N20[0], N20[1], N20[2],
            M10[0][0], M10[0][1], M10[0][2], M10[1][1], M10[1][2], M10[2][2],
            M11[0][0], M11[0][1], M11[0][2], M11[1][1], M11[1][2], M11[2][2]
        };
    }

    // ---------- write GRAD output ----------
    {
        std::ofstream out("grad_coeffs.dat");
        out << std::setprecision(16);
        out << "# "  << nB << " " << nQ << " " << nS << " "<< nT << " " << total << "\n";   


        for (long long idx = 0; idx < total; ++idx) {
            const auto& r = rows_grad[(size_t)idx];
            out << r.muB << " " << r.muQ << " " << r.muS << " " << r.T << " "
                << r.A_T << " " << r.A_E << " " << r.A_b << " " << r.A_q << " " << r.A_s << " "
                << r.CV_bb << " " << r.CV_bq << " " << r.CV_bs << " "
                << r.CV_qq << " " << r.CV_qs << " " << r.CV_ss << " "
                << r.CQ_b << " " << r.CQ_q << " " << r.CQ_s
                << "\n";
        }
    }

    // ---------- write CE output ----------
    {
        std::ofstream out("ce_moments.dat");
        out << std::setprecision(16);
        out << "# "  << nB << " " << nQ << " " << nS << " "<< nT << " " << total << "\n";        
        //out << "#T muB muQ muS  "
        //       "J30 J32  N20_b N20_q N20_s  "
        //       "M10_bb M10_bq M10_bs M10_qq M10_qs M10_ss  "
        //       "M11_bb M11_bq M11_bs M11_qq M11_qs M11_ss\n";

        for (long long idx = 0; idx < total; ++idx) {
            const auto& r = rows_ce[(size_t)idx];
            out << r.muB << " " << r.muQ << " " << r.muS << " " << r.T << " "
                << r.J30 << " " << r.J32 << " "
                << r.N20_b << " " << r.N20_q << " " << r.N20_s << " "
                << r.M10_bb << " " << r.M10_bq << " " << r.M10_bs << " "
                << r.M10_qq << " " << r.M10_qs << " " << r.M10_ss << " "
                << r.M11_bb << " " << r.M11_bq << " " << r.M11_bs << " "
                << r.M11_qq << " " << r.M11_qs << " " << r.M11_ss
                << "\n";
        }
    }

    std::cout << "Done. Wrote grad_coeffs.dat and ce_moments.dat\n";
    return 0;
}
