// main_cf21d_trap_pion.cpp
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <atomic>
#include <chrono>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "settings.h"
#include "surface.h"
#include "particle_system.h"

// ------------------------------ helpers ------------------------------
static inline double sqr(double x) { return x * x; }

// Bose/Fermi distribution:
// sign = -1 for bosons, +1 for fermions
static inline double feq(double arg_over_T, double sign) {
    const double ex = std::exp(arg_over_T);
    //skip if ex + sign <= 0 to avoid numerical issues (can happen for bosons at large negative arg_over_T)
    if (ex + sign <= 0.0) return 0.0;
    else return 1.0 / (ex + sign);
}

// IMPORTANT: you asked to include the g factor in f*fbar.
// Here we define:
//   fbar = 1 - sign*f
//   f*fbar (with g) = g * f * (1 - sign*f)
static inline double ffbar_with_g(double f, double sign, double g) {
    return  f * (1.0 - sign * f );
}

// Trapezoid integration on [a,b] with N intervals => N+1 points.
// returns two integrals simultaneously:
//   I0 = ∫ f(eta) d eta
//   I1 = ∫ cosh(eta - y) f(eta) d eta
template <class F>
static inline void trap_eta_integrals(
    double a, double b, int N, double y,
    const F& f_of_eta,
    double& I0, double& I1
) {
    const double h = (b - a) / double(N);
    I0 = 0.0;
    I1 = 0.0;

    for (int k = 0; k <= N; ++k) {
        const double eta = a + h * double(k);
        const double w = (k == 0 || k == N) ? 0.5 * h : h;

        const double f = f_of_eta(eta);
        I0 += w * f;
        I1 += w * std::cosh(eta - y) * f;
    }
}

// ------------------------------ main physics ------------------------------
// Computes d^3N/(dy dpT^2 dphi) on a (pT,phi) grid for one particle species.
// Parallelized over (pT,phi) tasks using OpenMP.
//
// Adds shear δf:
//   δf_shear = (f * fbar) * (π^{μν} p_μ p_ν) / (2 s)
// with contraction as in your screenshot:
//   π^{μν} p_μ p_ν =
//     m^2[ π^{00} cosh^2(η-y) + τ^4 π^{33} sinh^2(η-y) ]
//     + px^2 π^{11} + py^2 π^{22} + 2 px py π^{12}
static std::vector<std::vector<double>> compute_dN_dy_dpT2_dphi_trap_eta(
    const Surface& surface,
    const double m,
    const double g,
    const double sign,
    const double B,
    const double S,
    const double Q,
    const Settings& cfg,
    double y,
    const std::vector<double>& pT_grid,
    const std::vector<double>& phi_grid
) {
    // eta integration settings
    const double eta_min = cfg.has_key("eta_min") ? cfg.get_double("eta_min") : -5.0;
    const double eta_max = cfg.has_key("eta_max") ? cfg.get_double("eta_max") : +5.0;
    const int N_eta = cfg.has_key("n_eta_trap") ? cfg.get_int("n_eta_trap") : 40;

    // enable/disable shear delta-f from config
    const bool use_deltaf_shear = true;

    // Prefactor still has g as usual (you already had this)
    const double pref = g / std::pow(2.0 * M_PI * HBARC, 3);

    std::vector<std::vector<double>> out(
        pT_grid.size(), std::vector<double>(phi_grid.size(), 0.0)
    );

    // ----- progress settings -----
    const long long total_tasks = (long long)pT_grid.size() * (long long)phi_grid.size();
    std::atomic<long long> done{0};

    const long long progress_step_cfg =
        cfg.has_key("progress_step") ? (long long)cfg.get_int("progress_step") : -1;
    const long long progress_step =
        (progress_step_cfg > 0) ? progress_step_cfg : std::max(1LL, total_tasks / 100);

#ifdef _OPENMP
    const int nthreads_cfg = cfg.has_key("omp_threads") ? cfg.get_int("omp_threads") : 0;
    if (nthreads_cfg > 0) omp_set_num_threads(nthreads_cfg);
#endif

    const auto t0 = std::chrono::steady_clock::now();

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic, 1)
#endif
    for (int ipt = 0; ipt < (int)pT_grid.size(); ++ipt) {
        for (int iph = 0; iph < (int)phi_grid.size(); ++iph) {

            const double pT  = pT_grid[(size_t)ipt];
            const double mT  = std::sqrt(m * m + pT * pT);

            const double phi = phi_grid[(size_t)iph];
            const double px  = pT * std::cos(phi);
            const double py  = pT * std::sin(phi);

            double sum_cells = 0.0;

            for (int i = 0; i < surface.npoints; ++i) {

                // dsigma_mu (covariant components in your storage)
                const double dsig_tau = surface.dsigma_t[i];
                const double dsig_x   = surface.dsigma_x[i];
                const double dsig_y   = surface.dsigma_y[i];

                const double T   = surface.T[i];

                // particle chemical potential mu = B*muB + S*muS + Q*muQ
                const double mu = B * surface.muB[i] + S * surface.muS[i] + Q * surface.muQ[i];
                // (for now kept zero like your original)
                //const double mu = 0.0;
                
                // flow (assumed boost-invariant u^tau, u^x, u^y here)
                const double ut = surface.ut[i];
                const double ux = surface.ux[i];
                const double uy = surface.uy[i];

                // shear + entropy density + tau for this cell
                // Using your Surface variable names from earlier:
                //   π^{00} -> pitt
                //   π^{11} -> pixx
                //   π^{22} -> piyy
                //   π^{12} -> pixy
                //   π^{33} -> pinn   (Milne ηη component)
                const double tau_i = surface.tau[i];
                const double s_ent = surface.s[i];

                const double pitt = surface.shv_tt[i];
                const double pixx = surface.shv_xx[i];
                const double piyy = surface.shv_yy[i];
                const double pixy = surface.shv_xy[i];
                const double pi33 = surface.shv_etaeta[i];

                // f(eta) includes equilibrium + shear delta-f
                auto f_of_eta = [&](double eta) -> double {

                    const double d  = eta - y;
                    const double ch = std::cosh(d);
                    const double sh = std::sinh(d);

                    // p^tau(eta) = mT cosh(eta - y) = mT cosh(y-eta)
                    const double ptau = mT * ch;

                    // p·u = p^tau u^tau - p^x u^x - p^y u^y
                    const double p_dot_u = ptau * ut - px * ux - py * uy;

                    const double x = (p_dot_u - mu) / T;
                    //skip if mu> p·u to avoid numerical issues (can happen for bosons at large negative x)
                    if (mu > p_dot_u) return 0.0;
                    const double f0 = feq(x, sign);

                    if (!use_deltaf_shear || s_ent == 0.0) {
                        return f0;
                    }

                    // π^{μν} p_μ p_ν (your screenshot formula)
                    const double pip_p =
                        (m * m) * (pitt * ch * ch + std::pow(tau_i, 4) * pi33 * sh * sh)
                        + (px * px) * pixx
                        + (py * py) * piyy
                        + 2.0 * px * py * pixy;

                    // factor f * fbar INCLUDING g (as you requested)
                    const double f0f0bar_g = ffbar_with_g(f0, sign, g);

                    // δf_shear
                    const double deltaf = f0f0bar_g * (pip_p / (2.0 * s_ent));

                    // total distribution
                    //check if pdsigma < 0
                    double pdsigma = ptau * dsig_tau + px * dsig_x + py * dsig_y;
                    if (pdsigma < 0.0) {
                        return 0.0;
                    }
                    //if deltaf< 0, set it to zero
                    if (deltaf < 0.0) return f0;
                    double f = f0 + deltaf;

                    // (optional safety; uncomment if needed)
                    // ==if (f < 0.0) f = 0.0;

                    return f0;
                };

                double I0 = 0.0, Icosh = 0.0;
                trap_eta_integrals(eta_min, eta_max, N_eta, y, f_of_eta, I0, Icosh);

                // CF21D form:
                // mT * dsigma_tau * ∫ cosh(eta-y) f  +  (px dsigma_x + py dsigma_y) * ∫ f
                const double term_tau = mT * dsig_tau * Icosh;
                const double term_xy  = (px * dsig_x + py * dsig_y) * I0;

                sum_cells += (term_tau + term_xy);
            }

            out[(size_t)ipt][(size_t)iph] = pref * sum_cells;

            // ----- progress -----
            const long long local_done = ++done;

#ifdef _OPENMP
            if (omp_get_thread_num() == 0)
#endif
            {
                if (local_done == total_tasks || (local_done % progress_step == 0)) {
                    const auto t1 = std::chrono::steady_clock::now();
                    const double sec = std::chrono::duration<double>(t1 - t0).count();
                    const double pct = 100.0 * (double)local_done / (double)total_tasks;

                    std::cout << std::fixed << std::setprecision(2)
                              << "Progress: " << pct << "%  ("
                              << local_done << "/" << total_tasks
                              << ")  elapsed " << sec << " s\n";
                    std::cout.flush();
                }
            }
        }
    }

    return out;
}

// ------------------------------ executable ------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " config.yaml\n";
        return 1;
    }

    const std::string config_file = argv[1];
    Settings cfg(config_file);

    const std::string tables_path  = cfg.get_string("tables_path");
    const std::string surface_file = cfg.get_string("input_file");

    Surface surface(surface_file, cfg);
    surface.read_data();
    std::cout << "Surface loaded: npoints=" << surface.npoints << "\n";

    ParticleSystem ps;
    ps.read_particle_list(tables_path + "/pdg.dat");

    const int pion_pid = 211;
    int ipion = -1;
    for (int i = 0; i < ps.nparticles; ++i) {
        if (ps.pid[i] == pion_pid) { ipion = i; break; }
    }
    if (ipion < 0) throw std::runtime_error("PID 211 not found in pdg.dat");

    const double m = ps.mass[ipion];
    const double g = ps.spin_degeneracy[ipion];
    const double B = ps.baryon[ipion];
    const double S = ps.strange[ipion];
    const double Q = ps.charge[ipion];
    const double sign = ps.theta[ipion]; // -1 bosons, +1 fermions
    std::cout << "Using pion+ PID=211, m=" << m << ", g=" << g << "\n";

#ifdef _OPENMP
    std::cout << "OpenMP enabled. max_threads=" << omp_get_max_threads() << "\n";
#else
    std::cout << "OpenMP NOT enabled (compiled without -fopenmp). Running serial.\n";
#endif

    // --- grids ---
    const double y = cfg.has_key("y") ? cfg.get_double("y") : 0.0;

    const double pT_min = cfg.has_key("pT_min") ? cfg.get_double("pT_min") : 0.0;
    const double pT_max = cfg.has_key("pT_max") ? cfg.get_double("pT_max") : 3.0;
    const int n_pT      = cfg.has_key("n_pT") ? cfg.get_int("n_pT") : 60;
    const int n_phi     = cfg.has_key("n_phi") ? cfg.get_int("n_phi") : 48;

    std::vector<double> pT_grid(n_pT);
    for (int i = 0; i < n_pT; ++i) {
        const double a = double(i) / double(n_pT - 1);
        pT_grid[i] = pT_min + a * (pT_max - pT_min);
    }

    std::vector<double> phi_grid(n_phi);
    for (int i = 0; i < n_phi; ++i) {
        phi_grid[i] = (2.0 * M_PI) * (double(i) / double(n_phi));
    }

    // --- compute ---
    auto dN = compute_dN_dy_dpT2_dphi_trap_eta(surface, m, g, sign, B, S, Q, cfg, y, pT_grid, phi_grid);

    // --- output ---
    const std::string outname = cfg.has_key("output_file")
        ? cfg.get_string("output_file")
        : "dN_dy_dpT2_dphi_pid" + std::to_string(pion_pid) + ".dat";

    std::ofstream out(outname);
    out << std::setprecision(16);
    out << "# pid " << pion_pid << "  y " << y << "\n";
    out << "# columns: pT  phi  d3N_dy_dpT2_dphi\n";

    for (size_t ipt = 0; ipt < pT_grid.size(); ++ipt) {
        for (size_t iph = 0; iph < phi_grid.size(); ++iph) {
            out << pT_grid[ipt] << " " << phi_grid[iph] << " " << dN[ipt][iph] << "\n";
        }
        out << "\n";
    }

    std::cout << "Wrote: " << outname << "\n";
    return 0;
}
