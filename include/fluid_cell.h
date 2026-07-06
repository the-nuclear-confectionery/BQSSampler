#ifndef FLUID_CELL_H
#define FLUID_CELL_H

// ======================================================
// FluidCell
// ------------------------------------------------------
// Self-contained description of a single freeze-out cell.
// Holds every per-cell scalar the sampler reads when drawing
// particles, so external code can sample from one cell without
// constructing a full Surface object.
//
// All fields default to 0.0 so callers can brace-initialize only
// the quantities they have (e.g. an ideal cell leaves the
// dissipative members at zero).
// ======================================================
struct FluidCell {
    // --- spacetime coordinates ---
    double tau = 0.0, x = 0.0, y = 0.0, eta = 0.0;

    // --- fluid velocity (Milne components) ---
    double ut = 0.0, ux = 0.0, uy = 0.0, ueta = 0.0;

    // --- surface element ---
    double dsigma_t = 0.0, dsigma_x = 0.0, dsigma_y = 0.0, dsigma_eta = 0.0;

    // --- thermodynamics ---
    double T = 0.0, muB = 0.0, muS = 0.0, muQ = 0.0;
    double E = 0.0, P = 0.0, s = 0.0;   // energy density, pressure, entropy density

    // --- shear stress tensor (Milne components) ---
    double shv_tt = 0.0, shv_tx = 0.0, shv_ty = 0.0, shv_teta = 0.0;
    double shv_xx = 0.0, shv_xy = 0.0, shv_xeta = 0.0;
    double shv_yy = 0.0, shv_yeta = 0.0;
    double shv_etaeta = 0.0;

    // --- bulk pressure ---
    double bulk = 0.0;

    // --- diffusion currents ---
    double diff_B0 = 0.0, diff_Bx = 0.0, diff_By = 0.0, diff_Beta = 0.0;
    double diff_S0 = 0.0, diff_Sx = 0.0, diff_Sy = 0.0, diff_Seta = 0.0;
    double diff_Q0 = 0.0, diff_Qx = 0.0, diff_Qy = 0.0, diff_Qeta = 0.0;
};

#endif // FLUID_CELL_H
