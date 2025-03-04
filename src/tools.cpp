
#include "tools.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>


void Tools::read_2d_table(
    const std::string& filename,
    std::vector<double>& x_vals,
    std::vector<double>& y_vals,
    std::vector<double>& f_vals) {

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }

    double x, y, f;
    while (file >> x >> y >> f) {
        x_vals.push_back(x);
        y_vals.push_back(y);
        f_vals.push_back(f);
    }
    file.close();
}


double Tools::bilinear_interpolate(
    const std::vector<double>& x_vals,
    const std::vector<double>& y_vals,
    const std::vector<double>& f_vals,
    double x, double y) {

    int Nx = std::count_if(x_vals.begin(), x_vals.end(), [&](double val) { return val == x_vals[0]; });
    int Ny = x_vals.size() / Nx;
    double dx = x_vals[Nx] - x_vals[0];
    double dy = y_vals[1] - y_vals[0];

    // Find nearest lower indices
    int i0 = std::max(0, std::min(int((x - x_vals[0]) / dx), Nx - 2));
    int j0 = std::max(0, std::min(int((y - y_vals[0]) / dy), Ny - 2));

    int i1 = i0 + 1;
    int j1 = j0 + 1;

    double x0 = x_vals[i0 * Ny];
    double x1 = x_vals[i1 * Ny];
    double y0 = y_vals[j0];
    double y1 = y_vals[j1];


    int idx_00 = i0 * Ny + j0;
    int idx_01 = i0 * Ny + j1;
    int idx_10 = i1 * Ny + j0;
    int idx_11 = i1 * Ny + j1;

    double Q11 = f_vals[idx_00];
    double Q12 = f_vals[idx_01];
    double Q21 = f_vals[idx_10];
    double Q22 = f_vals[idx_11];

    // Bilinear interpolation
    double denom = (x1 - x0) * (y1 - y0);
    if (denom == 0) return Q11;

    return (Q11 * (x1 - x) * (y1 - y) +
            Q21 * (x - x0) * (y1 - y) +
            Q12 * (x1 - x) * (y - y0) +
            Q22 * (x - x0) * (y - y0)) / denom;
}
