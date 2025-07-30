#include "tables.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <tuple>

Table2D::Table2D(const std::string& filename) {
    load_from_file(filename);
}

void Table2D::load_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error: Could not open " + filename);
    }

    // Read first line to get Nx, Ny, and Npoints
    std::string header;
    int Npoints;
    file >> header >> Nx >> Ny >> Npoints;

    if (Nx <= 0 || Ny <= 0 || Npoints != Nx * Ny) {
        throw std::runtime_error("Error: Invalid grid dimensions in file.");
    }

    // Resize vectors based on Nx and Ny
    x_vals.resize(Nx);
    y_vals.resize(Ny);
    f_data.assign(Nx, std::vector<double>(Ny, 0.0));  // Store f_data[i][j]

    // Read data while ensuring the first column (`x`) remains fixed for each block
    for (int i = 0; i < Nx; ++i) {   // Loop over x values first (blocks)
        for (int j = 0; j < Ny; ++j) { // Then loop over y values
            double x, y, f;
            file >> x >> y >> f;

            if (j == 0) x_vals[i] = x;  // Store `x` once per block
            if (i == 0) y_vals[j] = y;  // Store `y` only on the first iteration of x

            f_data[i][j] = f;  // Store `f` in the correct position
        }
    }
    file.close();

    // Store min/max values for fast access
    x_min = x_vals.front();
    x_max = x_vals.back();
    y_min = y_vals.front();
    y_max = y_vals.back();

    // Compute uniform grid spacing
    dx = (x_max - x_min) / (Nx - 1);
    dy = (y_max - y_min) / (Ny - 1);
}

double Table2D::bilinear_interpolate(double x, double y) const {
    if (x < x_min || x > x_max || y < y_min || y > y_max) {
        std::cerr << "Warning: Interpolation point out of bounds. "
                  <<  "x: " << x << ", y: " << y
                  << ", valid range x: [" << x_min << ", " << x_max << "] y [" << y_min << ", " << y_max << "]\n";
        // Optionally, 
        throw std::out_of_range("Interpolation point out of table bounds.");
    }

    // Compute integer indices using floor
    int iL = std::floor((x - x_min) / dx);
    int jL = std::floor((y - y_min) / dy);
    int iR = iL + 1;
    int jR = jL + 1;

    // Ensure indices are within bounds
    if (iL < 0 || iR >= Nx || jL < 0 || jR >= Ny) {
        throw std::out_of_range("Index out of range in bilinear interpolation.");
    }

    // Access the required points
    double f_LL = f_data[iL][jL]; // Bottom-left
    double f_RL = f_data[iR][jL]; // Bottom-right
    double f_LR = f_data[iL][jR]; // Top-left
    double f_RR = f_data[iR][jR]; // Top-right

    // Compute interpolation weights
    double xL = x_min + iL * dx;
    double xR = xL + dx;
    double yL = y_min + jL * dy;
    double yR = yL + dy;

    // Bilinear interpolation formula
    double result = ((f_LL * (xR - x) + f_RL * (x - xL)) * (yR - y) +
                     (f_LR * (xR - x) + f_RR * (x - xL)) * (y - yL)) /
                    (dx * dy);

    return result;
}
