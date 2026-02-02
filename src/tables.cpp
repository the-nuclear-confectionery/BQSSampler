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





void Table4D::load_from_file(const std::string& filename,
                             int total_value_columns,
                             const std::vector<int>& columns_to_keep)
{
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open " + filename);
    }

    std::string header;
    long long Npoints;
    file >> header >> Nx >> Ny >> Nz >> Nw >> Npoints;

    if (Nx <= 1 || Ny <= 1 || Nz <= 1 || Nw <= 1)
        throw std::runtime_error("All dimensions must be > 1");

    if (Npoints != 1LL * Nx * Ny * Nz * Nw)
        throw std::runtime_error("Npoints mismatch");

    if (total_value_columns <= 0)
        throw std::runtime_error("total_value_columns must be > 0");

    // ------------------------------
    // Decide which columns to keep
    // ------------------------------
    std::vector<int> keep = columns_to_keep;
    if (keep.empty()) {
        keep.resize(total_value_columns);
        for (int c = 0; c < total_value_columns; ++c) keep[c] = c; // keep all
    }

    for (int idx : keep) {
        if (idx < 0 || idx >= total_value_columns)
            throw std::runtime_error("columns_to_keep index out of range");
    }

    n_fields = static_cast<int>(keep.size());

    x_vals.resize(Nx);
    y_vals.resize(Ny);
    z_vals.resize(Nz);
    w_vals.resize(Nw);

    data.resize(static_cast<std::size_t>(Npoints) * n_fields);

    std::vector<double> tmp(total_value_columns);

    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            for (int k = 0; k < Nz; ++k)
                for (int l = 0; l < Nw; ++l) {

                    double x,y,z,w;
                    file >> x >> y >> z >> w;

                    for (int c = 0; c < total_value_columns; ++c)
                        file >> tmp[c];
                    if (j==0 && k==0 && l==0) x_vals[i]=x;
                    if (i==0 && k==0 && l==0) y_vals[j]=y;
                    if (i==0 && j==0 && l==0) z_vals[k]=z;
                    if (i==0 && j==0 && k==0) w_vals[l]=w;

                    const int p = ((i * Ny + j) * Nz + k) * Nw + l;

                    for (int f = 0; f < n_fields; ++f) {
                        data[p*n_fields + f] = tmp[ keep[f] ];
                    }
                }

    x_min = x_vals.front(); x_max = x_vals.back();
    y_min = y_vals.front(); y_max = y_vals.back();
    z_min = z_vals.front(); z_max = z_vals.back();
    w_min = w_vals.front(); w_max = w_vals.back();

    dx = (x_max - x_min) / (Nx - 1);
    dy = (y_max - y_min) / (Ny - 1);
    dz = (z_max - z_min) / (Nz - 1);
    dw = (w_max - w_min) / (Nw - 1);
}


double Table4D::interpolate(double xP, double yP,
                            double zP, double wP,
                            int field) const
{

    //std::cout << "Interpolating field " << field << " at point ("
    //          << xP << ", " << yP << ", " << zP << ", " << wP << ")\n";
    //std::cout << "Table ranges: x[" << x_min << ", " << x_max << "], "
    //          << "y[" << y_min << ", " << y_max << "], "
    //          << "z[" << z_min << ", " << z_max << "], "
    //          << "w[" << w_min << ", " << w_max << "]\n";

    if (field < 0 || field >= n_fields)
        throw std::out_of_range("field out of range");

    // Bounds check
    if (xP < x_min || xP > x_max ||
        yP < y_min || yP > y_max ||
        zP < z_min || zP > z_max ||
        wP < w_min || wP > w_max)
        throw std::out_of_range("point out of bounds");

    // ------------------------------------------------------------
    // (1) O(1) index computation for uniform grids
    // i = floor((xP - x_min)/dx) clamped to [0, Nx-2]
    // ------------------------------------------------------------
    int i = static_cast<int>(std::floor((xP - x_min) / dx));
    int j = static_cast<int>(std::floor((yP - y_min) / dy));
    int k = static_cast<int>(std::floor((zP - z_min) / dz));
    int l = static_cast<int>(std::floor((wP - w_min) / dw));

    if (i < 0) i = 0; if (i > Nx-2) i = Nx-2;
    if (j < 0) j = 0; if (j > Ny-2) j = Ny-2;
    if (k < 0) k = 0; if (k > Nz-2) k = Nz-2;
    if (l < 0) l = 0; if (l > Nw-2) l = Nw-2;

    int ip = i+1, jp = j+1, kp = k+1, lp = l+1;

    // ------------------------------------------------------------
    // (2) Fractional coordinates inside the cell
    // Use x0 = x_min + i*dx 
    // ------------------------------------------------------------
    const double x0 = x_min + i * dx;
    const double y0 = y_min + j * dy;
    const double z0 = z_min + k * dz;
    const double w0 = w_min + l * dw;

    const double tx = (xP - x0) / dx;
    const double ty = (yP - y0) / dy;
    const double tz = (zP - z0) / dz;
    const double tw = (wP - w0) / dw;

    const double ax=1-tx, bx=tx;
    const double ay=1-ty, by=ty;
    const double az=1-tz, bz=tz;
    const double aw=1-tw, bw=tw;

    // ------------------------------------------------------------
    // (3) Explicit strides + accessor (like your LI_4D)
    // ------------------------------------------------------------
    const int stride_w  = 1;
    const int stride_z  = Nw;
    const int stride_y  = Nz * Nw;
    const int stride_x  = Ny * Nz * Nw;

    auto F = [&](int ii,int jj,int kk,int ll)->double {
        const int p = ii*stride_x + jj*stride_y + kk*stride_z + ll*stride_w;
        return data[p*n_fields + field];
    };

    // ------------------------------------------------------------
    // (4) 16-corner multilinear interpolation
    // ------------------------------------------------------------
    double result = 0.0;

    result += ax*ay*az*aw * F(i ,j ,k ,l );
    result += bx*ay*az*aw * F(ip,j ,k ,l );
    result += ax*by*az*aw * F(i ,jp,k ,l );
    result += bx*by*az*aw * F(ip,jp,k ,l );
    result += ax*ay*bz*aw * F(i ,j ,kp,l );
    result += bx*ay*bz*aw * F(ip,j ,kp,l );
    result += ax*by*bz*aw * F(i ,jp,kp,l );
    result += bx*by*bz*aw * F(ip,jp,kp,l );

    result += ax*ay*az*bw * F(i ,j ,k ,lp);
    result += bx*ay*az*bw * F(ip,j ,k ,lp);
    result += ax*by*az*bw * F(i ,jp,k ,lp);
    result += bx*by*az*bw * F(ip,jp,k ,lp);
    result += ax*ay*bz*bw * F(i ,j ,kp,lp);
    result += bx*ay*bz*bw * F(ip,j ,kp,lp);
    result += ax*by*bz*bw * F(i ,jp,kp,lp);
    result += bx*by*bz*bw * F(ip,jp,kp,lp);

    return result;
}
