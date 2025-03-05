#ifndef TABLES_H
#define TABLES_H

#include <vector>
#include <string>


class Table2D {
private:
    std::vector<double> x_vals;             
    std::vector<double> y_vals;             
    std::vector<std::vector<double>> f_data; 

    int Nx, Ny;   // Grid dimensions
    double x_min, x_max, y_min, y_max; // Min/max values for fast access
    double dx, dy; // Grid spacing

public:
    Table2D() = default;
    explicit Table2D(const std::string& filename);

    void load_from_file(const std::string& filename);
    double bilinear_interpolate(double x, double y) const;
};

#endif
