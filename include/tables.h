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




class Table4D {
private:
    std::vector<double> x_vals, y_vals, z_vals, w_vals;
    std::vector<double> data;   // flat storage

    int Nx=0, Ny=0, Nz=0, Nw=0;
    int n_fields=0;

    double x_min=0, x_max=0;
    double y_min=0, y_max=0;
    double z_min=0, z_max=0;
    double w_min=0, w_max=0;

    double dx=0, dy=0, dz=0, dw=0;

public:
    Table4D() = default;

    void load_from_file(const std::string& filename,
                        int total_value_columns,
                        const std::vector<int>& columns_to_keep = {});

    double interpolate(double xP, double yP,
                       double zP, double wP,
                       int field) const;
};

#endif
