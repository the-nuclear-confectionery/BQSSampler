#ifndef TOOLS_H
#define TOOLS_H
#include <iostream>
#include <vector>



class Tools {
public:
    static void read_2d_table(
        const std::string& filename,
        std::vector<double>& x_vals,
        std::vector<double>& y_vals,
        std::vector<double>& f_vals
    );

    static double bilinear_interpolate(
        const std::vector<double>& x_vals,
        const std::vector<double>& y_vals,
        const std::vector<double>& f_vals,
        double x, double y);
};

#endif
