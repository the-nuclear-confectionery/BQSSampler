#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <iostream>
#include <cmath>
#include "numerical_integrator.h"
#include "surface.h"



class Particle {
public:
    int pid;                 // Particle ID
    double mass;             // Mass of the particle
    double E;                // Energy of the particle
    double px;               // x-component of the momentum
    double py;               // y-component of the momentum
    double pz;               // z-component of the momentum
    double t;                // Time
    double x;                // x-coordinate
    double y;                // y-coordinate
    double z;                // z-coordinate

    // Constructor
    Particle(int pid, double mass, double E, double px, double py, double pz, double t, double x, double y, double z);

};

#endif // PARTICLE_H  