#ifndef SAMPLER_H
#define SAMPLER_H

#include <random>
#include <iostream>
#include <chrono>
#include <random>

#include "surface.h"
#include "particle.h"
#include "particle_system.h"
#include "lrf.h"


class Sampler {
public:
    // Constructor that takes a Surface object reference
    Sampler();


    // Method to sample 
    void sample(ParticleSystem& particle_system, const Surface& surface, const NumericalIntegrator& integrator);


private:

    unsigned sampler_seed;  // Store the random number
};



#endif // SAMPLER_H 