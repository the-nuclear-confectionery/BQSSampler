#include "particle.h"


// Constructor implementation
Particle::Particle(int pid, double mass, double E, double px, double py, double pz, double t, double x, double y, double z){
    this->pid = pid;
    this->mass = mass;
    this->E = E;
    this->px = px;
    this->py = py;
    this->pz = pz;
    this->t = t;
    this->x = x;
    this->y = y;
    this->z = z;  
}