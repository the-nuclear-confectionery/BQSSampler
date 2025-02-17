#ifndef _GLOBALS_H
#define _GLOBALS_H

#define HBARC 0.1973269718
#define TOL 1e-10

enum dir { tau = 0, x, y, eta, NUM_POSITION = 4 };
enum thermo { E = 0, T, P, s, NUM_thermo = 4 };
enum charge { B = 0, S, Q, NUM_CHARGES = 3 };
enum shv { tt = 0, tx, ty, teta, xx, xy, xeta, yy, yeta, etaeta, NUM_STRESS_TENSOR = 10 };



#endif // _GLOBALS_H