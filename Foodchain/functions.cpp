#include "functions.h"
#include "community.h"

#include "math.h"
#include <iostream>
#include <fstream>
#include "string"
using namespace std;

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration

double setAllometric(double a, double M, double s){
    return (a*pow(M,s));
}

 // INTEGRATION
 // System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params){
    Community *c = (Community *) params;

    for(int i=0; i<6; i++){
        f[i] = 0;
    }
    c->_flux[0] = 0;
    c->_flux[1] = 0;

    // PRIMARY PRODUCER
    f[0] = c->_r * (c->_K - y[0]) / c->_K // plant growth
            - y[4] * c->_y * pow(y[0],c->_FR) * y[1] / (c->_B0 + c->_c * y[1] * c->_B0 + pow(y[0],c->_FR))/c->_e[0]; // predation

    // HERBIVORE
    f[1] =  y[4] * c->_y * pow(y[0],c->_FR) * y[1] / (c->_B0 + c->_c * y[1] * c->_B0 + pow(y[0],c->_FR)) // predation
            - y[5] * c->_y * pow(y[1],c->_FR) * y[2] / (c->_B0 + c->_c * y[2] * c->_B0 + pow(y[1],c->_FR))/c->_e[1] // predation
            - y[4] * y[1]; // respiration

    // CARNIVORE
    f[2] =  y[5] * c->_y * pow(y[1],c->_FR) * y[2] / (c->_B0 + c->_c * y[2] * c->_B0 + pow(y[1],c->_FR)) // predation
            - y[5] * y[2]; // respiration

    // PRIMARY PRODUCER METABOLIC RATE
    f[3] = 0;

    // HERBIVORE METABOLIC RATE
    f[4] =  c->_X * y[4]
            * (c->_y * pow(y[0],c->_FR) / (c->_B0 + c->_c * y[1] * c->_B0 + pow(y[0],c->_FR)) // predation
            - 1 ); // respiration

    // CARNIVORE METABOLIC RATE
    f[5] =  c->_X * y[5]
            * (c->_y * pow(y[1],c->_FR) / (c->_B0 + c->_c * y[2] * c->_B0 + pow(y[1],c->_FR)) // predation
            - 1 ); // respiration

    // PRIMARY PRODUCTION
    c->_flux[0] = c->_r * (c->_K - y[0]) / c->_K;

    // SECONDARY PRODUCTION
    c->_flux[1] = y[4] * c->_y * pow(y[0],c->_FR) * y[1] / (c->_B0 + c->_c * y[1] * c->_B0 + pow(y[0],c->_FR))
                  + y[5] * c->_y * pow(y[1],c->_FR) * y[2] / (c->_B0 + c->_c * y[2] * c->_B0 + pow(y[1],c->_FR));

    return GSL_SUCCESS;
}
