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

    for(int i=0; i<c->_dim+c->_Dim*c->_Dim; i++){
        f[i] = 0;
    }

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

    // Linearised system
    if (c->_record){
        c->jacobian(y,f);
    }

    return GSL_SUCCESS;
}

void GSR(int dim, gsl_matrix *m, double norm[]){
    gsl_matrix * mGSR = gsl_matrix_calloc (dim, dim); // storage matrix
    gsl_matrix_memcpy(mGSR, m); // mGSR is a copy of m
    double scalProd(0); // scalar product of columns

    for (int i=0; i<dim; i++){ // columns
        for (int j=0; j<i; j++){ // columns of the precedent vectors
            scalProd = 0; // resetting
            for (int k=0; k<dim; k++){ // rows
                scalProd += gsl_matrix_get(m, k, i) * gsl_matrix_get(mGSR, k, j); // scalar product
            }
            for (int k=0; k<dim; k++){
                gsl_matrix_set(mGSR, k, i, gsl_matrix_get(mGSR, k, i) - scalProd * gsl_matrix_get(mGSR, k, j)); // remove the j component of the vector i
            }
        }
        norm[i] = 0;
        for (int k=0; k<dim; k++){ // rows
            norm[i] += pow(gsl_matrix_get(mGSR, k, i),2); // calculate the norm of the final vector
        }
        norm[i] = sqrt(norm[i]);
        for (int k=0; k<dim; k++){ // rows
            gsl_matrix_set(mGSR, k, i, gsl_matrix_get(mGSR, k, i)/norm[i]); // divide the vector by its norm
        }
    }
    gsl_matrix_memcpy(m, mGSR);
    gsl_matrix_free(mGSR); // Free the memory
}
