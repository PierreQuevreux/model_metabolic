#include "community.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
using namespace std;

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_blas.h> // for matrix operations
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// CREATION OF THE FOODWEB
void Community::setSystem(double ax, double s){
    // BODY MASSES
    _M[0] = 1;
    _M[1] = 100;
    _M[2] = 10000;

    // METABOLIC RATE
    _x[0] = 0;
    _x[1] = setAllometric(ax,_M[1]/_M[0],s);
    _x[2] = setAllometric(ax,_M[2]/_M[0],s);

    // ASSIMILATION EFFICIENCY
    _e[0] = 0.45;
    _e[1] = 0.85;
    _e[2] = 0.85;

    // INITIAL POPULATIONS
    _pop[0] = 1;
    _pop[1] = 0.5; //
    _pop[2] = 0.1;
}
int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, ofstream &fileChronic){
    //DEFINITIONS FOR GNU LIBRARY
    size_t dimension(6);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // integration method, type of the stepping function
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimension); // creation an instance of the stepping function
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0); // object keeping the local error on each step within an absolute error of 1e-6 and relative error of 0.0
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (dimension); // instance of an evolution function for a system of 1 dimensions
    gsl_odeiv_system sys = {func, NULL, dimension, this}; // the Jacobian is useless with this method : use NULL instead

/////////////////
// INTEGRATION //
/////////////////

    double y[6]; // array with the biomass
    for (int i=0; i<3; i++){
        y[i]=_pop[i]; // initialisation of biomasses
        y[i+3]=_x[i]; // initialisation of metabolic rates
    }

    double nPoint(0); // counter of the number of recorded points
    double tPoint(tRecord); // points when data are recorded

    while (t < tFinal){
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tFinal, &h, y); // integration
        if (status != GSL_SUCCESS)
        break;

        // Extinction
        for (int i=0; i<3; i++){
            if (y[i]<extinctionThreshold && y[i]!=0){
                y[i]=0; // extinction
                y[i+3]=0; // metabolic rate drops to zero
            }
        }
        // Metabolic rate limitation
        for (int i=3; i<6; i++){
            if (y[i]<_Xmin){
                y[i]=_Xmin; // lower limit
            }
            if (y[i]>_Xmax){
                y[i]=_Xmax; // upper limit
            }
        }

        ///////////////
        // RECORDING //
        ///////////////

        if (t>tPoint){

            // CHRONICS
            //cout << tPoint << endl;
            fileChronic << t << ";" << _K << ";" << _X << ";" << _flux[0] << ";" << _flux[1]; // write the values of PP, SP and t
            for (int i=0; i<6; i++){
                fileChronic << ";" << y[i]; // write the values of biomasses and metabolic rates
            }
            fileChronic << endl;

            tPoint += tStep ;
            nPoint++;
        }
    }

    // END !!!
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    return nPoint;
}

// constructor

Community::Community(
                    double ar,
                    double ay,
                    double c,
                    double B0,
                    double FR,
                    double K,
                    double X,
                    double Xmin,
                    double Xmax)
                     :
                    _r(ar)
                    ,_y(ay)
                    ,_c(c)
                    ,_B0(B0)
                    ,_FR(FR)
                    ,_K(K)
                    ,_X(X)
                    ,_Xmin(Xmin)
                    ,_Xmax(Xmax)
                    {}

// destructor
Community::~Community(){}
