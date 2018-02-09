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
    _pop[2] = 0.1; //0.1
}
int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, ofstream &fileLyapunov){
    int nIteration(static_cast<int>((tFinal-tRecord)/tStep)); // number of steps
    cout << "ok" << endl;
    //DEFINITIONS FOR GNU LIBRARY
    size_t dimension(_dim+_Dim*_Dim);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // integration method, type of the stepping function
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimension); // creation an instance of the stepping function
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0); // object keeping the local error on each step within an absolute error of 1e-6 and relative error of 0.0
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (dimension); // instance of an evolution function for a system of 1 dimensions
    gsl_odeiv_system sys = {func, NULL, dimension, this}; // the Jacobian is useless with this method : use NULL instead

/////////////////
// INTEGRATION //
/////////////////

    double y[_dim+_Dim*_Dim]; // array with the biomass
    for (int i=0; i<3; i++){
        y[i]=_pop[i]; // initialisation of biomasses
        y[i+3]=_x[i]; // initialisation of metabolic rates
    }
    for (int i=_dim; i<_dim+_Dim*_Dim; i++){
        y[i] = 0;
    }
    // Structures for Lyapunov exponent calculation
    double* norm = new double[_Dim]; // vector with the growth rate of the orthonormal vectors
    double* le = new double[_Dim]; // vector with the Lyapunov exponents (Lyapunov spectrum)
    for (int i=0; i<_Dim; i++){ // initialisation
        norm[i] = 0;
        le[i] = 0;
    }
    gsl_matrix * vectGSR = gsl_matrix_calloc (_Dim, _Dim); // matrix containing the vectors for the GSR

    double nPoint(0); // counter of the number of recorded points
    double tPoint(0); // points when data are recorded
    _record = false;

    while (t < tRecord){

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
    }
    tPoint = t + tStep;
    _record = true;

    // Initialise the matrix for the calculation of the Lyapunov exponent
    gsl_matrix_set_identity (vectGSR); // initialise the tangent trajectory
    for (int i=0; i<_Dim; i++){
        for (int j=0; j<_Dim; j++){
            y[_dim+i*_Dim+j] = gsl_matrix_get(vectGSR, i, j); // initialise the tangent trajectory in the y[] vector
        }
    }

///////////////
// RECORDING //
///////////////

    for (int i=0; i<nIteration; i++){
        gsl_odeiv_step_reset (s); // reset for the new integration
        gsl_odeiv_evolve_reset (e);
        // Integration
        while (t < tPoint){
            int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tPoint, &h, y); // integration
            if (status != GSL_SUCCESS)
            break;
        }
        // Matrix containing the final vectors calculated by the linearised system
            for (int i=0; i<_Dim; i++){
                for (int j=0; j<_Dim; j++){
                    gsl_matrix_set(vectGSR, i, j, y[_dim+i*_Dim+j]); // update the matrix
                }
            }
            GSR(_Dim, vectGSR, norm); // GSR - Gram-Schmidt reorthonormalization
            for (int i=0; i<_Dim; i++){
                for (int j=0; j<_Dim; j++){
                    y[_dim+i*_Dim+j] = gsl_matrix_get(vectGSR, i, j); // update the vector y[] with the reorthonormalised values
                }
            }
            // Norm for the calculation of the Lyapunov exponent
            for (int j=0; j<_Dim; j++){
                le[j] += log(norm[j]); // sum of the growth rates
            }

            tPoint += tStep ;
            nPoint++;
        }

        //fileChronic << t << ";" << y[0] << ";" << y[1] << ";" << y[2] << endl; // write the values of PP, SP and t

        nPoint++;
        tPoint = t + tStep;

    // Calculation of the Lyapunov spectrum
    for (int i=0; i<_Dim; i++){
        le[i] /= (nPoint*tStep);
    }
    fileLyapunov << _K << ";" << _X ; // write the values of PP, SP and t
    for (int i=0; i<_Dim; i++){
        fileLyapunov << ";" << le[i]; // write the values of biomasses and metabolic rates
    }
    fileLyapunov << endl;

    // END !!!
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    gsl_matrix_free(vectGSR);
    delete norm;
    delete le;

    return nPoint;
}
void Community::jacobian(const double y[], double f[]){
    // y{B1,B2,B3,x1,x2,x3}
    // y{ 0, 1, 2, 3, 4, 5}
    // Creation of multi-used terms
    _element = 0; // _element of the matrix
    _denomB2 = 0; // Denominator of the functional response of the herbivore
    _denomB3 = 0; // Denominator of the functional response of the carnivore
    _term1 = 0; // loss B1B1
    _term2 = 0; // loss B1B2
    _term3 = 0; // loss B2B2
    _term4 = 0; // loss B2B3
    _term5 = 0; // loss B1x2
    _term6 = 0; // loss B2x3
    // Calculation of the terms
    _denomB2 = _B0 + _c * y[1] * _B0 + y[0];
    _denomB3 = _B0 + _c * y[2] * _B0 + y[1];
    _term1 = y[4] * _y * _B0 * (1 + _c * y[1]) / pow(_denomB2,2);
    _term2 = y[4] * _y * y[0] * (_B0 + y[0]) / pow(_denomB2,2);
    _term3 = y[5] * _y * _B0 * (1 + _c * y[2]) / pow(_denomB3,2);
    _term4 = y[5] * _y * y[1] * (_B0 + y[1]) / pow(_denomB3,2);
    _term5 = _y * y[0] * y[1] / _denomB2;
    _term6 = _y * y[1] * y[2] / _denomB3;

    // (0,0) - dB1/dB1
    _element = _r * (_K - 2 * y[0]) / _K - y[1] * _term1 / _e[0];
    gsl_matrix_set (_jac, 0, 0, _element);
    // (0,1) - dB1/dB2
    _element = - _term2 / _e[0];
    gsl_matrix_set (_jac, 0, 1, _element);
    // (0,2) - dB1/dB3
    _element = 0;
    gsl_matrix_set (_jac, 0, 2, _element);
    // (0,3) - dB1/dx2
    _element = - _term5 / _e[0];
    gsl_matrix_set (_jac, 0, 3, _element);
    // (0,4) - dB1/dx3
    _element = 0;
    gsl_matrix_set (_jac, 0, 4, _element);

    // (1,0) - dB2/dB1
    _element = y[1] * _term1;
    gsl_matrix_set (_jac, 1, 0, _element);
    // (1,1) - dB2/dB2
    _element = _term2 - y[2] * _term3 / _e[1] - y[4];
    gsl_matrix_set (_jac, 1, 1, _element);
    // (1,2) - dB2/dB3
    _element = - _term4 / _e[1];
    gsl_matrix_set (_jac, 1, 2, _element);
    // (1,3) - dB2/dx2
    _element = _term5 - y[1];
    gsl_matrix_set (_jac, 1, 3, _element);
    // (1,4) - dB2/dx3
    _element = - _term6 / _e[1];
    gsl_matrix_set (_jac, 1, 4, _element);

    // (2,0) - dB3/dB1
    _element = 0;
    gsl_matrix_set (_jac, 2, 0, _element);
    // (2,1) - dB3/dB2
    _element = y[2] * _term3;
    gsl_matrix_set (_jac, 2, 1, _element);
    // (2,2) - dB3/dB3
    _element = _term4 - y[5];
    gsl_matrix_set (_jac, 2, 2, _element);
    // (2,3) - dB3/dx2
    _element = 0;
    gsl_matrix_set (_jac, 2, 3, _element);
    // (2,4) - dB3/dx3
    _element = _term6 - y[2];
    gsl_matrix_set (_jac, 2, 4, _element);

    // (3,0) - dx2/dB1
    _element = _X * _term1;
    gsl_matrix_set (_jac, 3, 0, _element);
    // (3,1) - dx2/dB2
    _element = - _c * y[4] * _X * _y * y[0] * _B0 / pow(_denomB2,2);
    gsl_matrix_set (_jac, 3, 1, _element);
    // (3,2) - dx2/dB3
    _element = 0;
    gsl_matrix_set (_jac, 3, 2, _element);
    // (3,3) - dx2/dx2
    _element = _X * (-1 + _y * y[0] / _denomB2);
    gsl_matrix_set (_jac, 3, 3, _element);
    // (3,4) - dx2/dx3
    _element = 0;
    gsl_matrix_set (_jac, 3, 4, _element);

    // (4,0) - dx3/dB1
    _element = 0;
    gsl_matrix_set (_jac, 4, 0, _element);
    // (4,1) - dx3/dB2
    _element = _X * _term3;
    gsl_matrix_set (_jac, 4, 1, _element);
    // (4,2) - dx3/dB3
    _element = - _c * y[5] * _X * _y * y[1] * _B0 / pow(_denomB3,2);
    gsl_matrix_set (_jac, 4, 2, _element);
    // (4,3) - dx3/dx2
    _element = 0;
    gsl_matrix_set (_jac, 4, 3, _element);
    // (4,4) - dx3/dx3
    _element = _X * (-1 + _y * y[1] / _denomB3);
    gsl_matrix_set (_jac, 4, 4, _element);

    // fill in the trajectory matrix
    for (int i=0; i<_Dim; i++){
        for (int j=0; j<_Dim; j++){
            gsl_matrix_set(_vect, i, j, y[_dim+i*_Dim+j]);
        }
    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _jac, _vect, 0.0, _prod); // matrix-matrix product - jac*vect

    // return the variations of the linearised system
    for (int i=0; i<_Dim; i++){
        for (int j=0; j<_Dim; j++){
            f[_dim+i*_Dim+j] = gsl_matrix_get(_prod, i, j);
        }
    }
}

// constructor

Community::Community(
                    int dim,
                    int Dim,
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
                    _dim(dim)
                    ,_Dim(Dim)
                    ,_r(ar)
                    ,_y(ay)
                    ,_c(c)
                    ,_B0(B0)
                    ,_FR(FR)
                    ,_K(K)
                    ,_X(X)
                    ,_Xmin(Xmin)
                    ,_Xmax(Xmax){
                        // Jacobian
                       _jac = new gsl_matrix;
                       _jac = gsl_matrix_calloc (Dim, Dim);
                       // matrix of vectors
                       _vect = new gsl_matrix;
                       _vect = gsl_matrix_calloc (Dim, Dim);
                       //gsl_matrix_set_identity(_vect); // _vect is the identity matrix (orthogonal tangent vectors)
                       // results of the jac*vect product
                       _prod = new gsl_matrix;
                       _prod = gsl_matrix_calloc (Dim, Dim);
                    }

// destructor
Community::~Community(){
    gsl_matrix_free(_jac); // free the memory
    gsl_matrix_free(_vect);
    gsl_matrix_free(_prod);
}
