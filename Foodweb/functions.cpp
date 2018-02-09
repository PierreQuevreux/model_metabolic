#include "functions.h"
#include "community.h"
#include "species.h"

#include "math.h"
#include <iostream>
#include <fstream>
#include "string"
#include <sstream> // to use istringstream and convert string to double
using namespace std;

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h> // lib glsblas to put in linkers

double setAllometric(double a, double M, double s){
    return (a*pow(M,s));
}

 // INTEGRATION
 // System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params){
    Community *c = (Community *) params;
    Species **p = c->_params;

    // Initialization to avoid strange values
    double C(0);
    double hPred[c->_diversity]; // array containing the denominator of the functional response for each consumer
    double Kplant(0); // shared growth limitation
    for (int i=0; i<c->_diversity; i++){
        f[i]=0;
        f[c->_diversity + i]=0;
        hPred[i]=0;
    }
    c->_flux[0]=0;
    c->_flux[1]=0;

    // FUNCTIONAL RESPONSE DENOMINATOR
    for (int i=0; i<c->_diversity; i++){ // index for each organism
        if(y[i] != 0 && p[i]->_W >0){ // test if the organism is alive and eats something
            for (int j=0; j<p[i]->_nbPrey; j++){ // index for each prey
                if(y[p[i]->_idPrey[j]] != 0){ // test if the prey is still alive
                    hPred[i] += pow( y[p[i]->_idPrey[j]] , c->_FR );
                }
            }
            hPred[i] /= p[i]->_W; // divide by the uniform preference
            hPred[i] += pow(c->_B0,c->_FR) * ( 1+c->_c*y[i]); // half saturation and predator interference
        }
        else if (p[i]->_type=="Plant"){
            Kplant += y[i];
        }
        else{hPred[i]=1;} // useless but to avoid division by zero...
    }
    Kplant = (c->_K - Kplant)/c->_K; // growth limitation by density dependency

    // PREDATION
    for (int i=0; i<c->_diversity; i++){ // index for each organism
        if (y[i] != 0 && p[i]->_W > 0){ // test if the organism is alive and eats something
            for (int j=0; j<p[i]->_nbPrey; j++){ // index for each prey
                C = y[i+c->_diversity] * c->_y * pow( y[p[i]->_idPrey[j]] , c->_FR ) / ( p[i]->_W * hPred[i]); // biomass uptake of the prey j
                f[i+c->_diversity] += C; // sum of biomass uptake for adaptive metabolic rate calculation
                f[p[i]->_idPrey[j]] -= C * y[i] / p[i]->_prey[j]->_e; // biomass loss for the prey j
                c->_flux[1] += C * y[i]; // secondary production
            }
            f[i] += f[i+c->_diversity] * y[i]; // total biomass gain of the consumer i
        }
    }

    // RESPIRATION, PLANT GROWTH AND METABOLIC RATE
    for (int i=0; i<c->_diversity; i++){
        if (p[i]->_type == "Consumer"){
            f[i] -= y[i] * y[i+c->_diversity] ; // respiration
            f[i+c->_diversity] -= y[i+c->_diversity] ; // metabolic rate
            f[i+c->_diversity] *= c->_X;
        }
        else{
            f[i] += p[i]->_R * y[i] * Kplant; // plant growth
            c->_flux[0] += p[i]->_R * y[i] * Kplant; // primary production
        }
    }

    return GSL_SUCCESS;
}

void trophicLevel(int diversity, string fileMatrixInteraction, double TL[]){
    // Matrix initialisation
    gsl_matrix * Q = gsl_matrix_calloc (diversity,diversity); // transition matrix between non-source compartments (initialised to zero)
    gsl_matrix_set_identity (Q); // set Q as the identity matrix
    gsl_matrix * N = gsl_matrix_calloc (diversity,diversity); //
    gsl_permutation * perm = gsl_permutation_alloc (diversity); // permutation for the inversion
    int s; // signum s (for LU decomposition)


    ifstream file(fileMatrixInteraction.c_str());
    if (file){}
    else {cout << "ERROR: unable to open file at path" << endl;}
    string line;
    string var[diversity]; // string array containing the string form of the values of the matrix of interactions
    string car; // separation of string values in the line
    double value; // converted value ready for the matrix
    int k(0); // counter
    istringstream iss; // variable for the conversion from string to double
    for (int i=0; i<diversity; i++){
        getline(file, line);
        k=0;
        for (int j=0; j<diversity; j++){
            var[j].clear(); // clear the string array
        }
        for (unsigned int j=0; j<line.size(); j++){
            car = line[j];
            if (car!=";"){
                var[k] += car; // write the data letter by letter
            }
            else {k++;}
        }
        for (int j=0; j<diversity; j++){
            iss.str( var[j].c_str() );
            iss >> value; // convert the string into a double
            value = gsl_matrix_get(Q, j, i) - value; // to have I-Q
            gsl_matrix_set (Q, j, i, value);
            iss.clear();
        }
    }

    // Inversion
    gsl_linalg_LU_decomp (Q, perm, &s); // LU decomposition of matrix Q
    gsl_linalg_LU_invert (Q, perm, N); // invert the matrix N

    // Out put
    for (int i=0; i<diversity; i++){
        TL[i]=0;
        for (int j=0; j<diversity; j++){
            TL[i] += gsl_matrix_get(N, i, j);
        }
    }

    // free the memory
    gsl_matrix_free (Q);
    gsl_matrix_free (N);
    gsl_permutation_free(perm);
}
