#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "string"
#include <gsl/gsl_matrix.h>
using namespace std;

class Community
{
private :
    // PARAMETERS
    int _dim; // dimension of the system
    int _Dim; // dimension of the linearised system
    bool _record; // start of the recording
    // parameters of the ecosystem
    double _r; // growth rate of the primary producer
    double _x[3]; // initial metabolic rate
    double _y; // maximum ingestion rate
    double _M[3]; // body masses of species
    double _e[3]; // assimilation efficiency
    double _c; // predator interference coefficient
    double _B0; // half saturation constant
    double _FR; // Hill exponent
    double _K; // primary producers carrying capacity
    double _X; // metabolic adaptability rate
    double _Xmin; // minimum metabolic rate
    double _Xmax; // minimum metabolic rate
    double _pop[3]; // initial biomasses
    // terms used to calculate the Jacobian matrix
    double _element; // Element of the matrix
    double _denomB2; // Denominator of the functional response of the herbivore
    double _denomB3; // Denominator of the functional response of the carnivore
    double _term1; // loss B1B1
    double _term2; // loss B1B2
    double _term3; // loss B2B2
    double _term4; // loss B2B3
    double _term5; // loss B1x2
    double _term6; // loss B2x3
    // Lyapunov exponents
    gsl_matrix * _jac; // Jacobian
    gsl_matrix * _vect; // matrix of vectors
    gsl_matrix * _prod; // results of the jac*vect product

public :
    // CREATION OF THE FOODWEB
    void setSystem(double ax, double s);
    //INTEGRATION
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, ofstream &fileLyapunov);
    // Return the Jacobian
    void jacobian(const double y[], double f[]);

    //constructor
    Community(
        int dim,
        int Dim,
        double r,
        double y,
        double c,
        double B0,
        double FR,
        double K,
        double X,
        double Xmin,
        double Xmax);
    // destructor
    ~Community();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE

};

#endif // COMMUNITY_H_INCLUDED
