#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "string"
#include <gsl/gsl_matrix.h>
using namespace std;

class Community
{
private :
    // PARAMETERS
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
    double _flux[2]; // array with the value of PP and SP at each step

public :
    // CREATION OF THE FOODWEB
    void setSystem(double ax, double s);
    //INTEGRATION
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold, ofstream &fileChronicofstream);

    //constructor
    Community(
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
