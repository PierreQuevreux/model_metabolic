#include "community.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream> // to use istringstream and convert string to double
#include "string"
#include "math.h"

using namespace std;

int main()
{
string path("data/");

int dim(6); // dimension of the system
int Dim(5); // dimension of the linearised system
// modified parameters
double X(0); // metabolic adaptability rate //2
double Km(9); // maximum K for the bifurcation
double Kstep(0.1); // step of the exploration of K
int n(static_cast<int>(Km/Kstep)+1); // number of steps
double K[n]; // carrying capacity
for(int i=0; i<n; i++){
    K[i] = 1 + i * Kstep;
}
// individual parameters
double c(0.5); // predator interference coefficient
double B0(0.5); // half saturation of predation
// allometric parameters
double ar(1); // growth rate scaling constant 4.89
double ax(0.314); // metabolic rate scaling constant
double ay(8); // maximum ingestion rate scaling constant
double s(-0.25); // scaling exponent
double FR(1); // Hill exponent
double Xmin(0.001); // minimum metabolic rate
double Xmax(1); // minimum metabolic rate
// integration parameters
double t = 0.0, tFinal = 10000; // time span of integration
double tRecord = 9000; // time from witch recording begins
double tStep = 0.01;
double h = 1e-6; // absolute accuracy
double extinctionThreshold (pow(10,-30)); // extinction biomass threshold

/////////////////////
// RECORDING FILES //
/////////////////////

// Lyapunov exponent calculation
string monfichier(path + "lyapunov.txt");
ofstream fileLyapunov (monfichier.c_str()); // creation of the file of populations
fileLyapunov << "K;X;LE1;LE2;LE3;LE4;LE5" << endl; // write the header

/////////////////
// SIMULATIONS //
/////////////////
for(int i=0; i<n; i++){
    Community FoodWeb(
            dim
            ,Dim
            ,ar
            ,ay
            ,c
            ,B0
            ,FR
            ,K[i]
            ,X
            ,Xmin
            ,Xmax);
    FoodWeb.setSystem(ax,s); // create the vectors of masses
    FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold,fileLyapunov);
}

    return 0;
}
