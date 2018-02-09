#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "string"
using namespace std;

double setAllometric(double a, double M, double s);

int func(double t, const double y[], double f[], void *params); // system of ODE

void trophicLevel(int diversity, string fileMatrixInteraction, double TL[]); // calculate the trophic levels from the interaction matrix

#endif // FUNCTIONS_H_INCLUDED
