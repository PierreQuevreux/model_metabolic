#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "string"
#include <gsl/gsl_matrix.h>

using namespace std;

double setAllometric(double a, double M, double s);

int func(double t, const double y[], double f[], void *params); // system of ODE

#endif // FUNCTIONS_H_INCLUDED
