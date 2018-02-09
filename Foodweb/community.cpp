#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
using namespace std;

#include <cstdlib> // to pick random numbers
#include <ctime>  // to pick random numbers

// INTEGRATION
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_randist.h> // for random distributions

// CREATION OF THE FOODWEB
void Community::setMass(){
    _M = new double[_diversity]; // vector of body mass
    _rank = new double[_diversity]; // vector of rank (position on a [0,1] axis)
    _rn = new double[_diversity]; // vector of niche width
    _cn = new double[_diversity]; // vector of niche centre

    const gsl_rng_type * T; // type of generator
    gsl_rng * r; // generator
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T); // instance of the generator
    gsl_rng_set(r, time(0)); // change the seed of the generator to have different values after each run

    for (int i=0; i<_diversity; i++){
        _rank[i] = gsl_rng_uniform (r); // pick random numbers
        _M[i] = pow(10,_rank[i] * (_Mmax - _Mmin) + _Mmin);
        _rn[i] = 1-pow(1-gsl_rng_uniform (r),2*_C/(1-2*_C)); // beta distribution
        _cn[i] = gsl_rng_uniform (r) * (_rank[i] - _rn[i]*0.5) + _rn[i]*0.5;
        //cout << _rank[i] << " - " << _M[i] << " - " << _cn[i] << " - " << _rn[i] << endl;
    }
    gsl_rng_free(r);
}
void Community::setObject(){
    _params = new Species*[_diversity];
        for (int i=0; i<_diversity; i++){
            _params[i] = NULL;
        }
    _pop = new double[_diversity * 2];
        for (int i=0; i<_diversity * 2; i++){
            _pop[i] = 0;
        }
    _TL = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _TL[i] = 0;
        }
    _tExt = new double[_diversity];
        for (int i=0; i<_diversity; i++){
            _tExt[i] = 0;
        }
    _biomass = new double[_diversity * 2];
        for (int i=0; i<_diversity * 2; i++){
            _biomass[i] = 0;
        }
    _biomassSQ = new double[_diversity * 2];
        for (int i=0; i<_diversity * 2; i++){
            _biomassSQ[i] = 0;
        }
    _CV = new double[_diversity * 2];
        for (int i=0; i<_diversity * 2; i++){
            _CV[i] = 0;
        }
    for (int i=0; i<2; i++){
        _flux[i] = 0;
        _sumFlux[i] = 0;
        _sumFluxSQ[i] = 0;
        _fluxCV[i] = 0;
    }
}
void Community::setSpecies(double ar, double s){
    double Mref(_M[0]);
    for (int i=1; i<_diversity; i++){
        if (_M[i]<Mref){
            Mref = _M[i];
            _iRef = i;
        }
    }
    for (int i=0; i<_diversity; i++){
        _params[i] = new Species (_M[i],_M[i]/Mref,_rank[i],_diversity,ar,s); // creation of primary producers
    }
    _reroll=0;
}
void Community::setInteractions(){
    _nbPlant=0; // reset the number of primary producers
    // COUNT THE NUMBER OF PREY/RESOURCE AND PREDATORS FOR EACH ORGANISM
    for (int i=0; i<_diversity; i++){
        if (i!=_iRef){ // the smallest species cannot be a predator
            for (int j=0; j<_diversity; j++){
                if (_rank[j] > (_cn[i]-_rn[i]*0.5) && _rank[j] < (_cn[i]+_rn[i]*0.5) ){
                    _params[i]->changeNbPrey(_params[i]->getNbPrey()+1); // increment the number of prey/resources
                    _params[j]->changeNbPred(_params[j]->getNbPred()+1); // increment the number of predator
                    if (i==j){
                        _params[i]->setCannib(1); // cannibalism
                    }
                }
            }
        }
    }
    // CREATE THE VECTORS
    for (int i=0; i<_diversity; i++){
        if (_params[i]->getNbPrey()>1){ // consumer with at least one non-cannibalistic prey
            _params[i]->setType("Consumer");
            _params[i]->setPrey();
            _params[i]->setPred();
        }
        else if ( _params[i]->getNbPrey()==1  && _params[i]->getCannib()==0 ){ // predator with one non-cannibalistic prey
            _params[i]->setType("Consumer");
            _params[i]->setPrey();
            _params[i]->setPred();
        }
        else if ( _params[i]->getNbPred()>0 && _params[i]->getNbPrey()==0 ){ // predator with one cannibalistic prey
            _params[i]->setType("Plant");
            _params[i]->changeNbPrey(0);
            _params[i]->setPred();
            _nbPlant++;
        }
        else if ( _params[i]->getNbPred()>1 && _params[i]->getNbPrey()==1  && _params[i]->getCannib()==1 ){ // predator with one cannibalistic prey
            _params[i]->setType("Plant");
            _params[i]->changeNbPrey(0);
            _params[i]->changeNbPred(_params[i]->getNbPred()-1); // remove the cannibalistic predator
            _params[i]->setPred();
            _nbPlant++;
        }
        else { // isolated species
            const gsl_rng_type * T; // type of generator
            gsl_rng * r; // generator
            gsl_rng_env_setup();
            T = gsl_rng_default;
            r = gsl_rng_alloc (T); // instance of the generator
            gsl_rng_set(r, time(0)); // change the seed of the generator to have different values after each run

            _rank[i] = gsl_rng_uniform (r); // pick random numbers
            _M[i] = pow(10,_rank[i] * (_Mmax - _Mmin));
            _rn[i] = 1-pow(1-gsl_rng_uniform (r),2*_C/(1-2*_C)); // beta distribution
            _cn[i] = gsl_rng_uniform (r) * (_rank[i] - _rn[i]*0.5) + _rn[i]*0.5;
            _reroll++;
            gsl_rng_free(r);
        }
    }
    // FILL IN THE VECTORS
    if (_reroll==0){
        for (int i=0; i<_diversity; i++){
            if (_params[i]->getType()=="Consumer"){
                for (int j=0; j<_diversity; j++){
                    if (_rank[j] > (_cn[i]-_rn[i]*0.5) && _rank[j] < (_cn[i]+_rn[i]*0.5) ){
                        _params[i]->addPrey(j,_params[j]); // add the prey
                        _params[j]->addPred(i,_params[i]); // add the predator
                    }
                }
            }
        }
        for (int i=0; i<_diversity; i++){
            if (_params[i]->getType()=="Consumer"){
                _params[i]->setID();
            }
        }
    }
}
int Community::checkFoodWeb(){
    if (_reroll==0){
        for (int i=0; i<_diversity-1; i++){
            if (_params[i]->getType()=="Consumer"){
                for (int j=i+1; j<_diversity; j++){
                    if (_params[j]->getType()=="Consumer" && _params[i]->get_preyID()==_params[j]->get_preyID() && _params[i]->get_predID()==_params[j]->get_predID() ){ // case of trophic identical species

                        const gsl_rng_type * T; // type of generator
                        gsl_rng * r; // generator
                        gsl_rng_env_setup();
                        T = gsl_rng_default;
                        r = gsl_rng_alloc (T); // instance of the generator
                        gsl_rng_set(r, time(0)); // change the seed of the generator to have different values after each run

                        _rank[i] = gsl_rng_uniform (r); // pick random numbers
                        _M[i] = pow(10,_rank[i] * (_Mmax - _Mmin));
                        _rn[i] = 1-pow(1-gsl_rng_uniform (r),2*_C/(1-2*_C)); // beta distribution
                        _cn[i] = gsl_rng_uniform (r) * (_rank[i] - _rn[i]*0.5) + _rn[i]*0.5;
                        _reroll++;
                        gsl_rng_free(r);
                    }
                }
            }
        }
    }
    if (_reroll==0){
        return 1;
    }
    else {
        for (int i=0; i<_diversity; i++){
            delete _params[i];
            _params[i] = NULL;
        }
        return 0;
    }
}
void Community::initialisation(double ax, double s){
    for (int i=0; i<_diversity; i++){
        _pop[i] = 0.1; // species biomass
        if(_params[i]->getType()=="Plant"){
            _pop[i+_diversity] = 0;
        }
        if(_params[i]->getType()=="Consumer"){
            _pop[i+_diversity] = setAllometric(ax,_params[i]->getMassNorm(),s);
        }
    }
}
int Community::Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold){
    //DEFINITIONS FOR GNU LIBRARY
    size_t dimension(_diversity*2);
    //size_t dimension(3);
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // integration method, type of the stepping function
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimension); // creation an instance of the stepping function
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0); // object keeping the local error on each step within an absolute error of 1e-6 and relative error of 0.0
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (dimension); // instance of an evolution function for a system of 1 dimensions
    gsl_odeiv_system sys = {func, NULL, dimension, this}; // the Jacobian is useless with this method : use NULL instead
    // INTEGRATIONS PARAMETERS
//    double t = 0.0, t1 = 3; // time span of integration
//    double tRecord = 2; // time from witch recording begins
//    double h = 1e-6; // absolute accuracy

/////////////////
// INTEGRATION //
/////////////////

    double y[_diversity*2]; // array with the biomass
    for (int i=0; i<_diversity*2; i++){
        y[i]=_pop[i]; // initialisation
    }
    _nPoint=0; // counter of the number of recorded points
    double tPoint(tRecord); // points when data are recorded


    while (t < tFinal){
        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tFinal, &h, y); // integration
        if (status != GSL_SUCCESS)
        break;

        // Extinction
        for (int i=0; i<_diversity; i++){
            if (y[i]<extinctionThreshold && y[i]!=0){
                y[i]=0; // extinction
                y[i+_diversity]=0; // metabolic rate drops to zero
                _TL[i]=0;
                _tExt[i]=t;
                if(_params[i]->getNbPred()>0){
                    for (int j=0; j<_params[i]->getNbPred(); j++){
                        _params[i]->getPred(j)->changeW(_params[i]->getPred(j)->getW()-1); // actualize the number of preys of the predators
                    }
                }
            }
        }
        // Metabolic rate limitation
        for (int i=_diversity; i<_diversity*2; i++){
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
            for (int i=0; i<_diversity*2; i++){
                _biomass[i]+=y[i]; // add the value of the biomass and metabolic rates
                _biomassSQ[i]+=pow(y[i],2); // sum of square of biomass and metabolic rates
            }

            // PRODUCTION
            _sumFlux[0] += _flux[0];
            _sumFlux[1] += _flux[1];
            _sumFluxSQ[0] += pow(_flux[0],2);
            _sumFluxSQ[1] += pow(_flux[1],2);

            tPoint += tStep ;
            _nPoint++;
        }
    }

    // END !!!
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    return 0;
}
// OUTPUT
void Community::matrixTL(string path){
    string monfichier(path + "matrixInteraction.txt");
    ofstream fileMatrix (monfichier.c_str()); // creation of the file containing the matrix of interactions

    double matrix[_diversity][_diversity];
    for (int i=0; i<_diversity; i++){
        if (_params[i]->getType()=="Consumer"){
            for (int j=0; j<_diversity; j++){
                matrix[j][i]=0;
                for (int k=0; k<_params[i]->getNbPrey(); k++){
                    matrix[_params[i]->getPreyID(k)][i] = (double)1/_params[i]->getW();
                }
            }
        }
        else {
            for (int j=0; j<_diversity; j++){
                matrix[j][i]=0;
            }
        }
    }

    for (int i=0; i<_diversity; i++){
        for (int j=0; j<_diversity-1; j++){
            fileMatrix << matrix[i][j] << ";" ;
        }
        fileMatrix << matrix[i][_diversity];
        fileMatrix << endl;
    }

    trophicLevel(_diversity, monfichier, _TL); // calculation of trophic levels
}
void Community::output(double data[], double species[], double CV[], double TL[]){

    // TROPHIC LEVEL
    for (int i=0; i<_diversity; i++){
        TL[i] = _TL[i];
        if (_TL[i]>_TLmax){
            _TLmax = _TL[i]; // keep the highest TL
        }
    }

    // BIOMASS
    for (int i=0; i<_diversity*2; i++){
        _biomass[i] /= _nPoint;
        species[i] = _biomass[i]; // RECORDING OF THE AVERAGE BIOMASS
    }

    // REMAINING SPECIES
    for (int i=0; i<_diversity; i++){
        if (_tExt[i]==0){
            _diversityFinal++;
        }
    }

    // REMAINING PRIMARY PRODUCERS
    _nbPlant = 0;
    for (int i=0; i<_diversity; i++){
        if (_tExt[i]==0 && _params[i]->getType()=="Plant"){
            _nbPlant++;
        }
    }

    // CONNECTIVITY & CONNECTANCE
    int S(0); // total number of interactions
    for (int i=0; i<_diversity; i++){
        S += _params[i]->getW();
    }
    _connectance = S/pow(_diversityFinal,2);

    // MEAN FLUXES AND COEFFICIENT OF VARIATION
    for (int i=0; i<2; i++){
        _sumFlux[i] /= _nPoint; // mean
        if (_sumFlux[i] != 0){
            _fluxCV[i] = pow( _sumFluxSQ[i]/_nPoint-pow(_sumFlux[i],2) ,0.5) / _sumFlux[i]; // coefficient of variation
        }
    }
    for (int i=0; i<_diversity; i++){
        if (_tExt[i] == 0){
            // species biomass
            _biomass[i] /= _nPoint; // mean
            _CV[i] = pow( _biomassSQ[i]/_nPoint-pow(_biomass[i],2) ,0.5) / _biomass[i]; // coefficient of variation
            CV[i] = _CV[i]; // record in the output vector
            // metabolic rate
            _biomass[i+_diversity] /= _nPoint; // mean
            _CV[i+_diversity] = pow( _biomassSQ[i+_diversity]/_nPoint-pow(_biomass[i+_diversity],2) ,0.5) / _biomass[i+_diversity]; // coefficient of variation
            CV[i+_diversity] = _CV[i+_diversity]; // record in the output vector
        }
    }
    for (int i=0; i<_diversity; i++){
        _meanCV += _CV[i];
        _meanCVx += _CV[i+_diversity];
    }
    _meanCV /= _diversityFinal;
    _meanCVx /= (_diversityFinal - _nbPlant); // plants have no adaptive metabolic rate

    // WRITING
    data[0] = _sumFlux[0]; // average primary production
    data[1] = _sumFlux[1]; // average secondary production
    data[2] = _fluxCV[0]; // primary production coefficient of variation
    data[3] = _fluxCV[1]; // secondary production coefficient of variation
    data[4] = _meanCV; // average coefficient of variation
    data[5] = _meanCV; // average coefficient of variation
    data[6] = _TLmax; // average trophic level
    data[8] = _connectance; // final connectance
    data[10] = _diversityFinal; // final species diversity
    data[11] = _nbPlant; // final number of primary producers

}

// constructor

Community* Community::getAdress(){
    return (this);
}

Community::Community(
                    // parameters of the ecosystem
                    int diversity,
                    double C,
                    int Mmin,
                    int Mmax,
                    double ay,
                    double c,
                    double B0,
                    double K,
                    double FR,
                    double X,
                    double Xmin,
                    double Xmax)
                     :
                    // parameters of the ecosystem
                    _diversity(diversity)
                    ,_nbPlant(0)
                    ,_C(C)
                    ,_Mmin(Mmin)
                    ,_Mmax(Mmax)
                    ,_y(ay)
                    ,_c(c)
                    ,_B0(B0)
                    ,_K(K)
                    ,_FR(FR)
                    ,_X(X)
                    ,_Xmin(Xmin)
                    ,_Xmax(Xmax)
                    ,_reroll(0)
                    ,_iRef(0)
                    ,_meanCV(0)
                    ,_meanCVx(0)
                    ,_connectance(0)
                    ,_TLmax(0)
                    ,_diversityFinal(0){
                        for (int i=0; i<2; i++){
                            _flux[i]=0;
                            _sumFlux[i]=0;
                            _sumFluxSQ[i]=0;
                            _fluxCV[i]=0;
                        }
                    }

// destructor
Community::~Community(){
    delete[] _rank;
    delete[] _rn;
    delete[] _cn;
    delete[] _M;
    for (int i=0; i<_diversity; i++){
        delete _params[i];
    }
    delete[] _params;
    //free(_params);
    delete[] _pop;
    delete[] _TL;
    delete[] _tExt;
    delete[] _biomass;
    delete[] _biomassSQ;
    delete[] _CV;
}
