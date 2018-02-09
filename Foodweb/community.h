#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "species.h"
#include "string"
using namespace std;

class Community
{
private :
    // PARAMETERS
    // parameters of the ecosystem
    int _diversity; // initial number of species
    int _nbPlant; // number of primary producers
    double _C; // initial connectance
    int _Mmin; // minimal body mass
    int _Mmax; // maximal body mass
    double _y; // maximum ingestion rate coefficient
    double _c; // predator interference coefficient
    double _B0; // half saturation constant
    double _K; // primary producers carrying capacity
    double _FR; // Hill exponent
    double _X; // metabolic adaptability rate
    double _Xmin; // minimum metabolic rate
    double _Xmax; // minimum metabolic rate

    int _reroll; // counter of species to rebuild
    int _iRef; // position of the reference species

    //STRUCTURES
    double *_M; // pointer toward the vector of body mass
    double *_rank; // pointer toward the linear size distribution of species  (position on a [0,1] axis)
    double *_rn; // pointer toward the vector of niche width
    double *_cn; // pointer toward the vector of niche centre

    //INTEGRATION STRUCTURES
    Species **_params; // pointer toward the vector of pointers towards the instanced classes
    double *_pop; // pointer toward the vector of biomasses
    double *_TL; // pointer toward the vector of trophic levels
    double *_tExt; // pointer toward the vector of extinction time of species
    int _nPoint; // number of recorded points counter

    //RECORDED OUT PUT
    double *_biomass; // sum of biomass for each species
    double *_biomassSQ; // sum of square of biomass for each species
    double *_CV; // coefficient of variation of the biomass of species
    double _meanCV; // mean CV across species
    double _meanCVx; // mean CV across species' metabolic rate
    double _flux[2]; // array with the value of PP and SP at each step
    double _sumFlux[2]; // array with the sum of the value of PP and SP
    double _sumFluxSQ[2]; // array with the sum of the square value of PP and SP
    double _fluxCV[2]; // array with the CV of PP and SP  at each step
    double _connectance; // final connectance
    double _TLmax; // maximal TL
    int _diversityFinal; // number of species at the end

public :
    // CREATION OF THE FOODWEB
    void setMass();
    void setObject(); // create the vectors of populations and parameters
    void setSpecies(double ar, double s); // create the species
    void setInteractions();
    int checkFoodWeb();
    void initialisation(double ax, double s); // initialize the vector of populations for the resources
    //INTEGRATION
    int Dynamic(double t, double tFinal, double tRecord, double tStep, double h, double extinctionThreshold);
    void matrixTL(string path); // build the interaction matrix to calculate the trophic levels
    void output(double data[], double species[], double CV[], double TL[]); // write the results in the recording structures

    Community* getAdress();

    //constructor
    Community(
        // parameters of the ecosystem
        int diversity,
        double C,
        int Mmin,
        int Mmax,
        double _y,
        double _c,
        double _B0,
        double _K,
        double _FR,
        double _X,
        double _Xmin,
        double _Xmax);
    // destructor
    ~Community();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE

};

#endif // COMMUNITY_H_INCLUDED
