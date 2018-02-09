#ifndef SPECIES_H_INCLUDED
#define SPECIES_H_INCLUDED

#include "string"
using namespace std;

class Species
{
protected :
    // GENERAL PARAMETERS
    double _M; // body mass
    double _Mnorm; // normalised body mass
    double _rank; // rank of the species
    string _type; // type of organism
    int _diversity; // number of species
    int _cannib; // cannibal ?

    // FUNCTIONAL PARAMETERS
    double _R; // growth rate
    double _x; // metabolic rate
    double _e; // assimilation efficiency by predators
    double _r; // width of the feeding niche
    double _c; // centre of the feeding niche

    // PREDATORS
    int _nbPred; // number of predators
    int _countPred; // counter to fill in the array of predators' pointers
    Species** _pred; // array containing the pointers toward the predators
    int* _idPred; // array containing the number of the predators

    // PREY/RESOURCES
    int _nbPrey; // initial number of prey/resource
    int _countPrey; // counter to fill in the array of prey' pointers
    int _W; // current number of preys
    Species** _prey; // array containing the pointers toward the prey/resource
    int* _idPrey; // array containing the number of the prey/resource

    string _preyID;
    string _predID;

public :
    // GENERAL FUNCTIONS
    Species* getAddress(); // return the pointer this
    double getMass(); // function returning the body mass
    double getMassNorm(); // function returning the normalised body mass
    void setType(string type); // function entering the type
    string getType(); // function returning the type
    void setID(); // convert the vectors _idPrey and _idPred into words what can be compared
    string get_preyID(); // function returning prey word
    string get_predID(); // function returning prey word
    void setCannib(int cannib); // establish the cannibalism
    int getCannib(); // function returning the cannibalism

    // FUNCTIONS RELATIVE TO PREDATORS
    int getNbPred(); // function returning the number of predators
    void changeNbPred(int nbPred); // change the total number of predators
    void setPred(); // create the vector of predators'pointers
    void addPred(int nPred,Species* pred); // add the pointer of the predator to the predators list
    Species* getPred(int i); // function returning the pointer toward the predator i

    // FUNCTIONS RELATIVE TO PREY/RESOURCES
    int getNbPrey(); // function returning the number of preys
    void changeNbPrey(int nbPrey); // change the total number of preys
    void setPrey(); // create the vector of predators'pointers
    void addPrey(int nPrey,Species* prey); // add the pointer of the predator to the predators list
    int getPreyID(int j); // function returning the position of the prey
    int getW(); // function returning _W
    void changeW(int W); // actualize the value of _W

    //constructor
    Species(double M
            ,double Mnorm
            ,double rk
            ,int diversity
            // allometric parameters
            ,double ar
            ,double s);
    // destructor
    ~Species();

    //INTEGRATION
    friend int func(double t, const double y[], double f[], void *params); // system of ODE
};

#endif // SPECIES_H_INCLUDED
