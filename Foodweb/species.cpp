#include "species.h"
#include "functions.h"

#include <iostream>
#include "string"
#include "math.h"
#include <cstdlib> // to use calloc
#include <stdio.h>
#include <gsl/gsl_randist.h> // for beta distribution
using namespace std;

#include <sstream>
template<typename T>
string to_string( const T & Value ){
    ostringstream oss; // utiliser un flux de sortie pour créer la chaîne
    oss << Value; // écrire la valeur dans le flux
    return oss.str(); // renvoyer une string
}

// GENERAL FUNCTIONS
Species* Species::getAddress(){
    return (this);
}
double Species::getMass(){
    return(_M);
}
double Species::getMassNorm(){
    return(_Mnorm);
}
void Species::setType(string type){
    _type = type;
    if (type=="Consumer"){
        _e = 0.85;
    }
    else {
        _e = 0.45;
    }
}
string Species::getType(){
    return(_type);
}
void Species::setID(){
    if (_type == "Consumer"){
        _preyID = to_string(_idPrey[0]);
        for (int i=1; i<_nbPrey; i++){
            _preyID += to_string(_idPrey[i]);
        }
        if(_nbPred==0){
            _predID="-1";
        }
        else{
            _predID = to_string(_idPred[0]);
            for (int i=1; i<_nbPred; i++){
                _predID += to_string(_idPred[i]);
            }
        }
        //cout << _type << " | " << _nbPred << " - " << _nbPrey << " <-> " << _predID << " - " << _preyID << endl;
    }
}
string Species::get_preyID(){
    return(_preyID);
}
string Species::get_predID(){
    return(_predID);
}
void Species::setCannib(int cannib){
    _cannib=cannib;
}
int Species::getCannib(){
    return(_cannib);
}

// FUNCTIONS RELATIVE TO PREDATORS
int Species::getNbPred(){
    return(_nbPred);
}
void Species::changeNbPred(int nbPred){
    _nbPred=nbPred;
}
void Species::setPred(){
    if(_nbPred>0){
        _pred = new Species*[_nbPred];
        for (int i=0; i<_nbPred; i++){
            _pred[i] = NULL;
        }
        //_pred = (Species**)calloc(_nbPred, sizeof(Species*)); // create a vector of pointers initialized at 0
        _idPred = new int[_nbPred]; // create the vector of predators' numbers
    }
}
void Species::addPred(int nPred, Species* pred){
    _pred[_countPred]=pred; // pointer toward the predator
    _idPred[_countPred]=nPred; // position of the predator
    _countPred++; // counter of predators
}
Species* Species::getPred(int i){
    return (_pred[i]);
}

// FUNCTIONS RELATIVE TO PREY/RESOURCES
int Species::getNbPrey(){
    return(_nbPrey);
}
void Species::changeNbPrey(int nbPrey){
    _nbPrey=nbPrey;
    _W=nbPrey;
}
void Species::setPrey(){
    _prey = new Species*[_nbPrey];
        for (int i=0; i<_nbPrey; i++){
            _prey[i] = NULL;
        }
    //_prey = (Species**)calloc(_nbPrey, sizeof(Species*)); // create a vector of pointers initialized at 0
    _idPrey = new int[_nbPrey]; // create the vector of preys' numbers
}
void Species::addPrey(int nPrey, Species* prey){
    _prey[_countPrey]=prey; // pointer toward the prey
    _idPrey[_countPrey]=nPrey; // position of the prey
    _countPrey++; // counter of prey
}
int Species::getPreyID(int j){
    return(_idPrey[j]);
}
int Species::getW(){
    return(_W);
}
void Species::changeW(int W){
    _W=W;
}

// constructor
Species::Species(double M
                 ,double Mnorm
                 ,double rk
                 ,int diversity
                 // allometric parameters
                 ,double ar
                 ,double s)
                 :_M(M)
                 ,_Mnorm(Mnorm)
                 ,_rank(rk)
                 ,_type("VOID")
                 ,_diversity(diversity)
                 ,_cannib(0)
                 ,_e(0)
                 ,_nbPred(0)
                 ,_countPred(0)
                 ,_nbPrey(0)
                 ,_countPrey(0)
                 ,_W(0){
                     _R = setAllometric(ar,Mnorm,s);
                 }
// destructor
Species::~Species(){
    if(_nbPred != 0 && _type!="VOID"){
        delete[] _pred;
        //free(_pred);
        delete[] _idPred;
    }
    if(_nbPrey != 0 && _type!="VOID"){
        delete[] _prey;
        //free(_prey);
        delete[] _idPrey;
    }
}
