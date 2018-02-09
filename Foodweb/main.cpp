#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "string"
#include "math.h"

using namespace std;

int main()
{
string path("data/"); // path to the folder where results are saved
stringstream out; // to convert int into string
string N; // string version of the int

// VARIABLES
int nSim(0); // number written on the file
// modified parameters
int diversity (40); // initial number of species
int nK(91); // number of K values tested
double K[91]; // carrying capacity
for (int i=0; i<nK; i++){
	K[i] = 1+i*0.1;
}
int nX(61); // number of X values tested
double X[61]; // metabolic adaptability rate
for (int i=0; i<nX; i++){
	X[i] = i*0.0001;
}
int replicate(10); // number of replicates

// PARAMETERS
// parameters of the ecosystem
double C(0.15); // initial connectance
int Mmax(3); // maximal body mass
int Mmin(-3); // minimal body mass
// individual parameters
double c(0.5); // predator interference coefficient
double B0(0.5); // half saturation of predation
double FR(1); // Hill exponent for the other consumers
// allometric parameters
double ar(1); // growth rate scaling constant 4.89
double ax(0.314); // metabolic rate scaling constant
double ay(8); // maximum ingestion rate scaling constant
double s(-0.25); // growth rate scaling exponent
double Xmin(0.001); // minimum metabolic rate
double Xmax(1); // minimum metabolic rate
// integration parameters
double t = 0.0, tFinal = 10000; // time span of integration
double tRecord = 9000; // time from witch recording begins
double tStep = 1;
double h = 1e-6; // absolute accuracy
double extinctionThreshold (pow(10,-30)); // extinction biomass threshold

int ndata(15); // number of recorded variables
double data[15]; // recording variables
double species[diversity*2]; // array with the average biomass of species
double CV[diversity*2]; // array with the average coefficient of variation of species
double TL[diversity]; // array with the average trophic level of species
int nbSimu(0); // total number of simulations done

// FILE WITH GENERAL DATA
out << nSim;
N = out.str();

// FILE WITH GENERAL DATA
    string monfichier(path + "data" + N + ".txt");
    ofstream fileGeneral (monfichier.c_str()); // creation of the file of populations
    fileGeneral << "meanPP;";
    fileGeneral << "meanSP;";
    fileGeneral << "cvPP;";
    fileGeneral << "cvSP;";
    fileGeneral << "cvSpecies;";
    fileGeneral << "cvMetabolicRate;";
    fileGeneral << "TLmax;";
    fileGeneral << "connectance;";
    fileGeneral << "connectanceFinal;";
    fileGeneral << "diversity;";
    fileGeneral << "diversityFinal;";
    fileGeneral << "nbPlant;";
    fileGeneral << "simulation;";
    fileGeneral << "X;";
    fileGeneral << "K";
    fileGeneral << endl;

// FILE WITH SPECIES AVERAGE BIOMASSES
monfichier = path + "species" + N + ".txt";
ofstream filespecies (monfichier.c_str()); // creation of the file of populations
filespecies << "nbSimu;X;K;x1" ; // write the header
for (int i=2; i<=diversity; i++){
    filespecies << ";x" << i ;
}
for (int i=1; i<=diversity; i++){
    filespecies << ";m" << i ;
}
filespecies << endl;

// FILE WITH SPECIES AVERAGE COEFFICIENT OF VARIATION
monfichier = path + "CV" + N + ".txt";
ofstream fileCV (monfichier.c_str()); // creation of the file of populations
fileCV << "nbSimu;X;K;x1" ; // write the header
for (int i=2; i<=diversity; i++){
    fileCV << ";x" << i ;
}
for (int i=1; i<=diversity; i++){
    fileCV << ";m" << i ;
}
fileCV << endl;

// FILE WITH SPECIES AVERAGE TROPHIC LEVEL
monfichier = path + "TL" + N + ".txt";
ofstream fileTL (monfichier.c_str()); // creation of the file of populations
fileTL << "nbSimu;X;K;x1" ; // write the header
for (int i=2; i<=diversity; i++){
    fileTL << ";x" << i ;
}
fileTL << endl;


for (int i1=0; i1<nK; i1++){ // K
    for (int i2=0; i2<nX; i2++){ // X
        for (int i3=0; i3<replicate; i3++){ // replicates
            data[7] = C;
            data[9] = diversity;
            data[12] = nSim * nK * nX * replicate + nbSimu;
            data[13] = X[i2];
            data[14] = K[i1];
            Community FoodWeb(
                    // parameters of the ecosystem
                    diversity
                    ,C
                    ,Mmin
                    ,Mmax
                    ,ay
                    ,c
                    ,B0
                    ,K[i1]
                    ,FR
                    ,X[i2]
                    ,Xmin
                    ,Xmax);
            FoodWeb.setMass(); // create the vectors of masses
            FoodWeb.setObject(); // create the vectors of populations and parameters
            int test(0);
            while (test==0){
                FoodWeb.setSpecies(ar, s); // create the species
                FoodWeb.setInteractions(); // build the food web
                test = FoodWeb.checkFoodWeb(); // remove redundant species
            }
            FoodWeb.initialisation(ax, s); // initialisation of the densities
            FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold);
            FoodWeb.matrixTL(path); // create the matrix of interactions
            FoodWeb.output(data, species, CV, TL);

             // FILE WITH GENERAL DATA
            for (int i=0; i<ndata-1; i++){
                fileGeneral << data[i] << ";";
            }
            fileGeneral << data[ndata-1] << endl;

            // Writing in the species file
            filespecies << data[12] << ";"; // number of the simulation
            filespecies << data[13] << ";"; // X
            filespecies << data[14] << ";"; // K
            for (int k=0; k<diversity*2-1; k++){ // species' biomass
                filespecies << species[k] << ";";
            }
            filespecies << species[diversity*2-1] << endl;

            // Writing in the species file
            fileCV << data[12] << ";"; // number of the simulation
            fileCV << data[13] << ";"; // X
            fileCV << data[14] << ";"; // K
            for (int k=0; k<diversity*2-1; k++){ // species' CV
                fileCV << CV[k] << ";";
            }
            fileCV << CV[diversity*2-1] << endl;

            // Writing in the species file
            fileTL << data[12] << ";"; // number of the simulation
            fileTL << data[13] << ";"; // X
            fileTL << data[14] << ";"; // K
            for (int k=0; k<diversity-1; k++){ // species' TL
                fileTL << TL[k] << ";";
            }
            fileTL << TL[diversity-1] << endl;

            nbSimu ++;
            cout << nbSimu << "/" << nK * nX * replicate << endl;
        }
    }
}

fileGeneral.close();
    return 0;
}
