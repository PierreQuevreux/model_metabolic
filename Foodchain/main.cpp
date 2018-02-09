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

// modified parameters
double X(2); // metabolic adaptability rate
double Km(9); // maximum K for the bifurcation
double Kstep(0.1); // step of the exploration of K
int n(static_cast<int>(Km/Kstep)+1); // number of steps
double K[n]; // carrying capacity
int nPoint[n]; // number of simulations
for(int i=0; i<n; i++){
    K[i] = 1 + i * Kstep;
    nPoint[i] = 0;
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
double t = 0.0, tFinal = 2000; // time span of integration
double tRecord = 1000; // time from witch recording begins
double tStep = 1;
double h = 1e-6; // absolute accuracy
double extinctionThreshold (pow(10,-30)); // extinction biomass threshold

/////////////////////
// RECORDING FILES //
/////////////////////

string monfichier(path + "chronic.txt");
ofstream fileChronic (monfichier.c_str()); // creation of the file of populations
fileChronic << "t;K;X;PP;SP" ; // write the header
for (int i=0; i<3; i++){
    fileChronic << ";Species" << i+1 ; // names of resources and species
}
for (int i=0; i<3; i++){
    fileChronic << ";x" << i+1 ; // names of resources and species
}
fileChronic << endl;

/////////////////
// SIMULATIONS //
/////////////////
for(int i=0; i<n; i++){
    Community FoodWeb(
            ar
            ,ay
            ,c
            ,B0
            ,FR
            ,K[i]
            ,X
            ,Xmin
            ,Xmax);
    FoodWeb.setSystem(ax,s); // create the vectors of masses
    nPoint[i] = FoodWeb.Dynamic(t,tFinal,tRecord,tStep,h,extinctionThreshold,fileChronic);
}
fileChronic.close();
/////////////////
// BIFURCATION //
/////////////////
int N(0); // total number of recorded points
for(int i=0; i<n; i++){
    N += nPoint[i];
}
cout << N << endl;
double** chronic = new double*[N]; // table with the data (lines)
for (int i=0; i<N; i++){
    chronic[i] = new double[8]; // columns
}

// OPENING OF THE FILE AND STORAGE OF THE INFORMATION IN THE TABLE
ifstream file(monfichier.c_str());
if (file){}
else {cout << "ERROR: unable to open file at path" << endl;}
string line;
string car;
getline(file, line); // skip the header
int k(0); // counter for the table of names
string var[11]; // string array containing the string form of the parameters
istringstream iss; // variable for the conversion from string to double
for(int i=0; i<N; i++){
    getline(file, line);
    k=0;
    for (int j=0; j<11; j++){
        var[j].clear(); // clear the string array
    }
    for (unsigned int j=0; j<line.size(); j++){
        car = line[j];
        if (car!=";"){
            var[k] += line[j]; // write the data letter by letter
        }
        else {k++;}
    }
    iss.str( var[1].c_str() ); // K
    iss >> chronic[i][0];
    iss.clear();
    iss.str( var[2].c_str() ); // X
    iss >> chronic[i][1];
    iss.clear();
    iss.str( var[5].c_str() ); // Species1
    iss >> chronic[i][2];
    iss.clear();
    iss.str( var[6].c_str() ); // Species2
    iss >> chronic[i][3];
    iss.clear();
    iss.str( var[7].c_str() ); // Species3
    iss >> chronic[i][4];
    iss.clear();
    iss.str( var[8].c_str() ); // Species1 metabolic rate
    iss >> chronic[i][5];
    iss.clear();
    iss.str( var[9].c_str() ); // Species2 metabolic rate
    iss >> chronic[i][6];
    iss.clear();
    iss.str( var[10].c_str() ); // Species3 metabolic rate
    iss >> chronic[i][7];
    iss.clear();
}

// FINDING OF EXTREMA AND WRITING IN THE BIFURCATION FILE
monfichier = path + "bifurcation.txt";
ofstream fileBifurcation (monfichier.c_str()); // creation of the file of populations
fileBifurcation << "K;X;value;species" << endl ; // write the header
k=0;
for(int i=1; i<nPoint[0]-1; i++){
    if( (chronic[i][2]<=chronic[i-1][2] && chronic[i][2]<=chronic[i+1][2]) || (chronic[i][2]>=chronic[i-1][2] && chronic[i][2]>=chronic[i+1][2]) ){
        fileBifurcation << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][2] << ";" << "Species1" << endl;
    }
    if( (chronic[i][3]<=chronic[i-1][3] && chronic[i][3]<=chronic[i+1][3]) || (chronic[i][3]>=chronic[i-1][3] && chronic[i][3]>=chronic[i+1][3]) ){
        fileBifurcation << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][3] << ";" << "Species2" << endl;
    }
    if( (chronic[i][4]<=chronic[i-1][4] && chronic[i][4]<=chronic[i+1][4]) || (chronic[i][4]>=chronic[i-1][4] && chronic[i][4]>=chronic[i+1][4]) ){
        fileBifurcation << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][4] << ";" << "Species3" << endl;
    }
}
for(int j=1; j<n; j++){
    k += nPoint[j-1];
    for(int i=k+1; i<k+nPoint[j]-1; i++){
        if( (chronic[i][2]<=chronic[i-1][2] && chronic[i][2]<=chronic[i+1][2]) || (chronic[i][2]>=chronic[i-1][2] && chronic[i][2]>=chronic[i+1][2]) ){
            fileBifurcation << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][2] << ";" << "Species1" << endl;
        }
        if( (chronic[i][3]<=chronic[i-1][3] && chronic[i][3]<=chronic[i+1][3]) || (chronic[i][3]>=chronic[i-1][3] && chronic[i][3]>=chronic[i+1][3]) ){
            fileBifurcation << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][3] << ";" << "Species2" << endl;
        }
        if( (chronic[i][4]<=chronic[i-1][4] && chronic[i][4]<=chronic[i+1][4]) || (chronic[i][4]>=chronic[i-1][4] && chronic[i][4]>=chronic[i+1][4]) ){
            fileBifurcation << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][4] << ";" << "Species3" << endl;
        }
    }
}

monfichier = path + "bifurcationMR.txt";
ofstream fileBifurcationMR (monfichier.c_str()); // creation of the file of populations
fileBifurcationMR << "K;X;value;species" << endl ; // write the header
k=0;
for(int i=1; i<nPoint[0]-1; i++){
    if( (chronic[i][5]<=chronic[i-1][5] && chronic[i][5]<=chronic[i+1][5]) || (chronic[i][5]>=chronic[i-1][5] && chronic[i][5]>=chronic[i+1][5]) ){
        fileBifurcationMR << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][5] << ";" << "Species1" << endl;
    }
    if( (chronic[i][6]<=chronic[i-1][6] && chronic[i][6]<=chronic[i+1][6]) || (chronic[i][6]>=chronic[i-1][6] && chronic[i][6]>=chronic[i+1][6]) ){
        fileBifurcationMR << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][6] << ";" << "Species2" << endl;
    }
    if( (chronic[i][7]<=chronic[i-1][7] && chronic[i][7]<=chronic[i+1][7]) || (chronic[i][7]>=chronic[i-1][7] && chronic[i][7]>=chronic[i+1][7]) ){
        fileBifurcationMR << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][7] << ";" << "Species3" << endl;
    }
}
for(int j=1; j<n; j++){
    k += nPoint[j-1];
    for(int i=k+1; i<k+nPoint[j]-1; i++){
        if( (chronic[i][5]<=chronic[i-1][5] && chronic[i][5]<=chronic[i+1][5]) || (chronic[i][5]>=chronic[i-1][5] && chronic[i][5]>=chronic[i+1][5]) ){
            fileBifurcationMR << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][5] << ";" << "Species1" << endl;
        }
        if( (chronic[i][6]<=chronic[i-1][6] && chronic[i][6]<=chronic[i+1][6]) || (chronic[i][6]>=chronic[i-1][6] && chronic[i][6]>=chronic[i+1][6]) ){
            fileBifurcationMR << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][6] << ";" << "Species2" << endl;
        }
        if( (chronic[i][7]<=chronic[i-1][7] && chronic[i][7]<=chronic[i+1][7]) || (chronic[i][7]>=chronic[i-1][7] && chronic[i][7]>=chronic[i+1][7]) ){
            fileBifurcationMR << chronic[i][0] << ";" << chronic[i][1] << ";" << chronic[i][7] << ";" << "Species3" << endl;
        }
    }
}

for (int i=0; i<N; i++){
    delete chronic[i];
}
delete[] chronic;

    return 0;
}
