#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <valarray>

using namespace std;

double m_t (double B, double T,  double Jo, double So){

	const double kB = 1.0; 
	const double Beta = 1.0/(kB*T);

	double B_eff = B + So * Jo;

        double x = Beta* B_eff;

        double coth_x = 1.0/tanh( x );

        double m = coth_x - 1.0/x;

        return m;

}



int main(){


    int N = 10;                         //lattice dimension, length of single side
    double B = 0.1;                     //external magnetic field
    
    double * J = new double [2];
    J[0] =2.0; J[1] = 1.5;              //coupling constant
    valarray <double> Jn (J,2);
    
    double T = 0.0;
    double M = 1.0;
  
    std::ostringstream mf;
    mf <<"mft_"<<N<<".txt" ;
    std::string mft = mf.str();

    ofstream file; file.open(mft.c_str());

    //loop over T values
    for (T =0.0001; T <= 15.0; T+=0.05)
    {

      file << T << "\t" << M  << "\n";      
      M  = m_t (B, T, Jn.sum(), M);
    }

 
     
return 0;

}

