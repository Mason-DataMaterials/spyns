#include <iomanip>
#include <random>
#include <map>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;


double z_N (int N, double beta, double J){

	double Z = 0.0;
	Z = (4.0*M_PI / beta*J) * sinh(beta*J);
        double ZN = pow(Z, (double)N);
        cout << sinh(beta*J) << "\t" << Z << "\t" << ZN << "\n";

	return ZN;
}



int main (int argc, char *argv[])
{
    int N;
    double * J = new double [2];
    
    J[0] = 13.1; //meV
    J[1] = 13.7; //meV
    
    double JSum;
    
    double temperature = 100;
    
    double kB = 8.617330350E-5; //eVK^-1
    double kT = 1000*kB*temperature; //meV 
    
    int num_1st_nn,num_2nd_nn;
    int * num_nn;
    num_nn = new int [2];
    
    
    double ** first_nn;
    double ** second_nn;
    
    ifstream infile2 ("neighbor_lists.txt");
    if (infile2.is_open())
    {
        
        infile2>>N;
        infile2>>num_1st_nn>>num_2nd_nn;
        
        num_nn[0] = num_1st_nn;
        num_nn[1] = num_2nd_nn;
        
        //populate neighbor lists
        first_nn = new double * [N];
        second_nn = new double * [N];
        
        for (int i = 0; i < N; i++)
        {
            first_nn[i] = new double [num_1st_nn];
            second_nn[i] = new double [num_2nd_nn];
        }  
        
        for (int i = 0; i < N; i++){
            for (int j = 0; j < num_1st_nn; j++){
                infile2 >> first_nn[i][j]; 
            }
            
            for (int j = 0; j < num_2nd_nn; j++){
                infile2 >> second_nn[i][j]; 
            }
        }
    }
        
        infile2.close(); 
        
        
        JSum = J[0]*num_nn[0] + J[1] * num_nn[1];
        
        double Beta = 1.0/kT;
        double BJ = Beta*JSum;
        
        double const exp_dist_mean   = BJ ;
        double const exp_dist_lambda = 1.0 / exp_dist_mean;
        
        
        std::random_device rd; 
        std::exponential_distribution<> rng (exp_dist_lambda);
        std::mt19937 rnd_gen (rd ());
        
        
        
        double z = z_N( N, Beta, JSum);
        z = 1.0/z;
        
        double cos_theta_list[N];
        for (int i =0; i < N; ++i){
            cos_theta_list[i] = rng (rnd_gen) ; 
            
            
          cout << i << "\t" << z << "\t" << cos_theta_list[i]/Beta*JSum << "\n";
        }
        
        
        
        return 0;
    }
    
