//
//  main.cpp
//  IsingModel
//
//  Created by Swabir Silayi on 12/2/17.
//  Copyright © 2017 Swabir Silayi. All rights reserved.
//

//#include "ising_funcs.hpp"

#include <stdio.h>
#include <sstream>
#include <string>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


#include <omp.h>
#include <fftw3.h>


#define REAL 0
#define IMAG 1


using namespace std;

/* Random number generator
 * using the Mersenne Twister
 * http://stackoverflow.com/a/22935167
 */
random_device rd;
mt19937 engine(rd());
uniform_real_distribution<double> dist(0, 1);



void initialize(int N, double& L, int key, double ** si, double a){
    
    for (int i = 0; i < N; i++)
        si[i] = new double [4];
    
    
    double u = 1.0;
    double d = -1.0;
    
    
    // find M large enough to fit N atoms bcc lattice
    int M = 1;
    while (2 * M * M * M < N)
        ++M;
    L = a*M;           // lattice constant of conventional cell
    
    // 2 atomic positions in bcc unit cell
    double dx[2] = {0.0, 0.5};
    double dy[2] = {0.0, 0.5};
    double dz[2] = {0.0, 0.5};
    
    int n = 0;
    
    for (int x = 0; x < M; x++){
        for (int y = 0; y < M; y++){
            for (int z = 0; z < M; z++)
                for (int p = 0; p < 2; p++){
                    if (n < N){
                        si[n][0] = (x + dx[p]) * a;
                        si[n][1] = (y + dy[p]) * a;
                        si[n][2] = (z + dz[p]) * a;
                        
                        double p  = dist(engine);
                        if (p > 0.5) si[n][3] = u;
                        else si[n][3] = d;
                        /*
                        if (key == 0)//nm
                        {
                            si[n][3] = 0.0;
                        }
                        
                        
                        if (key == 1)//fm
                        {
                            si[n][3] = u;
                        }
                        
                        
                        if (key == 2)//afm_1
                        {
                            if (n%2 == 0 )
                                si[n][3] = u;
                            else if (n%2 == 1 )
                                si[n][3] = d;
                        }
                        
                        if (key == 3)//afm_2
                        {
                            if(n==0) si[n][3] = u;
                            else if(n==1) si[n][3] = d;
                            
                            else if (n%2 == 0 && n != 0)
                                si[n][3] = -si[n-2][3];
                            else if (n%2 == 1 && n > 1)
                                si[n][3] = -si[n-2][3];
                        }
                        */
                        ++n;
                    }
                }
        }
    }
    return;
}

void neighbors ( int N, double L, double a, double * point, double ** si, double * sj)
{
    
    double nnbcc[14] = {0.0, sqrt(3.0)/2.0, 1.0, sqrt(2.0), sqrt(11.0/4.0),
        sqrt(3.0), 2.0, sqrt(19.0/4.0), sqrt(5.0), sqrt(6.0),
        sqrt(27.0/4.0), sqrt(8.0), sqrt(35.0/4.0), 3.0};
    
    sj[0] = sj[1] = 0.0;
    
    for (int i = 0; i < N; i++){
        
        
        double rSqd =0.0;
        double dr[3] = {0.0};
        double dist =0.0;
        
        dr[0] = -si[i][0] + point[0];
        dr[1] = -si[i][1] + point[1];
        dr[2] = -si[i][2] + point[2];
        
        
        //min-image convention
        for (int d = 0; d < 3; d++)
        {
            if (dr[d] >= 0.5*L)
                dr[d] = dr[d]- L;
            if (dr[d] < -0.5*L)
                dr[d] = dr[d]+L;
            
        }
        
        rSqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        dist = sqrt(rSqd);
        
        //cout << nnbcc[1]*a << "\t" << dist << "\n";
        if ( round(dist-nnbcc[1]*a) == 0.0){
            sj[0]+=si[i][3];
        }
        
        if ( round(dist-nnbcc[2]*a) == 0.0){
            sj[1]+=si[i][3];
            
        }
        
    }
    
    return;
    
}
/* Magnetization(N, si)
 * sums the spins of each site and returns total
 * spin/magnetization
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 * returns: M - sum of all spins
 */
double Magnetization(int N, double ** si){
    
    double M = 0.0;
    
    for (int i = 0; i < N; i++)
    {
                M = M+si[i][3];
    }
    
    return M;
    
}
/* Energy(N, si, h, J)
 * calculates total energy of the grid
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 * returns: H - total energy of grid
 */
double Energy( int N, double L, double a, double ** si, double h, double * J){
    
    double H = 0.0;
    
    double * sj = new double [2];
    sj[0] = sj[1] = 0.0;
    
    for (int i = 0; i < N; i++)
    {
                neighbors ( N, L,a, si[i], si, sj);
                H = H + ( -(si[i][3]*(J[0]*sj[0] + J[1]*sj[1])) - h*si[i][3] );
    }
    

    free (sj);
    
    return H / 2.0;
}

/* deltaE(N, si, h, J, point)
 * calculates the change in energy if a spin site is flipped
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 *           point - site of interest
 * returns: dE - change in energy from flipping point
 */
double deltaE (int N, double L, double a, double ** si, double h, double * J, int i)
{
    double * sj = new double [2];
    sj[0] = sj[1] = 0.0;
    
    neighbors ( N, L,a, si[i], si, sj);
    double dE = 2.0 * ( (si[i][3]*(J[0]*sj[0] +J[1]*sj[1]))
                       + (h*si[i][3]) );
    
    
    free(sj);
    
    return dE;
}

/* MC_move(N, si, h, J, M, H, T)
 * single step Monte Carlo move to determine
 * spin flip of a random site and update the Magnetization, M,
 * and energy, H.
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 *           M - total magnetization
 *           H - total energy
 *           T - temperature
 */
void MC_move(int N, double L, double a, double  ** si, double h, double * J,
             double &M, double &H, double T)
{
    
    double beta = 1.0/T;
    int * point = new int [3];
    
    for (int sweep = 0; sweep <N; sweep++)
    {
        //random spin site
        uniform_int_distribution<int> nsites_gen(0, N - 1);
        int i = nsites_gen(engine);
        
        
        double dE = deltaE(N, L, a, si, h, J, i);
        double dm = 0.0, dh = 0.0;
        
        double p = dist(engine);
        
        if (dE < 0.0)
        {
            
            si[i][3] *= -1;
            dm = 2.0*si[i][3];
            dh = dE;
        }
        else if (p < exp(-dE*beta))
        {
            si[i][3] *= -1;
            dm = 2.0*si[i][3];
            dh = dE;
        }
        
        M = M+dm;
        H = H+dh;
        
    }
    
    
    free (point);
    return;
}



/* ave_var(cc, Avg, Var, U)
 * calculates and upates the instantaneous average (Avg)
 * and variance (Var) values of a quantity (U) at a time
 * step (cc).
 * takes in: cc  - time step
 *           U   - quantity being averaged (M,H, M^2, H^2)
 *           Avg - variable for the average value
 *           Var - variable for the variance
 */
void ave_var(int cc, double &Avg, double &Var, double U)
{
    
    double old_Avg = Avg;
    
    Avg = cc == 0 ? U : Avg + ((U - Avg) / cc);
    
    Var = cc == 0 ? 0 : Var * cc + (U - old_Avg) * (U - Avg);
    Var /= (cc + 1);
}

/* save_data(N, nT, data, h, J)
 * saves calculated quantities to file
 *
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           nT - number of temperature steps
 *           data - calculated values and averages
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 */

void save_data(int N, double * data)
{
    
    double T = data[8];
    
    /*open file for writing*/
    ofstream hfile;
    std::ostringstream hf;
    hf <<"data_"<<N<<".txt" ;
    std::string hmf = hf.str();
    hfile.open(hmf.c_str(), ios_base::app);
    
    if (T == 1){
        hfile << "T\t<M>\tM_Var\t<M^2>\tM^2_Var\t"
        << "<E>\tE_Var\t<E^2>\tE^2_Var\t"
        << "X\tC_v\n" ;
    }
    
    double N3 = (double)N*N*N;
    
    
    double MAvg = data[0], HAvg=data[4], MVar= data[1], HVar = data[5];
    double M2Avg = data[2], H2Avg=data[6], M2Var= data[3], H2Var=data[7];
    
    
    double beta = 1.0/T;
    double beta2 = beta*beta;
    
    //Susceptibility
    double X = beta*N3*( M2Avg - (MAvg*MAvg));
    
    //Specific Heat
    double Cv = (beta2/N3)*(H2Avg - (HAvg*HAvg));
    
    //write to file
    hfile << T << "\t"
    << MAvg << "\t" << MVar << "\t" << M2Avg << "\t" << M2Var << "\t"
    << HAvg << "\t" << HVar << "\t" << H2Avg << "\t" << H2Var << "\t"
    << X    << "\t" << Cv   << "\n" ;
    
    
    
    hfile.close();
    
    return;
}
void fft(fftw_complex *in, fftw_complex *out, int L)
{
    
    
    fftw_plan plan = fftw_plan_dft_1d(L, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();
    
    return;
    
}

void ifft(fftw_complex *in, fftw_complex *out, int L)
{
    
    fftw_plan plan = fftw_plan_dft_1d(L, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();
    
    
    for (int i = 0; i < L; ++i) {
        out[i][REAL] /= L;
        out[i][IMAG] /= L;
    }
}


/* aucf(si_hist, MAvg, MCSteps, T)
 * calculates the magnetization autocorrelation  at
 * each simulation temperature and saves to file
 *
 * takes in: si_hist - list of magnetization values at every timestep
 *           MAvg - Average calculated Magnetizatization
 *           MCSteps - number of Monte Carlo steps
 *           T - simulation temperature
 */

void aucf(double * si_hist, double MAvg, int MCSteps, double T)
{
 
    ofstream acfile;
    std::ostringstream acf;
    acf <<"aucf"<<T<<".txt" ;
    std::string aucf = acf.str();
    acfile.open(aucf.c_str());
    
    
    fftw_complex Mp[MCSteps];
    fftw_complex Mp_w[MCSteps];
    fftw_complex Xw[MCSteps];
    fftw_complex Xt[MCSteps];
    
    int i,t;
#pragma omp parallel for private (i) shared (MAvg)    
    for (i = 0; i < MCSteps; i++){
        Mp[i][REAL] = si_hist[i] - MAvg;
        Mp[i][IMAG] = 0;
    }

    fft(Mp, Mp_w, MCSteps);
#pragma omp parallel for     
    for (i = 0; i < MCSteps; i++){
        Xw[i][REAL]= pow(abs(Mp_w[i][REAL]), 2.);
        Xw[i][IMAG] = 0;
    }
    
    ifft(Xw, Xt, MCSteps);
    
    for (t = 0; t < MCSteps; t++){
    
       acfile << t  << "\t" << Xt[t][REAL] << "\t" << Xt[t][REAL] / Xt[0][REAL] << "\n"; 
        
    }
    
 return;
    
}

int main(int argc, const char * argv[])
{
    
    int N = 1024;
    
    //external magnetic field
    double h = 0.0;
    
    //coupling constant
    double * J = new double [2];
    
    //in eV 
    J[0] = 2.04731206e-01;   J[1] = 6.94636080e-02;
   
    //in Ry
    J[0] = J[0]/13.605698065894 ;
    J[1] = J[1]/13.605698065894 ;
    //temperature value (K)
    double T;
    
    //boltzmann constant(eV K^-1)
    //double kB = 8.617652e-5;
    
    //boltzmann constant(Ry K^-1)
    double kB = 0.0000063306;
    
    
    //number of Monte Carlo moves
    int MCSteps = 50000;
    
    //number of equilibration steps
    int eqSteps = 0.2*MCSteps;

    double ** si = new double * [N];
    double a = 5.3970578; //a.u.
    double L;
    
    initialize(N, L , 1, si, a);
    
    double * si_hist = new double [MCSteps+eqSteps];
   
    
    for (int i = 0; i < eqSteps+MCSteps; i++) 
    {
        si_hist[i] = 0.0;
        
    }
    
    
    double * data = new double  [10];
    for (int i =0; i < 10; i++) data[i] = 0.0;
    
    //initial Energy and Magnetization
    double M = Magnetization(N, si);
    double H = Energy( N, L, a, si, h, J);
    
    
    //loop over T values
    for (T =500; T <= 1200; T+=50)
    {
    
    
        
        //Average and Variance
        double MAvg = 0.0, HAvg = 0.0, MVar=0.0, HVar=0.0;
        double M2Avg=0.0, H2Avg=0.0, M2Var=0.0, H2Var=0.0;

        double kT = kB*T;
        
        
        //Equilibration steps (burn-in)
        for (int step = 0; step < eqSteps; step++)
        {
            MC_move(N, L, a, si, h, J, M, H, kT);
           
        }
        
        
        //Production steps
        for (int step = 0; step < MCSteps; step++)
        {
            /*each move is a complete lattice sweep - N*N*N steps */
            MC_move(N, L, a, si, h, J, M, H, kT);
            
            ave_var(step, MAvg, MVar, M);
            ave_var(step, HAvg, HVar, H);
            ave_var(step, M2Avg, M2Var, M*M);
            ave_var(step, H2Avg, H2Var, H*H);
            
            si_hist[step] = M;
            
            // cout << step << "\t" <<M <<"\t" << MAvg << "\t" << HAvg << "\n";
        
        }//end loop over MCSteps
        
        //accumulate data
        data[0] = MAvg; data[1] = MVar;
        data[2] = M2Avg; data[3] = M2Var;
        data[4] = HAvg; data[5] = HVar;
        data[6] = H2Avg; data[7] = H2Var;
        data[8] = kT;
        
        cout << "saving to file ... \n";
        save_data(N, data);
        
        
        
        aucf(si_hist, MAvg, MCSteps, T);
        
        
        cout << "Done : T = " << T << "\t M = " << M << "\t H = " << H << "\n";
        
        
    }//end loop over T
    


    return 0;
}



