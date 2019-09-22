#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <omp.h>
#include <fftw3.h>
#include <vector>


#define REAL 0
#define IMAG 1


using namespace std;


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
    acf <<"aucf_files/aucf_"<<T<<".txt" ;
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
    
    for (t = 0; t < int(0.5*MCSteps); t++){
    
       acfile << t  << "\t" << Xt[t][REAL] << "\t" << Xt[t][REAL] / Xt[0][REAL] << "\n"; 
        
    }
    
 return;
    
}


int main(int argc, char *argv[]){
    
    
     ifstream infile("aucf.in");
     char * temp = new char[30];
     std::string dir = "aucf_files/";
     
     
     int N;
     int MCSteps;
     double min_temp  =  0.0;
     double max_temp  =  0.0;
     double temp_step = 0.0;
    
    
    if (infile){
        infile>>temp>>temp>>temp>>temp>>N;
        infile>>temp>>temp>>temp>>temp>>MCSteps;
            
            infile>>temp>>temp>>temp>>min_temp;
            infile>>temp>>temp>>temp>>max_temp;
            infile>>temp>>temp>>temp>>temp_step;
            
            infile.close();
    }else {cout << "Found it, NOT! \n"; return 0;}
    
    
    int nsteps = int ((max_temp - min_temp)/temp_step); 
    double * MAvg = new double [nsteps+1];
    double * HAvg = new double [nsteps+1];
    double * X = new double [nsteps+1];
    double * Cv = new double [nsteps+1];
    double * temp_values= new double [nsteps+1];
    
     for (int i = 0; i <= nsteps; i++){
        MAvg[i] = 0.0;
        HAvg[i] = 0.0;
        X[i] = 0.0;
        Cv[i] = 0.0;
        temp_values[i] = 0.0;
    }
    
    
    //read in averages
    
    ifstream avfile("data.txt");
    double hold = 0.0;
    if(avfile){
        int i = 0;
        while(!avfile.eof())
        {
            avfile>>temp_values[i]>>MAvg[i]>>hold>>hold>>hold>>HAvg[i]>>hold>>hold>>hold>>X[i]>>Cv[i];
            i++;
        }
    }
   
    double * M_hist = new double [MCSteps+1];
    double * H_hist = new double [MCSteps+1];
    
    for (int i = 0; i <= MCSteps; i++){
        M_hist[i] = 0.0;
        H_hist[i] = 0.0;
    }
    
    double T;
    int i = 0;
    for (T = min_temp; T <= max_temp; T+=temp_step)
    {
    
        //if (T==1200) {i++;continue;}
        
        ostringstream t;
        t<<T;
        std::string filename = dir + t.str();

        int step;
        ifstream afile(filename);
        if(afile){
            while(!afile.eof()){
                afile>>step>>M_hist[step]>>H_hist[step];
                if (step >=100000) break;
            }
            
            aucf(M_hist, MAvg[i], step, T);

            cout << "Done: "<<temp_values[i] << "\n";
        }
        
        else {cout << filename << " not found! \n";}
         
        
        
        i++;
        
    }
    
     
    return 0;
}
