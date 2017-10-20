/*ising 2.4
 * compile and run with:
 * >>g++ ising2.cpp -Wall -lfftw3 -lm -fopenmp -o ising
 * >>./ising
 * 
*/ 
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
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

/* initalize(N, si)
 * Generates initial random spin configuration
 * for an N*N*N lattice grid. Each lattice site
 * is addressed as si[x][y][z] with 0 <= x,y,z < N
 * takes in: N  - length of one side
 *           si - 3D lattice of spin sites
 *
 */
void initialize(int N, int *** si){
    
    for (int x = 0; x < N; x++){
        si[x] = new int * [N];
        
        for (int y = 0; y < N; y++){
            si[x][y] = new int [N];
            
            for (int z = 0; z < N; z++){
                double p  = dist(engine);
                if (p > 0.5) si[x][y][z] = 1;
                else si[x][y][z] = -1;
            }
        }
    }
    
    
    
    return;
}

/* neighbors(N, si, point)
 * calculates the sum of spins for the
 * 6 neighboring sites of a point in the
 * grid.
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           point - site of interest
 * returns: sj - sum of the spins of neighbors of point
 */
void neighbors(int N, int *** si, int * point, int * sj){
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    int xp = (x+N-1) % N;
    int xm = (x+1) % N;
    int yp = (y+N-1) % N;
    int ym = (y+1) % N;
    int zp = (z+N-1) % N;
    int zm = (z+1) % N;
    
    int sj1 = si[xp][y][z] + si[xm][y][z] + si[x][yp][z] + si[x][ym][z] + si[x][y][zp] + si[x][y][zm];
    int sj2 = si[xp][y][zp] + si[x][yp][zp] + si[xm][y][zp] + si[x][ym][zp]
    + si[xp][y][zm] + si[x][yp][zm] + si[xm][y][zm]+ + si[x][ym][zm];
    
    sj[0] = sj1;
    sj[1] = sj2;
    
    return;
}

/* Magnetization(N, si)
 * sums the spins of each site and returns total
 * spin/magnetization
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 * returns: M - sum of all spins
 */
double Magnetization(int N, int *** si){
    
    double M = 0.0;
    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                M = M+(double)si[i][j][k];
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
double Energy( int N, int *** si, double h, double * J){
    
    double H = 0.0;
    int * point = new int [3];
    int * sj = new int [2];
    sj[0] = sj[1] = 0;
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                point[0] = i; point[1] = j; point[2]=k;
                
                neighbors(N, si, point, sj);
                H = H + ( -(si[i][j][k]*(J[0]*sj[0] + J[1]*sj[1])) - h*si[i][j][k] );
            }
    }
    
    
    free (point);
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
double deltaE (int N, int *** si, double h, double * J, int * point)
{
    int * sj = new int [2];
    sj[0] = sj[1] = 0;
    
    neighbors(N, si , point, sj);
    double dE = 2.0 * ( (si[point[0]][point[1]][point[2]]*(J[0]*sj[0] +J[1]*sj[1]))
    + (h*si[point[0]][point[1]][point[2]]) );
    
    
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
void MC_move(int N, int *** si, double h, double * J,
             double &M, double &H, double T)
{
    
    double beta = 1.0/T;
    int * point = new int [3];
    
    for (int sweep = 0; sweep <N*N*N; sweep++)
    {
        //random spin site
        uniform_int_distribution<int> nsites_gen(0, N - 1);
        int i = nsites_gen(engine);
        int j = nsites_gen(engine);
        int k = nsites_gen(engine);
        
        
        point[0] = i; point[1] = j; point[2]=k;
        
        double dE = deltaE(N, si, h, J, point);
        double dm = 0.0, dh = 0.0;
        
        double p = dist(engine);
        
        if (dE < 0.0)
        {
            
            si[i][j][k] *= -1;
            dm = 2.0*si[i][j][k];
            dh = dE;
        }
        else if (p < exp(-dE*beta))
        {
            si[i][j][k] *= -1;
            dm = 2.0*si[i][j][k];
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

/* configuration(N, si, filename)
 * saves a snapshot of the configuration
 * to file named "filename".
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           filename - name of file to save to
 */
void configuration(int N, int *** si, string filename){
    ofstream cfile;
    cfile.open(filename.c_str());
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                cfile << i <<"\t" << j << "\t" << k << "\t"
                << si[i][j][k] << "\n";
            }
    }
    
    cfile.close();
        
    return;
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


int main(int argc, char *argv[]){
    
    int N = 20;                 //lattice dimension, lenght of single side
    double h = 1.0;             //external magnetic field
    
    double * J = new double [2];
    J[0] =1.0; J[1] = 0.5;              //coupling constant
    
    double T;                   //temperature value
    
    
    int MCSteps = 10000;    //number of Monte Carlo moves
    int eqSteps = 0.2*MCSteps;  //number of equilibration steps
   
    
    //spin sites
    int *** si = new int ** [N];
    initialize(N, si);
    
    double * si_hist = new double [MCSteps+eqSteps];
   
    
    for (int i = 0; i < eqSteps+MCSteps; i++) 
    {
        si_hist[i] = 0.0;
        
    }
    
    
    double * data = new double  [10];
    for (int i =0; i < 10; i++) data[i] = 0.0;
    
    //initial Energy and Magnetization
    double M = Magnetization(N, si);
    double H = Energy(N, si, h, J);
    
    //Average and Variance
    double MAvg, HAvg, MVar, HVar;
    double M2Avg, H2Avg, M2Var, H2Var;
    
    
    //loop over T values
    for (T =1; T <= 20; T+=0.5)
    {
        //Equilibration steps (burn-in)
        for (int step = 0; step < eqSteps; step++)
        {
            MC_move(N, si, h, J, M, H, T);
            //cout << step << "\t" <<M <<"\t" << MAvg << "\t" << HAvg << "\n";
        }
        
        
        //Production steps
        for (int step = 0; step < MCSteps; step++)
        {
            /*each move is a complete lattice sweep - N*N*N steps */
            MC_move(N, si, h, J, M, H, T);
            
            ave_var(step, MAvg, MVar, M);
            ave_var(step, HAvg, HVar, H);
            ave_var(step, M2Avg, M2Var, M*M);
            ave_var(step, H2Avg, H2Var, H*H);
            
            si_hist[step] = M;
            
            
        }//end loop over MCSteps
        
        //accumulate data
        data[0] = MAvg; data[1] = MVar;
        data[2] = M2Avg; data[3] = M2Var;
        data[4] = HAvg; data[5] = HVar;
        data[6] = H2Avg; data[7] = H2Var;
        data[8] = T;
        
        cout << "saving to file ... \n";
        save_data(N, data);
        
        aucf(si_hist, MAvg, MCSteps, T);
        
        
        cout << "Done : T = " << T << "\t M = " << M << "\t H = " << H << "\n";
        
        
    }//end loop over T
    
    
    
    
    return 0; 
}





