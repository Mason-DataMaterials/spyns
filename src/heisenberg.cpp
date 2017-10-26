
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
void initialize(int N, double **** si){
    
    for (int x = 0; x < N; x++)
    {
        si[x] = new double ** [N];
        
        for (int y = 0; y < N; y++){
            
            si[x][y] = new double * [N];
            
            for (int z = 0; z < N; z++){
                
                si[x][y][z] = new double [5];
                
                //(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
                
                double theta = 1.0/cos(1.0-2.0 * (double) dist(engine));
                
                double phi = M_PI*2.0* (double)dist(engine);
                
                //theta
                si[x][y][z][3] = theta;
    
                //phi
                si[x][y][z][4] = phi;
                
                //Sx
                si[x][y][z][0] = sin(si[x][y][z][3])*cos(si[x][y][z][4]);
                //Sy
                si[x][y][z][1] = sin(si[x][y][z][3])*sin(si[x][y][z][4]);
                //Sz
                si[x][y][z][2] = cos(si[x][y][z][3]);
                
                
                
                
                //cout << x << "\t" << y << "\t" << z << "\t" << 
                //si[x][y][z][0] << "\t"<< si[x][y][z][1] <<"\t"<< si[x][y][z][2]  << "\n";
                
            }
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
double Magnetization(int N, double **** si){
    
    double M = 0.0;
    
    double * m = new double [3];
    for(int i = 0; i < 3; i++) m[i] = 0.0;
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                for(int p = 0; p <3; p++)
                {
                    m[p] += si[i][j][k][p];
                }
            }
        }
    }
    
    for (int i = 0; i < 3; i++) m[i] /= (double) N;
    M = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        
    return M;
    
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
void neighbors(int N, double **** si, int * point, int i, int * sj){
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    int xp = (x+N-1) % N;
    int xm = (x+1) % N;
    int yp = (y+N-1) % N;
    int ym = (y+1) % N;
    int zp = (z+N-1) % N;
    int zm = (z+1) % N;
    
    int sj1 = si[xp][y][z][i] + si[xm][y][z][i] + si[x][yp][z][i] 
                   + si[x][ym][z][i] + si[x][y][zp][i] + si[x][y][zm][i];
                   
    int sj2 = si[xp][y][zp][i] + si[x][yp][zp][i] + si[xm][y][zp][i] + si[x][ym][zp][i]
                  + si[xp][y][zm][i] + si[x][yp][zm][i] + si[xm][y][zm][i] +  si[x][ym][zm][i];
    
    sj[0] = sj1;
    sj[1] = sj2;
    
   
    return;
}

/* Energy(N, si, h, J)
 * calculates total energy of the grid
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 * returns: H - total energy of grid
 */
double Energy(int N, double **** si, double h, double * J)
{
    double H  = 0.0;
    
    double * Hi = new double [3];
    for (int i =0; i <3; i++) Hi[i] = 0.0;
    
    int * point = new int [3];
    
    int * sj_x = new int [2];
    sj_x[0] = sj_x[1] = 0;
    
    int * sj_y = new int [2];
    sj_y[0] = sj_y[1] = 0;
    
    int * sj_z = new int [2];
    sj_z[0] = sj_z[1] = 0;
            
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                point[0] = i; point[1] = j; point[2]=k;
                
                neighbors(N, si, point, 0, sj_x);
                neighbors(N, si, point, 1, sj_y);
                neighbors(N, si, point, 2, sj_z);
            
                Hi[0] = Hi[0] + ( (si[i][j][k][0]*-(J[0]*sj_x[0] + J[1]*sj_x[1]) ) - h*si[i][j][k][0] );
                Hi[1] = Hi[1] + ( (si[i][j][k][1]*-(J[0]*sj_y[0] + J[1]*sj_y[1]) ) - h*si[i][j][k][1] );
                Hi[2] = Hi[2] + ( (si[i][j][k][2]*-(J[0]*sj_z[0] + J[1]*sj_z[1]) ) - h*si[i][j][k][2] );
            }
    }
    
    H = Hi[0]+Hi[1]+Hi[2];
    
    return H/2.0;
    
}
 


 
 /* dE_X(N, si, h, J, point, X)
 * calculates the change in energy if a spin site is flipped
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 *           point - site of interest
 * returns: dE - change in energy from changing value X - theta or phi
 */

double dE_theta (int N, double **** si, double h, double * J, int * point, double new_theta)
{
    int * sj_x = new int [2];
    sj_x[0] = sj_x[1] = 0;
    
    int * sj_y = new int [2];
    sj_y[0] = sj_y[1] = 0;
    
    int * sj_z = new int [2];
    sj_z[0] = sj_z[1] = 0;
    
    neighbors(N, si, point, 0, sj_x);
    neighbors(N, si, point, 1, sj_y);
    neighbors(N, si, point, 2, sj_z);
    
    double theta = si[point[0]][point[1]][point[2]][3];
    double phi = si[point[0]][point[1]][point[2]][4];
    
    
    double sx  = sin(new_theta)*cos(phi) - sin(theta)*cos(phi);
    double sy =  sin(new_theta)*sin(phi) - sin(theta)*sin(phi);
    double sz =  cos(new_theta) - cos(theta);
    
    
    double dE_x =  (sx * -(J[0]*sj_x[0] + J[1]*sj_x[1]) ) - h*sx ;
    
    double dE_y =  (sy * -(J[0]*sj_y[0] + J[1]*sj_y[1]) ) - h*sy ;
    
    double dE_z =  (sz * -(J[0]*sj_z[0] + J[1]*sj_z[1]) ) - h*sz ;
    
    double dE = dE_x + dE_y + dE_z;
    
    
    free (sj_x);
    free (sj_y);
    free (sj_z);
    
    return dE;
}

double dE_phi (int N, double **** si, double h, double * J, int * point, double new_phi)
{
    int * sj_x = new int [2];
    sj_x[0] = sj_x[1] = 0;
    
    int * sj_y = new int [2];
    sj_y[0] = sj_y[1] = 0;
    
    int * sj_z = new int [2];
    sj_z[0] = sj_z[1] = 0;
    
    neighbors(N, si, point, 0, sj_x);
    neighbors(N, si, point, 1, sj_y);
    neighbors(N, si, point, 2, sj_z);
    
    double theta = si[point[0]][point[1]][point[2]][3];
    double phi = si[point[0]][point[1]][point[2]][4];
    
    
    double sx  = sin(theta)*cos(new_phi) - sin(theta)*cos(phi);
    double sy =  sin(theta)*sin(new_phi) - sin(theta)*sin(phi);
    
    
    
    double dE_x =  (sx * -(J[0]*sj_x[0] + J[1]*sj_x[1]) ) - h*sx ;
    
    double dE_y =  (sy * -(J[0]*sj_y[0] + J[1]*sj_y[1]) ) - h*sy ;
    
    
    double dE = dE_x + dE_y;
    
    
    free (sj_x);
    free (sj_y);
    free (sj_z);
    
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
void MC_move(int N, double **** si,  double h, double *J, double &M, double &H, double T)
{
    double phi, theta, p;
    
    double beta = 1.0/T;
    int * point = new int [3];
    
    for (int sweep = 0; sweep < N*N*N; sweep++)
    {
        //random spin site
        uniform_int_distribution<int> nsites_gen(0, N - 1);
        int i = nsites_gen(engine);
        int j = nsites_gen(engine);
        int k = nsites_gen(engine);
        
        
        point[0] = i; point[1] = j; point[2]=k;
        
        
        //theta update, phi constant
        theta = 1.0/cos(1.0-2.0 * (double)dist(engine));
        
        double dE = dE_theta(N,si, h, J, point, theta);
    
        if (dE <= 0.0)
        {
                //theta
                si[i][j][k][3] = theta;
                //Sx
                si[i][j][k][0] = sin(si[i][j][k][3])*cos(si[i][j][k][4]);
                //Sy
                si[i][j][k][1] = sin(si[i][j][k][3])*sin(si[i][j][k][4]);
                //Sz
                si[i][j][k][2] = cos(si[i][j][k][3]);
            
                H = H+dE;
        }
        else {
            
            p = dist(engine);
            
            if (p <= exp(-dE*beta))
            {
                
                si[i][j][k][3] = theta;
                //Sx
                si[i][j][k][0] = sin(si[i][j][k][3])*cos(si[i][j][k][4]);
                //Sy
                si[i][j][k][1] = sin(si[i][j][k][3])*sin(si[i][j][k][4]);
                //Sz
                si[i][j][k][2] = cos(si[i][j][k][3]);
            
                H = H+dE;
                
            }
            
            
        } //end if for theta
         
        
        //phi update, theta constant
        phi = M_PI*2.0* (double)dist(engine);
        dE = dE_phi(N,si, h, J, point, phi);
        
        
        if (dE <= 0.0)
        {
            
                //phi
                si[i][j][k][4] = phi;
                //Sx
                si[i][j][k][0] = sin(si[i][j][k][3])*cos(si[i][j][k][4]);
                //Sy
                si[i][j][k][1] = sin(si[i][j][k][3])*sin(si[i][j][k][4]);
               
            
                
        }
        else {
            p = dist(engine);
            
            if (p < exp(-dE*beta))
            {
                //phi
                si[i][j][k][4] = phi;
                //Sx
                si[i][j][k][0] = sin(si[i][j][k][3])*cos(si[i][j][k][4]);
                //Sy
                si[i][j][k][1] = sin(si[i][j][k][3])*sin(si[i][j][k][4]);
               
               
                
            }
        }//end if for phi
        
    }//end sweep
    
  
    
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

/* configuration(N, si, filename)
 * saves a snapshot of the configuration
 * to file named "filename".
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           filename - name of file to save to
 */
void configuration(int N, double **** si, string filename){
    ofstream cfile;
    cfile.open(filename.c_str());
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                cfile << i <<"\t" << j << "\t" << k << "\t"
                << si[i][j][k][3] << "\t" << si[i][j][k][4] << "\n";
            }
    }
    
    cfile.close();
        
    return;
}




int main(int argc, char *argv[]){

    
    int N = 12;                         //lattice dimension, lenght of single side
    double h = 1.0;                     //external magnetic field
    
    double * J = new double [2];
    J[0] =1.0; J[1] = 0.5;              //coupling constant
    
    
    int MCSteps = 50000;                //number of Monte Carlo moves
    int eqSteps = 0.2*MCSteps;          //number of equilibration steps
    
    
    double T;  
    
    
    
    double **** si = new double *** [N];//spin configuration
    initialize(N, si);
    
    double * si_hist = new double [MCSteps+eqSteps];
        
    for (int i = 0; i < eqSteps+MCSteps; i++) 
    {
        si_hist[i] = 0.0;
    }
    
    
    double * data = new double  [10];
    for (int i =0; i < 10; i++) data[i] = 0.0;
    
    double M = Magnetization(N, si);
    double H = Energy(N, si, h, J);
    
    //Average and Variance
    double MAvg, HAvg, MVar, HVar;
    double M2Avg, H2Avg, M2Var, H2Var;
    
    cout << "Running Heisenberg for N = " << N << ": \n";
    
    //loop over T values
    for (T =1; T <= 15; T+=1)
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
            
            M = Magnetization(N, si);
            H = Energy(N, si, h, J);
            
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
