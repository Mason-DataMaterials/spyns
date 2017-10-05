//ising 2.2
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <omp.h>

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
int neighbors2(int N, int *** si, int * point){
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    int xp = (x+N-1) % N;
    int xm = (x+1) % N;
    int yp = (y+N-1) % N;
    int ym = (y+1) % N;
    int zp = (z+N-1) % N;
    int zm = (z+1) % N;
    
    int sj = si[xp][y][zp] + si[x][yp][zp] + si[xm][y][zp] + si[x][ym][zp] + si[xp][y][zm] + si[x][yp][zm] + si[xm][y][zm]+ + si[x][ym][zm];
    
    return sj;
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
int neighbors(int N, int *** si, int * point){
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    int xp = (x+N-1) % N;
    int xm = (x+1) % N;
    int yp = (y+N-1) % N;
    int ym = (y+1) % N;
    int zp = (z+N-1) % N;
    int zm = (z+1) % N;
    
    int sj = si[xp][y][z] + si[xm][y][z] + si[x][yp][z] + si[x][ym][z] + si[x][y][zp] + si[x][y][zm];
    
    return sj;
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
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                M = M+(double)si[i][j][k];
            
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
double Energy( int N, int *** si, double h, double J){
    
    double H = 0.0;
    int * point = new int [3];
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                point[0] = i; point[1] = j; point[2]=k;
                
                int sj = neighbors(N, si, point);
                H = H + ( -(J*si[i][j][k]*sj) - h*si[i][j][k] );
            }
            
            free (point);
        
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
double deltaE (int N, int *** si, double h, double J, int * point)
{
    
    int sj = neighbors(N, si , point);
    double dE = -2.0 * ( (J*si[point[0]][point[1]][point[2]]*sj)
    + (h*si[point[0]][point[1]][point[2]]) );
    
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
void MC_move(int N, int *** si, double h, double J,
             double &M, double &H, double T)
{
    
    double beta = 1.0/T;
    uniform_int_distribution<int> nsites_gen(0, N - 1);
    
    int i = nsites_gen(engine);
    int j = nsites_gen(engine);
    int k = nsites_gen(engine);
    
    int * point = new int [3];
    point[0] = i; point[1] = j; point[2]=k;
    
    double dE = deltaE(N, si, h, J, point);
    double dm = 0.0, dh = 0.0;
    
    double p = T*log(dist(engine));
    
    if (dE > p)
    {
        
        si[i][j][k] *= -1;
        dm = 2.0*si[i][j][k];
        dh = -0.5*dE;
    }
    
    M = M+dm;
    H = H+dh;
    
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
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                cfile << i <<"\t" << j << "\t" << k << "\t"
                << si[i][j][k] << "\n";
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

void save_data(int N, double * data, double h, double J)
{
    
    double T = data[8];
    
    /*open file for writing*/
    ofstream hfile;
    std::ostringstream hf;
    hf <<T<<".txt" ;
    std::string hmf = hf.str();
    hfile.open(hmf.c_str());
    
    hfile << "T\t<M>\tM_Var\t<M^2>\tM^2_Var\t"
    << "<E>\tE_Var\t<E^2>\tE^2_Var\t"
    << "X\tC_v\n" ;
    
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

/* aucf(N, MAvg, si_hist, MCSteps, T)
 * calculates spin time autocorrelation funtion and saves
 * to file.
 * takes in: N - length of one side
 *           MAvg - Average value of M at temperature T
 *           si_hist - list of M values at each simulation step
 *           MCSteps - number of MC simulation steps
 *           T - temperature values
 */
void aucf(int N, double MAvg, double ** si_hist, int MCSteps, double T){
    
    int start = 0.0;
    int stop = MCSteps;
    int sample_rate = 10;
    int scnt = 0;
    
    
    ofstream acfile;
    std::ostringstream acf;
    acf <<"aucf"<<T<<".txt" ;
    std::string aucf = acf.str();
    acfile.open(aucf.c_str());
    
    int tmax = MCSteps;
    double * chi = new double [tmax];
    double chi_Avg, chi_Var;
    
    for (int t = 0; t < tmax; t++){
        
        double dT = 1./(double)(tmax - t);
        
        double term1 = 0.0;
        double term2 = 0.0;
        double term3 = 0.0;
        
        for (int tp = 0; tp < (tmax-t); tp++)
        {
            
            term1 += dT * si_hist[tp][0]  * si_hist[tp+t][0] ;
            term2 += dT * si_hist[tp][0];
            term3 += dT * si_hist[tp+t][0];
            
        }
        
        chi[t] = term1 - (term2 * term3);
        ave_var(t, chi_Avg, chi_Var, chi[t]);
        
        acfile << t << "\t" << chi[t] << "\t" << chi_Avg << "\t" << chi_Var << "\t" << term1-(MAvg*MAvg) << "\n";
        
    }
    
    acfile.close();
    
    return;
}


int main(int argc, char *argv[]){
    
    int N = 10;                 //lattice dimension, lenght of single side
    double h = 1.0;             //external magnetic field
    double J = 2.0;             //coupling constant
    
    double T;                   //temperature value
    
    if ( argc != 5 ) {
        cout<<"usage: "<< argv[0] <<" N T h J \n";
        return 0;  
    }
    else { 
        N = (int) atof(argv[1]);
        T = atof(argv[2]);
        h = atof(argv[3]);
        J = atof(argv[4]);
        
    }
    
    
    
    int MCSteps = 1000*N*N*N;    //number of Monte Carlo moves
    int eqSteps = 0.5*MCSteps;  //number of equilibration steps
    int sweep = N*N*N;          //number of steps for a complete lattice sweep
    
    
 /*   
    cout << "Ising Model MC Simulation \n";
    cout << "PARAMETERS: N = " << N << "\n";
    cout << "            T = " <<fixed << T << "\n";
    cout << "            B = " <<fixed << h << "\n";
    cout << "            J = " <<fixed << J << "\n";
    cout << "No. of MC Steps = " << MCSteps <<"\n\n"; 
    
    cout << "Initializing ... \n"; 
*/
    //spin sites
    int *** si = new int ** [N];
    initialize(N, si);
    
    int nt = omp_get_max_threads();
    
    double ** si_hist = new double * [MCSteps+eqSteps];
    for (int i = 0; i < eqSteps+MCSteps; i++) si_hist[i] = new double[3];
    for (int i = 0; i < (eqSteps+MCSteps); i++)
    {   si_hist[i][0] = 0.0;
        si_hist[i][1] = 0.0;
        si_hist[i][2] = 0.0;
    }
    
    uniform_int_distribution<int> nsites_gen(0, N - 1);
    
    int i = nsites_gen(engine);
    int j = nsites_gen(engine);
    int k = nsites_gen(engine);
    
    int * point = new int [3];
    point[0] = 1; point[1] = 2; point[2]=3;
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    int xp = (x+N-1) % N;
    int xm = (x+1) % N;
    int yp = (y+N-1) % N;
    int ym = (y+1) % N;
    int zp = (z+N-1) % N;
    int zm = (z+1) % N;
    
    cout << x << "\t" << y << "\t" << z << "\t" << 0 << "\n";
    cout << xp << "\t" << y << "\t" << z << "\t" << 1 << "\n";
    cout << xm << "\t" << y << "\t" << z << "\t" << 1 << "\n";
    cout << x << "\t" << yp << "\t" << z << "\t" << 1 << "\n";
    cout << x << "\t" << ym << "\t" << z << "\t" << 1 << "\n";
    cout << x << "\t" << y << "\t" << zp << "\t" << 1 << "\n";
    cout << x << "\t" << y << "\t" << zm << "\t" << 1 << "\n";
    
    cout << xp << "\t" << y << "\t" << zp << "\t" << 2 << "\n";
    cout << x << "\t" << yp << "\t" << zp << "\t" << 2 << "\n";
    cout << xm << "\t" << y << "\t" << zp << "\t" << 2 << "\n";
    cout << x << "\t" << ym << "\t" << zp << "\t" << 2 << "\n";
    cout << xp << "\t" << y << "\t" << zm << "\t" << 2 << "\n";
    cout << x << "\t" << yp << "\t" << zm << "\t" << 2 << "\n";
    cout << xm << "\t" << y << "\t" << zm << "\t" << 2 << "\n";
    cout << x << "\t" << ym << "\t" << zm << "\t" << 2 << "\n";
  /*  
    int sj = si[xp][y][z] + si[xm][y][z] + si[x][yp][z] + si[x][ym][z] + si[x][y][zp] + si[x][y][zm];
    
     int sj = 
     si[xp][y][zp] + 
     si[x][yp][zp] + 
     si[xm][y][zp] + 
     si[x][ym][zp] + 
     si[xp][y][zm] + 
     si[x][yp][zm] + 
     si[xm][y][zm] + 
     si[x][ym][zm];
    
  */  
 /*   
    double * data = new double  [10];
    for (int i =0; i < 10; i++) data[i] = 0.0;
    
    //initial Energy and Magnetization
    double M = Magnetization(N, si);
    double H = Energy(N, si, h, J);
    
    //Average and Variance
    double MAvg, HAvg, MVar, HVar;
    double M2Avg, H2Avg, M2Var, H2Var;
    
    cout << "Initial M = " << M << "\n"; 
    cout << "Initial H = " << H << "\n"; 
    
    cout << "\nBurn In ... \t";
    //Equilibration steps (burn-in)
    for (int step = 0; step < eqSteps; step++)
    {
        MC_move(N, si, h, J, M, H, T);
        cout << step << "\t" << M << "\n";
    }
 /*   
    cout << "\nProdcution ... \n";
    //Production steps
    for (int step = 0; step < MCSteps; step++){
        
        MC_move(N, si, h, J, M, H, T);
        
        
        /*average after every complete lattice sweep - N*N*N steps */
/*        
        if (step % sweep == 0)
        {
           ave_var(step, MAvg, MVar, M);
            ave_var(step, HAvg, HVar, H);
            ave_var(step, M2Avg, M2Var, M*M);
/*            ave_var(step, H2Avg, H2Var, H*H);
        }
        
        
        si_hist[step][0] = M;
        si_hist[step][1] = MAvg;
        si_hist[step][2] = M2Avg;
        
        
    }//end loop over MCSteps
    
 /*   
    //accumulate data
    data[0] = MAvg; data[1] = MVar;
    data[2] = M2Avg; data[3] = M2Var;
    data[4] = HAvg; data[5] = HVar;
    data[6] = H2Avg; data[7] = H2Var;
    data[8] = T;
    
    //compute time correlation at each temperature T
//    aucf(N, MAvg, si_hist, MCSteps, T);
    
    cout << "saving to file ... \n";
    save_data(N, data, h, J);
    
    cout << "Done for T = " << T << "\n" ;
    cout << "Final M = " << M << "\n"; 
    cout << "Final H = " << H << "\n"; 
    

    free (si);
    free (data);
    free (si_hist);
*/
 
    return 0;

}
