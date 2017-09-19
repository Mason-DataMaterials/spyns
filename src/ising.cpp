//ising 2.0
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>


using namespace std;

/* rnd(a,b)
 * Random number generator 
 * using the Mersenne Twister
 * takes in: a - lower limit, real number
 *           b - upper limit, real number
 *
 * returns:  x - random number between a & b
 */
double rnd (double a, double b){
    
    /*http://stackoverflow.com/a/22935167*/
    random_device rd{};    
    mt19937 engine{rd()};
    uniform_real_distribution<double> dist{a, b};
    double x = dist(engine);
    return x;
}

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
                double p  = rnd(0.0, 1.0);
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
    
    int sj = si[xp][y][z] + si[xm][y][z] + si[x][yp][z] 
                      + si[x][ym][z] + si[x][y][zp] + si[x][y][zm];
    
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
    
    int i = (int) rnd(0.0, (double)N)%N; 
    int j = (int) rnd(0.0, (double)N)%N; 
    int k = (int) rnd(0.0, (double)N)%N;
    
    int * point = new int [3];
    point[0] = i; point[1] = j; point[2]=k;
    
    double dE = deltaE(N, si, h, J, point);
    double dm = 0.0, dh = 0.0;
    
    double p = T*log(rnd(0.0, 1.0));
    
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

/*main program to calculate M,H values*/
int main(){
    
    int N = 10;                  //lattice dimension, lenght of single side
    int MCSteps = 500*N*N*N;    //number of Monte Carlo moves
    int eqSteps = 0.5*MCSteps;  //number of equilibration steps
    
    
    double h = 1.0;             //external magnetic field
    double J = 2.0;             //coupling constant
    double T;                   //temperature value
    
    //spin sites
    int *** si = new int ** [N];
    initialize(N, si);
    
    //initial Energy and Magnetization
    double M = Magnetization(N, si);
    double H = Energy(N, si, h, J);
    
    //Average and Variance 
    double MAvg, HAvg, MVar, HVar;
    double M2Avg, H2Avg, M2Var, H2Var;
    
    //reassurance
    cout << "Running h = " << h << "\n";
    
    /*open file for writing*/
    ofstream tfile;
    std::ostringstream tf;
    tf <<"h="<<h<<"_mvst.txt" ;
    std::string tmf = tf.str();
    tfile.open(tmf.c_str());
    tfile << "T  \t  <M>  \t  M_Var  \t  <M^2>  \t  M^2_Var \t" 
    << "<E>  \t  E_Var  \t  <E^2>  \t  E^2_Var \t" 
    << "X    \t  C_v   \n" ;
    
    /*loop over T values*/
    for (T=1.0; T <= 20.; T+=1){
        
        
        //Equilibration steps (burn-in)
        for (int step = 0; step < eqSteps; step++)
            MC_move(N, si, h, J, M, H, T);
        
        //Production steps
        for (int step = 0; step < MCSteps; step++){
            
            MC_move(N, si, h, J, M, H, T);
            
            //running average and variance calculations 
            ave_var(step, MAvg, MVar, M);
            ave_var(step, HAvg, HVar, H);
            ave_var(step, M2Avg, M2Var, M*M);
            ave_var(step, H2Avg, H2Var, H*H);
            
            
        }//end loop over MCSteps
        
        
        //Quantities
        
        double beta = 1.0/T;
        double beta2 = beta*beta;
        double N3 = (double)N*N*N;
        
        //Susceptibility
        double X = beta*N3*( M2Avg - (MAvg*MAvg));
        
        //Specific Heat
        double Cv = (beta2/N3)*(H2Avg - (HAvg*HAvg));
        
        /*write to file*/
        tfile << T << "\t"  
        << MAvg << "\t" << MVar << "\t" << M2Avg << "\t" << M2Var << "\t" 
        << HAvg << "\t" << HVar << "\t" << H2Avg << "\t" << H2Var << "\t" 
        << X    << "\t" << Cv   << "\n" ;
        
        //reassure      
        cout << "Done h = " << h << ", T = " << T <<"\n";
        
    }//end for loop over T
    
    tfile.close();
    
    //reassure
    cout << "Done h = " << h << "\n";
    
    free (si);
    
    return 0;
    
}
