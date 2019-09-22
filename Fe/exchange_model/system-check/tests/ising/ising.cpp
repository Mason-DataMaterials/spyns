//ising 2.4
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
void initialize(int N, int * si){
    
    for (int x = 0; x < N; x++)
    {
                double p  = dist(engine);
                if (p > 0.5) si[x] = 1;
                else si[x] = -1;
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
void neighbors(int N, int * si, int * num_nn, int * first_nn, int * second_nn, int * sj){
    
   
    
    int sj1 = 0;
    int sj2 = 0;
    
    for (int i = 0; i < num_nn[0]; i++) sj1+=si[first_nn[i]];
    for (int i = 0; i < num_nn[1]; i++) sj2+=si[second_nn[i]];
    
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
double Magnetization(int N, int * si){
    
    double M = 0.0;
    
    for (int i = 0; i < N; i++)
    {
                M = M+(double)si[i];
    }       
    return M / double (N);
        
}

/* Energy(N, si, h, J)
 * calculates total energy of the grid
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 * returns: H - total energy of grid
 */
double Energy( int N, int * si, int * num_nn, int ** first_nn, int ** second_nn, double h, double * J){
    
    double H = 0.0;
    
    int * sj = new int [2];
    sj[0] = sj[1] = 0;
    
    for (int i = 0; i < N; i++)
    {         
        neighbors(N, si, num_nn, first_nn[i], second_nn[i], sj);
        H = H + ( -(si[i]*(J[0]*sj[0] + J[1]*sj[1])) - h*si[i] );
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
double deltaE (int N, int * si, int * num_nn, int ** first_nn, int ** second_nn, double h, double * J, int point)
{
    int * sj = new int [2];
    sj[0] = sj[1] = 0;
    
    neighbors(N, si, num_nn , first_nn[point], second_nn[point], sj);
    double dE =  2.0 * ( (si[point]* -(J[0]*sj[0] +J[1]*sj[1])) + ( -h*si[point]) );
    
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
void MC_move(int N, int * si, int * num_nn, int ** first_nn, int ** second_nn,
             double h, double * J,
             double &M, double &H, double T)
{
    
    double beta = 1.0/T;

    
    for (int sweep = 0; sweep < N; sweep++)
    {
        //random spin site
        uniform_int_distribution<int> nsites_gen(0, N - 1);
        int point = nsites_gen(engine);
    
        
        
        double dE = deltaE(N, si, num_nn, first_nn, second_nn, h, J, point);
        double dm = 0.0, dh = 0.0;
        
        double p = dist(engine);
        
        if (dE < 0.0)
        {
            si[point] *= -1;
            dm = -2.0*si[point];
            dh = dE;
        }
        else if (p < exp(-dE*beta))
        {
            si[point] *= -1;
            dm = -2.0*si[point];
            dh = dE;
        }
        
        
        
        M = M+dm;
        H = H+dh;
        
    }
    

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
void configuration(int N, int * si, string filename){
    ofstream cfile;
    cfile.open(filename.c_str());
    
    for (int i = 0; i < N; i++)
    {
        cfile << i <<"\t" << si[i] << "\n";
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
    
    
    
    
    double MAvg = data[0], HAvg=data[4], MVar= data[1], HVar = data[5];
    double M2Avg = data[2], H2Avg=data[6], M2Var= data[3], H2Var=data[7];
    
    
    double beta = 1.0/T;
    double beta2 = beta*beta;
    
    //Susceptibility
    double X = beta*double (N) *( M2Avg - (MAvg*MAvg));
    
    //Specific Heat
    double Cv = (beta2/double (N) )*(H2Avg - (HAvg*HAvg));
    
    //write to file
    hfile << T << "\t"
    << MAvg << "\t" << MVar << "\t" << M2Avg << "\t" << M2Var << "\t"
    << HAvg << "\t" << HVar << "\t" << H2Avg << "\t" << H2Var << "\t"
    << X    << "\t" << Cv   << "\n" ;
    
    
    
    hfile.close();
    
    return;
}

void aucf(double ** si_hist, double ** auc, double M, double H, int nsave, int step){
    
    int i = ((step - 1) % nsave) + 1;;
    int index, t;
    
    if (step > nsave) 
    {
        index = i;
        
        for (t=nsave; t >= 1; t--)
        {
            auc[t][0] = auc[t][0] + M*si_hist[index][0];
            auc[t][1] = auc[t][1] + H*si_hist[index][1];
            
            index = index + 1;
            if(index > nsave) index = 1;
            
        }
    }
    
    si_hist[i][0] = M;
    si_hist[i][1] = H;
    
    
    return;
}


void aucf_norm(double ** auc, int nsave, 
               double MAvg, double HAvg, 
               double M2Avg, double H2Avg, 
               double T, int MCSteps)
{
    
    ofstream acfile;
    std::ostringstream acf;
    acf <<"aucf"<<T<<".txt" ;
    std::string aucf = acf.str();
    acfile.open(aucf.c_str());
    
    double NM = 1.0/(MCSteps-nsave);
    double MM = MAvg*MAvg;
    double HH = HAvg*HAvg;
    
    auc[0][0] = M2Avg - MM;
    auc[0][1] = H2Avg - HH;
    
    for (int t = 1;t < nsave;t++)
    {
        auc[t][0] = (auc[t][0]*NM - MM)/auc[0][0];
        auc[t][1] = (auc[t][1]*NM - HH)/auc[0][1];
        
        acfile << T  << "\t" << t << "\t" << auc[t][0] << "\t" << auc[t][1] << "\n";
    }
    
    return;
    
}


int main(int argc, char *argv[]){
    
    
    
    int N ;                 //lattice dimension, lenght of single side
    double h ;               //external magnetic field
    double * J = new double [2];
    double T;                   //temperature value
    
    int MCSteps = 10000;        //number of Monte Carlo moves
    int eqSteps = 0.2*MCSteps;  //number of equilibration steps
    
    int nsave;
    
    char * temp = new char[30];

    ifstream infile("hinput");
    if (infile)
    {
        infile>>temp>>temp>>T;
        infile>>temp>>temp>>h;
        infile>>temp>>temp>>J[0];
        infile>>temp>>temp>>J[1];
        infile>>temp>>temp>>MCSteps;
        infile>>temp>>temp>>nsave;

        infile.close();
    }
    else{ cout << "hinput not found!\n"; return 0; }

    
    int num_1st_nn,num_2nd_nn;
    int * num_nn;
    num_nn = new int [2];


    int ** first_nn;
    int ** second_nn;
    
       ifstream infile2 ("neighbor_lists.txt");
    if (infile2.is_open())
    {

        infile2>>N;
        infile2>>num_1st_nn>>num_2nd_nn;

        num_nn[0] = num_1st_nn;
        num_nn[1] = num_2nd_nn;

        //populate neighbor lists
        first_nn = new int * [N];
        second_nn = new int * [N];

        for (int i = 0; i < N; i++)
        {
            first_nn[i] = new int [num_1st_nn];
            second_nn[i] = new int [num_2nd_nn];
        }  
        

        for (int i = 0; i < N; i++){
             
            for (int j = 0; j < num_1st_nn; j++){
                infile2 >> first_nn[i][j]; 
              }
            for (int j = 0; j < num_2nd_nn; j++){
                infile2 >> second_nn[i][j]; 
              }
          }

        infile2.close();
    }
    else{ 
        cout << "'neighbor_lists.txt' not found! \n"; 
        return 0; 
    }
 
    

    int sweep = N;          //number of steps for a complete lattice sweep
    
    
    
    //spin sites
    int * si = new int [N];
    initialize(N, si);
    
    /*
    for (int i = 0; i < N; i++){
        
        cout << si[i] << "\t";
             for (int j = 0; j < num_nn[0] ; j++)
                 cout << first_nn[i][j] << "\t";
             for (int j = 0; j < num_nn[1] ; j++)
                 cout << second_nn[i][j] << "\t";
             
             cout << "\n";
    }
    */
    
    //initial Configuration, Energy and Magnetization
    
    configuration( N, si, "init.txt" );
    double M = Magnetization(N, si);
    double H = Energy(N, si, num_nn, first_nn, second_nn, h, J);
    
   // cout << M << "\t" << H << "\n";
    
    
    
   // return 0;
    
    double ** si_hist = new double * [MCSteps+eqSteps];
    double ** auc = new double * [MCSteps+eqSteps];
    
    for (int i = 0; i < eqSteps+MCSteps; i++) 
    {
        si_hist[i] = new double[3];
        auc[i] = new double [2];
    }
    
    //averages, Cv, X
    double * data = new double  [10];
    
  
    
    //Average and Variance
    double MAvg, HAvg, MVar, HVar;
    double M2Avg, H2Avg, M2Var, H2Var;
    
    
    
    for (T = 0.05; T <=15.0; T+=0.5)
    {
        
        for (int i = 0; i < (eqSteps+MCSteps); i++)
        {   
            si_hist[i][0] = 0.0; si_hist[i][1] = 0.0; si_hist[i][2] = 0.0;
            auc[i][0] = 0.0; auc[i][1] = 0.0; 
        }
        
        
        for (int i =0; i < 10; i++) data[i] = 0.0;
        
        
        //Equilibration steps (burn-in)
        for (int step = 0; step < eqSteps; step++)
        {
            MC_move(N, si, num_nn, first_nn, second_nn, h, J, M, H, T);
        }
 
 
        ofstream con_file;
        std::ostringstream cf;
        cf <<"conf_"<<T<<".txt" ;
        std::string cnf = cf.str();
            

        //Production steps
        for (int step = 0; step < MCSteps; step++)
        {
            
            //Monte Carlo move over N*N*N atoms
            MC_move(N, si, num_nn, first_nn, second_nn, h, J, M, H, T);
            
            //averaging updated M and H values
            ave_var(step, MAvg, MVar, M);
            ave_var(step, HAvg, HVar, H);
            ave_var(step, M2Avg, M2Var, M*M);
            ave_var(step, H2Avg, H2Var, H*H);
            
            
            //cout << step << "\t" << T<<"\t" << M << "\t" << H << "\n";
            
            //collect autocorrelation snapshop
            //           aucf(si_hist, auc, M, H, nsave, step);
            
          if (step % nsave == 0){  
              configuration( N, si, cnf.c_str() );  
          }
            
        }//end loop over MCSteps
        
        //accumulate data
        data[0] = MAvg;  data[1] = MVar;
        data[2] = M2Avg; data[3] = M2Var;
        data[4] = HAvg;  data[5] = HVar;
        data[6] = H2Avg; data[7] = H2Var;
        data[8] = T;
        
        
        //calculate quantities and store
        save_data(N, data);
        
        //normalize and store autocorrelation function
        //        aucf_norm(auc, nsave, MAvg, HAvg, M2Avg, H2Avg, T, MCSteps);
        
        
        
   }
    
    
    free (si);
    free(num_nn);
    free(first_nn);
    free (second_nn);
    free (si_hist);
    free (data);
    free(auc);
    
    
    return 0;
}






