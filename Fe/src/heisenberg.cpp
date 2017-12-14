
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



/* initialize(N, L, a, si)
 * calculates the box side length and initializes
 * particles arranged in a bcc lattice structure
 * with randomized spins in the x, y and z directions
 * within a box of side L.
 *
 * takes in: N - number of particles
 *           L - length of box side
 *           a - lattice constant
 *           si - coordinates and x,y,z spin values 
 *                                for each particle
 *
 */


void initialize(int N, double& L, double a, double ** si ){
    
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
                        
                        
                        //(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
                        double theta = 1.0/cos(1.0-2.0 * (double) dist(engine));
                        double phi = M_PI*2.0* (double)dist(engine);
                        
                        //theta
                        si[n][3] = theta;
                        
                        //phi
                        si[n][4] = phi;
                        
                        
                        //Sx
                        si[n][5] = sin( si[n][3] )*cos( si[n][4] );
                        //Sy
                        si[n][6] = sin( si[n][3] )*sin( si[n][4] );
                        //Sz
                        si[n][7] = cos( si[n][3] );
                        
                        ++n;
                    }
                }
        }
    }
    
    /*
    cout << N << "\n";
    cout << " " << "\n";
    
    for (int i = 0; i< N; i++){
    
        cout << "Fe" << "\t" << si[i][0] << "\t" << si[i][1] << "\t" << si[i][2] << "\t"
                             << si[i][5] << "\t" << si[i][6] << "\t" << si[i][7] << "\n";
    }
    */
    
    
    
    return;
}

/* Magnetization(N, si)
 * Calculates total magnetization of the system. 
 * takes in: N - number of particles
 *           si - coordinates and x,y,z spin values 
 *                                for each particle
 *
 */

double Magnetization(int N, double ** si){
    
    double M = 0.0;
    
    double * m = new double [3];
    for(int i = 0; i < 3; i++) m[i] = 0.0;
    
    for (int i = 0; i < N; i++)
    {
        int k = 0;
        for(int p = 5; p <8; p++)
        {
            m[k] += si[i][p];
        }
    }
    
    for (int i = 0; i < 3; i++) m[i] /= (double) N;
    M = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
    
    return M;
    
}


/* 
 * bcc_neighbor_lists(N, L, a, si, first_nn, second_nn)
 * records the first (nn1 = 8) and second (nn2 = 6) 
 * neighbors of every point in the bcc grid.
 * takes in: N - number of particles
 *           L - length of box side
 *           a - lattice constant
 *           si - coordinates and x,y,z spin values 
 *                                for each particle
 *           first_nn - empty array for lists of
 *                        first nearest neighbors
 *           second_nn - empty array for lists of
 *                         second nearest neighbors
 *
 */

void bcc_neighbor_lists(int N, double L, double a, double ** si, double ** first_nn, double ** second_nn)
{
    
    double nn_bcc_distances[14] = {0.0, sqrt(3.0)/2.0, 1.0, sqrt(2.0), sqrt(11.0/4.0),
    sqrt(3.0), 2.0, sqrt(19.0/4.0), sqrt(5.0), sqrt(6.0),
    sqrt(27.0/4.0), sqrt(8.0), sqrt(35.0/4.0), 3.0};


    for (int p = 0; p < N; p++)
       { 
            int nn1_count = 0;
            int nn2_count = 0;
            
            for (int i = 0; i < N; i++)
	    {
                                
                double rSqd =0.0;
                double dr[3] = {0.0};
                double dist =0.0;
                
                dr[0] = si[p][0]-si[i][0] ;
                dr[1] = si[p][1]-si[i][1] ;
                dr[2] = si[p][2]-si[i][2] ;
                
                
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
               
                if (p != i)
		{
			
			double diff1 = round(dist-nn_bcc_distances[1]*a) ;
                
               		if ( diff1 == 0.0)//1e-4 && dist > nn_bcc_distances[0]*a  )
                	{
                   
                   		first_nn[p][nn1_count] = i;
                    		nn1_count+=1;
                    
                	}
                
                	double diff2 = round(dist-nn_bcc_distances[2]*a) ;
                 
                	if ( diff2 == 0.0) //1e-4 && dist > nn_bcc_distances[1]*a )
                	{
                    		second_nn[p][nn2_count] = i;
                    		nn2_count+=1;
                         
                	} 
		}
               
            }//end for i
           
       }//end for p
 
 return;    
 
}


/* neighbors(N, si, point)
 * calculates the sum of spins for the
 * first (nn1 = 8) and second (nn2 = 6) 
 * neighbor sites of a point in 
 * the bcc grid.
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           point - site of interest
 * returns: sj - sum of the spins of neighbors of point
 */
void neighbors(int N, int point, int xyz, double * sj, double ** si, double ** first_nn, double ** second_nn)
{    
    
   
    
    for (int i =0; i<8; i++){
        int index = int (first_nn[point][i]);
        sj[0] += si[index][xyz];
    }
        
    for (int i =0; i<6; i++){
        int index = int (second_nn[point][i]);
        sj[1] += si[index][xyz];
    }
        
   
    return;
}

/* Energy(N, si, first_nn, second_nn, h, J)
 * Calculates total energy of the system. 
 * takes in: N - number of particles
 *           si - coordinates and x,y,z spin values 
 *                                for each particle
 *           first_nn - empty array for lists of
 *                        first nearest neighbors
 *           second_nn - empty array for lists of
 *                         second nearest neighbors
 * returns: energy, H
 */


double Energy(int N, double ** si, double ** first_nn, double ** second_nn, double h, double * J)
{
    double H  = 0.0;
    
    double * Hi = new double [3];
    for (int i =0; i <3; i++) Hi[i] = 0.0;
    
    
    double * sj_x = new double [2];
    sj_x[0] = sj_x[1] = 0;
    
    double * sj_y = new double [2];
    sj_y[0] = sj_y[1] = 0;
    
    double * sj_z = new double [2];
    sj_z[0] = sj_z[1] = 0;
    
    
    for (int i = 0; i < N; i++)
    {
        
        int point = i;
        
        neighbors(N, point, 5, sj_x, si, first_nn, second_nn );
        neighbors(N, point, 6, sj_y, si, first_nn, second_nn );
        neighbors(N, point, 7, sj_z, si, first_nn, second_nn );
         
        Hi[0] = Hi[0] + ( (si[i][5]*-(J[0]*sj_x[0] + J[1]*sj_x[1]) ) - h*si[i][5] );
        Hi[1] = Hi[1] + ( (si[i][6]*-(J[0]*sj_y[0] + J[1]*sj_y[1]) ) - h*si[i][6] );
        Hi[2] = Hi[2] + ( (si[i][7]*-(J[0]*sj_z[0] + J[1]*sj_z[1]) ) - h*si[i][7] );
        
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

double dE_theta (       int N, 
                        double ** si, 
                        double ** first_nn, 
                        double ** second_nn, 
                        double h, 
                        double * J, 
                        int point, 
                        double new_theta 
                )
{
    double * sj_x = new double [2];
    sj_x[0] = sj_x[1] = 0;
    
    double * sj_y = new double [2];
    sj_y[0] = sj_y[1] = 0;
    
    double * sj_z = new double [2];
    sj_z[0] = sj_z[1] = 0;
    
    neighbors(N, point, 5, sj_x, si, first_nn, second_nn );
    neighbors(N, point, 6, sj_y, si, first_nn, second_nn );
    neighbors(N, point, 7, sj_z, si, first_nn, second_nn );
    
    double theta = si[point][3];
    double phi   = si[point][4];
    
    
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

double dE_phi (         int N, 
                        double ** si, 
                        double ** first_nn, 
                        double ** second_nn, 
                        double h, 
                        double * J, 
                        int point, 
                        double new_phi
              )
{
    double * sj_x = new double [2];
    sj_x[0] = sj_x[1] = 0;
    
    double * sj_y = new double [2];
    sj_y[0] = sj_y[1] = 0;
    
    double * sj_z = new double [2];
    sj_z[0] = sj_z[1] = 0;
    
    neighbors(N, point, 5, sj_x, si, first_nn, second_nn );
    neighbors(N, point, 6, sj_y, si, first_nn, second_nn );
    neighbors(N, point, 7, sj_z, si, first_nn, second_nn );
    
    double theta = si[point][3];
    double phi   = si[point][4];
    
    
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

 
 /* MC_move(N, si, first_nn, second_nn, h, J, M, H, T)
 *  carries out the monte carlo  move at each time step
 *  to update the xyz spin values for all atoms and
 *  update the energy, H, and magnetization, M, values.
 * takes in: N - length of one side
 *           si - 3D lattice of spin sites
 *           first_nn - empty array for lists of
 *                        first nearest neighbors
 *           second_nn - empty array for lists of
 *                         second nearest neighbors
 *           h - external manetic field
 *           J - nearest neighbor coupling constant
 *           M - magnetization
 *           H - energy
 *           T - Temperature multiplied by the Boltzmann const.
 */


void MC_move(int N, double ** si,  double ** first_nn, double ** second_nn, 
                            double h, double *J, double &M, double &H, double T)
{
    double phi, theta, p;
    
    double beta = 1.0/T;
    int * point = new int [3];
    
    for (int sweep = 0; sweep < N; sweep++)
    {
        //random spin site
        uniform_int_distribution<int> nsites_gen(0, N - 1);
        int point = nsites_gen(engine);
        
       
        
        //theta update, phi constant
        theta = 1.0/cos(1.0-2.0 * (double)dist(engine));
        
        double dE = dE_theta (N, si, first_nn, second_nn, h, J, point, theta );
    
        int i = sweep;
        
        if (dE <= 0.0)
        {
                //theta
                si[i][3] = theta;
                //Sx
                si[i][5] = sin(si[i][3])*cos(si[i][4]);
                //Sy
                si[i][6] = sin(si[i][3])*sin(si[i][4]);
                //Sz
                si[i][7] = cos(si[i][3]);
            
                H = H+dE;
        }
        else {
            
            p = dist(engine);
            
            if (p <= exp(-dE*beta))
            {
                
                si[i][3] = theta;
                //Sx
                si[i][5] = sin(si[i][3])*cos(si[i][4]);
                //Sy
                si[i][6] = sin(si[i][3])*sin(si[i][4]);
                //Sz
                si[i][7] = cos(si[i][3]);
            
                H = H+dE;
                
            }
            
            
        } //end if for theta
         
        
        //phi update, theta constant
        phi = M_PI*2.0* (double)dist(engine);
        dE = dE_phi(N, si, first_nn, second_nn, h, J, point, phi);
        
        
        if (dE <= 0.0)
        {
            
                //phi
                si[i][4] = phi;
                //Sx
                si[i][5] = sin(si[i][3])*cos(si[i][4]);
                //Sy
                si[i][6] = sin(si[i][3])*sin(si[i][4]);
               
            
                
        }
        else {
            p = dist(engine);
            
            if (p < exp(-dE*beta))
            {
                //phi
                si[i][4] = phi;
                //Sx
                si[i][5] = sin(si[i][3])*cos(si[i][4]);
                //Sy
                si[i][6] = sin(si[i][3])*sin(si[i][4]);
                
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
    double kB = data[9];
    double kT = kB*T;
    
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
    
    
    double beta = 1.0/kT;
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

void aucf(double * si_hist, double MAvg, int MCSteps, int T)
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
     int N = 1024; 
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
    
    acfile.close();
    
    
 return;
    
}
/* plot_script (T_min, T_max, t_steps)
 * Generates script for visualization with gnuplot
 * takes in: T_min - minimum simulation temperature
 *           T_max - maximum simulaiton temperature
 *           t_steps - step size between temperatures
 *                      during simulation.
*/                       
void plot_script(double T_min, double T_max, double t_steps)
{
    ofstream gfile; gfile.open("aucfplot.gnu", ios_base::app);
    
    double dT = (T_max - T_min)/ t_steps;
    
    int t_count = 1;
    for (double T = T_min; T <=T_max; T+=dT)
    {
        if (T == T_min)
            gfile << "plot 'aucf" <<t_count<<".txt' u 1:3 w l t 'T="<<T<<"',\\"<<"\n";
        else if (T == T_max)
            gfile << "     'aucf" <<t_count<<".txt' u 1:3 w l t 'T="<<T<<"'" <<"\n";
        else 
            gfile << "     'aucf" <<t_count<<".txt' u 1:3 w l t 'T="<<T<<"',\\"<<"\n";
            
    }
    
}


int main(int argc, char *argv[]){

    
    //lattice dimension, lenght of single side
     int N = 250; 
    
    //external magnetic field
    double h = 0.0;                     
    
    //coupling constant (Ry)
    double * J = new double [2];
    J[0] = 7.52373031e-03;   J[1] = 2.55273958e-03;
    
    int MCSteps = 50000;                //number of Monte Carlo moves
    int eqSteps = 0.2*MCSteps;          //number of equilibration steps
    
    
    double T;  
    
    //boltzmann constant(Ry K^-1)
    double kB = 0.0000063306;
    
    
    //lattice constant for Fe
    double a = 5.3970578; //a.u.
    double L;
    
    
    double **si = new double * [N];//spin configuration
    for (int i = 0; i < N; i++)
        si[i] = new double [8];
    
      
    initialize(N, L, a, si );

     
    double **first_nn = new double * [N];
    double **second_nn = new double * [N];
     for (int i = 0; i < N; i++)
        {
            first_nn[i] = new double [8];
            second_nn[i] = new double [6];
        }
       
    bcc_neighbor_lists(N, L, a, si, first_nn, second_nn);
    
    
    double * si_hist = new double [MCSteps+eqSteps];
        
    for (int i = 0; i < eqSteps+MCSteps; i++) 
    {
        si_hist[i] = 0.0;
    }
    
    
    double * data = new double  [10];
    for (int i =0; i < 10; i++) data[i] = 0.0;
    

    double M = Magnetization(N, si);
   
    double H = Energy(N, si, first_nn, second_nn, h, J);
    
     
     //Average and Variance
    double MAvg, HAvg, MVar, HVar;
    double M2Avg, H2Avg, M2Var, H2Var;
    
    int t_count = 1;
    
    double T_min = 400.0;
    double T_max = 1200.0;
    double t_steps = 50.0;
    double dT = (T_max - T_min) / t_steps;
    
    //loop over T values
    for (T = T_min; T <= T_max; T+=dT)
    {
        
        double kT = kB*T;
        
        //Equilibration steps (burn-in)
        for (int step = 0; step < eqSteps; step++)
        {
            MC_move(N, si, first_nn, second_nn, h, J, M, H, kT);
            //cout << step << "\t" <<M <<"\t" << MAvg << "\t" << HAvg << "\n";
        }
        
        
        //Production steps
        for (int step = 0; step < MCSteps; step++)
        {
            /*each move is a complete lattice sweep - N*N*N steps */
            MC_move(N, si, first_nn, second_nn, h, J, M, H, kT);
            
            M = Magnetization(N, si);
            H = Energy(N, si, first_nn, second_nn, h, J);
            
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
        data[8] = T; data[9] = kB;
        
        cout << "saving to file ... \n";
        save_data(N, data);
        
        aucf(si_hist, MAvg, MCSteps, t_count);
        
        
        cout << "Done : T = " << T << "\t M = " << M << "\t H = " << H << "\n";
        
        t_count+=1;
        
    }//end loop over T
    
    

    plot_script(T_min, T_max, t_steps);
        
    
    return 0;
}
