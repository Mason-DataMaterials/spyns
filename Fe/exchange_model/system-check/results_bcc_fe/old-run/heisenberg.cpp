#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <omp.h>

using namespace std;


random_device rd;
mt19937 engine(rd());
uniform_real_distribution<double> dist(0.0, 1.0);




//Random theta [0,pi]
double rng_theta(void)
{
    double t = acos( 1.0 - 2.0*dist(engine) );
    return t;
}

//Random phi [0,2*pi]
double rng_phi(void)
{
    double p = 2.0*M_PI*dist(engine);
    return p;
}

//Initialize the lattice/spin array
void initialize(int N, double ** lattice)
{
    int i,j,k; //Counters
    //Populate initialization lattice
    for (i = 0 ; i < N; i++)
    {
        
        lattice[i][0] = rng_theta();
        lattice[i][1] = rng_phi();
        //printf("Init latice: %f\t%f\n",lattice[i][0],lattice[i][1]);
    }
}


double magnetization(int N, double **si){
    
    //MAGNETIZATION
    //m_rms = M = (1/N)sqrt(sum(S_i*S_i))
    //NOTE: (Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
    double Mxsum = 0.0;
    double Mysum = 0.0;
    double Mzsum = 0.0;
    
    
    int i,j,k; //Counters
    //Populate initialization lattice
    for (i = 0 ; i < N; i++)
    {
        //Current angles
        double theta = si[i][0];
        double phi = si[i][1];
        
        double X = sin(theta)*cos(phi);
        double Y = sin(theta)*sin(phi);
        double Z = cos(theta);
        
        Mxsum += X;
        Mysum += Y;
        Mzsum += Z;
    }
    
    double Mx = Mxsum/(double)N;
    double My = Mysum/(double)N;
    double Mz = Mzsum/(double)N;
    double M = sqrt(Mx*Mx + My*My + Mz*Mz);
    
    return M;
    
}


double energy(int N, double ** si, double ** first_nn, double ** second_nn, int * num_nn, double * J){
    
    //ENERGY
    //Energy sweep - calculate energy of lattice in single sweep
    //Reset e_storage magnetization calc variables
   
    int i,j,k;
    double e_total = 0.0;
    double spinX, spinY,spinZ;
    double theta, phi, n_theta, n_phi, nn_theta, nn_phi;
    
    
       
    double * spinSum = new double [3];
    for (i =0; i<3; i++){
        spinSum[i] = 0.0;
    }
          
    
    for (i = 0 ; i < N; i++)
    {
        
        //Current angles
        theta = si[i][0];
        phi = si[i][1];
        //cout << i << "\t" << j << "\t" << k << "\t" << theta << "\t" << phi << "\n";
        
        //Left side of dot product
        spinX = sin(theta)*cos(phi);
        spinY = sin(theta)*sin(phi);
        spinZ = cos(theta);
        
        
        
        //Dot product calculation
        //E = -J (Si_current)*(S1 + S2+ S3+ S4+ S5+ S6)
        //(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
        
        for (int j = 0; j < num_nn[0]; j++)
        {
            n_theta  = si[(int)first_nn[i][j]][0];
            n_phi = si[(int)first_nn[i][j]][1];
            
            spinSum[0] = spinSum[0] + sin(n_theta)*cos(n_phi);
            spinSum[1] = spinSum[1] + sin(n_theta)*sin(n_phi);
            spinSum[2] = spinSum[2] + cos(n_theta);
        }
        
        double dotProduct_n = spinX*spinSum[0] + spinY*spinSum[1] + spinZ*spinSum[2];
        
        spinSum[0] = 0.0;
        spinSum[1] = 0.0;
        spinSum[2] = 0.0;
        
        for (int j = 0; j < num_nn[1]; j++)
        {
            nn_theta  = si[(int)second_nn[i][j]][0];
            nn_phi = si[(int)second_nn[i][j]][1];
            
            spinSum[0] = spinSum[0] + sin(nn_theta)*cos(nn_phi);
            spinSum[1] = spinSum[1] + sin(nn_theta)*sin(nn_phi);
            spinSum[2] = spinSum[2] + cos(nn_theta);
        }
        
        double dotProduct_nn = spinX*spinSum[0] + spinY*spinSum[1] + spinZ*spinSum[2];
        
        //Energy summation calculation
        e_total = e_total - J[0]*dotProduct_n - J[1]*dotProduct_nn;
        
    }
    
    
    //Calculate energy
    double e = e_total/2.0;
    
    free (spinSum);
    
    return e;
    
}

void mc_step(int N, double ** si, double ** first_nn, double ** second_nn, int * num_nn, double * J, 
             double temperature, double &accept, double &test){
    
    
    int i,j,k;
    double * currentSpin = new double [2];
    double * oldSpin = new double [2];
    double * newSpin = new double [2];
    for (i =0; i<2; i++){
        currentSpin[i] = 0.0;
        oldSpin[i] = 0.0;
        newSpin[i] = 0.0;
    }
    
    double theta, phi, thetaPrime, phiPrime, n_theta, n_phi, nn_theta, nn_phi;
    double * spinDiff = new double [3];
    double * spinSum = new double [3];
    for (i =0; i<3; i++){
        spinDiff[i] = 0.0;
        spinSum[i] = 0.0;
    }
    
    
    double beta = 1.0/temperature;
    
  //Sweeps
    for (int sweep = 0; sweep < N; sweep++)
    {
       //random spin site
        uniform_int_distribution<int> nsites_gen(0, N - 1);
        i = nsites_gen(engine);
        j = nsites_gen(engine);
        k = nsites_gen(engine); 
        
       
        //Current spin value
        theta = si[i][0];
        phi = si[i][1];
        
        //cout << s << "\t"<< i << "\t" << j << "\t" << k << "\t" << theta << "\t" << phi << "\n";
        //Propose new spin value
        thetaPrime  = rng_theta();
        phiPrime = phi;
        
    
        //checking theta, keep phiPrime the same
        
        //Left side of dot product
        spinDiff[0] = sin(thetaPrime)*cos(phiPrime) - sin(theta)*cos(phi);
        spinDiff[1] = sin(thetaPrime)*sin(phiPrime) - sin(theta)*sin(phi);
        spinDiff[2] = cos(thetaPrime) - cos(theta);
        
                
        for (int j = 0; j < num_nn[0]; j++)
        {
            n_theta  = si[(int)first_nn[i][j]][0];
            n_phi = si[(int)first_nn[i][j]][1];
            
            spinSum[0] = spinSum[0] + sin(n_theta)*cos(n_phi);
            spinSum[1] = spinSum[1] + sin(n_theta)*sin(n_phi);
            spinSum[2] = spinSum[2] + cos(n_theta);
        }
        
        double dotProduct_n = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
        
        spinSum[0] = 0.0;
        spinSum[1] = 0.0;
        spinSum[2] = 0.0;
        
        for (int j = 0; j < num_nn[1]; j++)
        {
            nn_theta  = si[(int)second_nn[i][j]][0];
            nn_phi = si[(int)second_nn[i][j]][1];
            
            spinSum[0] = spinSum[0] + sin(nn_theta)*cos(nn_phi);
            spinSum[1] = spinSum[1] + sin(nn_theta)*sin(nn_phi);
            spinSum[2] = spinSum[2] + cos(nn_theta);
        }
        
        double dotProduct_nn = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
        
        //Change in energy
        double deltaE = -J[0]*dotProduct_n - J[1]*dotProduct_nn;
        
        
        //Acceptance check
        test++;
        if (deltaE <= 0)
        {
            si[i][0] = thetaPrime;
            si[i][1] = phiPrime;
            accept++;
        }
        
        else //deltaE > 0
        {
            double prob = (sin(thetaPrime)/sin(theta))*exp(-beta*deltaE);
            double randNum = dist(engine);
            
            if (randNum <= prob)
            {
                si[i][0] = thetaPrime;
                si[i][1] = phiPrime;
                accept++;
            }
            
        }
        
        //Current spin value
        theta = si[i][0];
        phi = si[i][1];
        
        
       //checking phi, keep thetaPrime the same 
        
        //Propose new spin value
        thetaPrime  = theta;
        phiPrime = rng_phi();
        
        
        //Left side of dot product
        spinDiff[0] = sin(thetaPrime)*cos(phiPrime) - sin(theta)*cos(phi);
        spinDiff[1] = sin(thetaPrime)*sin(phiPrime) - sin(theta)*sin(phi);
        spinDiff[2] = cos(thetaPrime) - cos(theta);
   
    
                        
        for (int j = 0; j < num_nn[0]; j++)
        {
            n_theta  = si[(int)first_nn[i][j]][0];
            n_phi = si[(int)first_nn[i][j]][1];
            
            spinSum[0] = spinSum[0] + sin(n_theta)*cos(n_phi);
            spinSum[1] = spinSum[1] + sin(n_theta)*sin(n_phi);
            spinSum[2] = spinSum[2] + cos(n_theta);
        }
        
        dotProduct_n = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
        
        spinSum[0] = 0.0;
        spinSum[1] = 0.0;
        spinSum[2] = 0.0;
        
        for (int j = 0; j < num_nn[1]; j++)
        {
            nn_theta  = si[(int)second_nn[i][j]][0];
            nn_phi = si[(int)second_nn[i][j]][1];
            
            spinSum[0] = spinSum[0] + sin(nn_theta)*cos(nn_phi);
            spinSum[1] = spinSum[1] + sin(nn_theta)*sin(nn_phi);
            spinSum[2] = spinSum[2] + cos(nn_theta);
        }
        
        dotProduct_nn = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
        
        //Change in energy
        deltaE = -J[0]*dotProduct_n - J[1]*dotProduct_nn;
        

        //Acceptance check
        test++;
        if (deltaE <= 0)
        {
            si[i][0] = thetaPrime;
            si[i][1] = phiPrime;
            accept++;
        }
        
        else //deltaE > 0
        {
            double prob = (sin(thetaPrime)/sin(theta))*exp(-beta*deltaE);
            double randNum = dist(engine);
            
            if (randNum <= prob)
            {
                si[i][0] = thetaPrime;
                si[i][1] = phiPrime;
                accept++;
            }
        }
    }
    
    
    
    free (currentSpin);
    free (oldSpin);
    free (newSpin);
    free (spinDiff);
    free (spinSum);
    
    return;
}

//Temperature array initialization
void temp_init(double * temps)
{
    int t; //Counters
    double MAX_TEMP_STEPS = 13.0;
    double MAX_TEMP = 1300.0;
    //Temperature array initialization
    for (t = 0; t < MAX_TEMP_STEPS; t++)
    {
        double dT = MAX_TEMP/MAX_TEMP_STEPS;
        temps[t] = MAX_TEMP-t*dT;
    }
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

void save_data(int N, int MCSteps,double * data)
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
    
 /*   
    hfile << "T\t<M>\tM_Var\t<M^2>\tM^2_Var\t"
    << "<E>\tE_Var\t<E^2>\tE^2_Var\t"
    << "X\tC_v\n" ;
*/ 
    
    double MAvg = data[0], HAvg=data[4], MVar= data[1], HVar = data[5];
    double M2Avg = data[2], H2Avg=data[6], M2Var= data[3], H2Var=data[7];
    double M4Avg = data[10],  M4Var= data[11];
    
    
    double beta = 1.0/kT;
    double beta2 = beta*beta;
    double MN = 1.0/double (MCSteps*N);
    
    //Susceptibility
    double X = (beta*MN)*( M2Avg - (MAvg*MAvg));
    
    //Specific Heat
    double Cv = (beta2*MN)*(H2Avg - (HAvg*HAvg));
    
    double b_par = 1.0 - ( (1.0/3.0) * (M4Avg)/ (M2Avg*M2Avg) );
    //write to file
    hfile << T << "\t"
    << MAvg << "\t" << MVar << "\t" << M2Avg << "\t" << M2Var << "\t"
    << HAvg << "\t" << HVar << "\t" << H2Avg << "\t" << H2Var << "\t"
    << X    << "\t" << Cv   << "\t" 
    << M4Avg << "\t" << M4Var << "\t" << b_par << "\n" ;
    
    
    hfile.close();
    
    return;
}

void temp_step( int N, double ** si, double ** first_nn, double ** second_nn, int * num_nn, double * J, 
             double temperature, int MCSteps, int sampling_freq){
    
    double accept=0.0, test=0.0;
    
    double H = 0.0, M = 0.0;
    
    //Average and Variance
    double MAvg=0.0, HAvg=0.0, MVar=0.0, HVar=0.0;
    double M2Avg=0.0, H2Avg=0.0, M2Var=0.0, H2Var=0.0;
    double M4Avg=0.0, H4Avg=0.0, M4Var=0.0, H4Var=0.0;
    
    //hold final calculated quantities
    double * data = new double  [14];
    for (int i =0; i < 14; i++) data[i] = 0.0;
    
    int step = 0;
    int cc = 0;
    
    double kB = 8.617330350E-5; //eVK^-1
    double kT = 1000*kB*temperature; //meV
    
    cout << "Equilibrating for T = " << temperature << "\n";
    //equlibration
    for (step = 0; step < 25000; step++){
        
        mc_step( N, si, first_nn, second_nn, num_nn, J, kT, accept, test );
    }
    
    cout << "Calculating for T = " << temperature << "\n";
    //measurement
    for (step = 0; step < MCSteps; step++){
        
       mc_step( N, si, first_nn, second_nn, num_nn, J, kT, accept, test );
        
        
       //cout << step <<"\t" <<  H << "\t" << M << "\t" << accept/test << "\n";
        
        if (step % sampling_freq == 0 )
        { 
            M = magnetization (N, si);
            H = energy( N, si, first_nn, second_nn, num_nn, J); 
            
            ave_var(cc, MAvg, MVar, M);
            ave_var(cc, HAvg, HVar, H);
            ave_var(cc, M2Avg, M2Var, M*M);
            ave_var(cc, H2Avg, H2Var, H*H);
            ave_var(cc, M4Avg, M4Var, M*M*M*M);
            
           //cout << cc <<"\t" <<  HAvg << "\t" << MAvg  << "\n";
            cc++;
        }
        
        
    }
    
    //accumulate data
    data[0] = MAvg; data[1] = MVar;
    data[2] = M2Avg; data[3] = M2Var;
    data[4] = HAvg; data[5] = HVar;
    data[6] = H2Avg; data[7] = H2Var;
    data[8] = temperature; data[9] = kT; 
    data[10] = M4Avg; data[11] = M4Var;
    
    
    save_data(N,(int) ((double)step/(double)sampling_freq),data);
    
    free(data);
    
    
    return;
}


int main(int argc, char *argv[]){
    
    int N;
    double * J = new double [2];
    
    J[0] = 13.1; //meV
    J[1] = 13.7; //meV
    
    double temperature = 0.1;
    int MCSteps = 75000;
    int sampling_freq = 50;
    
    
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
        
        infile2.close();
    }
    else{ 
        cout << "'neighbor_lists.txt' not found! \n";
        return 0; 
    }
    
    
    //initialize spin configuration
    double ** si = new double * [N];
    for (int i = 0; i < N; i++){
        si[i] = new double [2];
    }
    
    initialize(N, si);
    
    double M = magnetization (N, si);
    double H = energy( N, si, first_nn, second_nn, num_nn, J); 
    
    ofstream ifile; ifile.open("init.txt");
    ifile << N << "\n";
    
    //print initial config to file
    for (int i = 0 ; i < N; i++)
    {
        
        //Current angles
        double theta = si[i][0];
        double phi = si[i][1];
        ifile << i << "\t" << theta << "\t" << phi << "\n";
        
    }
    ifile.close();
    
    double * temps = new double [15];
    temp_init(temps);
    

    for (int t = 0; t < 13; t++){
        
        temperature = temps[t];
        temp_step( N, si, first_nn, second_nn, num_nn, J, temperature, MCSteps, sampling_freq );
        
        cout << "done T = " << temperature << "\n";
    }

    
    
    free(si);
    
    return 0;
    
}
