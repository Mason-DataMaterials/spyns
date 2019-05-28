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


//temperature values array initialization
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

//Initialize the lattice/spin array
void initialize(int N, double ** si)
{
    int i,j,k; //Counters
    double theta, phi;
    //Populate initialization lattice
    for (i = 0 ; i < N; i++)
    {
        
        theta = rng_theta();
        phi = rng_phi();
        si[i][0] = theta;
        si[i][1] = phi;
        si[i][2] = sin(theta)*cos(phi);
        si[i][3] = sin(theta)*sin(phi);
        si[i][4] = cos(theta);
       // printf("Initial lattice: %f\t%f\n",si[i][0],si[i][1]);
    }
    
    return;
}


double magnetization(int N, double **si){
    
    
    double Mxsum = 0.0;
    double Mysum = 0.0;
    double Mzsum = 0.0;
    
    
    int i,j,k; //Counters
    //Populate initialization lattice
    for (i = 0 ; i < N; i++)
    {
        //Current angles
        
        double spinX = si[i][2];
        double spinY = si[i][3];
        double spinZ = si[i][4];
        
        
        Mxsum += spinX;
        Mysum += spinY;
        Mzsum += spinZ;
    }
    
    double Mx = Mxsum/(double)N;
    double My = Mysum/(double)N;
    double Mz = Mzsum/(double)N;
    
    double M = sqrt(Mx*Mx + My*My + Mz*Mz);
    
    return M;
    
}


double site_energy(int site, double ** si, double ** neighbor_list, int num_neighbors, double * J){
    
    int i,j,k;
    double e_site = 0.0;
    double spinX, spinY,spinZ;
    
    
    double * spinSum = new double [3];
    for (i =0; i<3; i++){
        spinSum[i] = 0.0;
    }
    
    
    //Dot product calculation
    //E = -J (Si_current)*(S1 + S2 + S3 + S4 + S5 + S6)
    
    //lookup neighbors and sum
    for (int j = 0; j < num_neighbors; j++)
    {
        //sum ( Sj_k * J ) , k=x,y,z
        spinSum[0] = spinSum[0] + si[(int)neighbor_list[site][j]][2]*J[j]; //J * Sj_x
        spinSum[1] = spinSum[1] + si[(int)neighbor_list[site][j]][3]*J[j]; //J * Sj_y
        spinSum[2] = spinSum[2] + si[(int)neighbor_list[site][j]][4]*J[j]; //J * Sj_z
    }
    
    spinX = si[site][2];
    spinY = si[site][3];
    spinZ = si[site][4];
    
    //( Si_k * Sj_k * J )
    double dotProduct_J = spinX*spinSum[0] + spinY*spinSum[1] + spinZ*spinSum[2];
    
    //Energy summation calculation: sum(-Jj * Si_k * Sj_k) over j 
    e_site = - dotProduct_J;
    
    free (spinSum);
    
    return e_site;
    
}

double total_energy( int N, double ** si, double ** neighbor_list, int num_neighbors, double * J){
    
    double e_total = 0.0;
    
    for (int i = 0; i < N; i++){
        
       e_total += site_energy(i, si, neighbor_list, num_neighbors,  J);
    }
    
    return e_total/2.0;
        
}

void get_trial_spin_for_site(double * trial_spin){
    
        double theta = rng_theta();
        double phi = rng_phi();
        trial_spin[0] = theta;
        trial_spin[1] = phi;
        trial_spin[2] = sin(theta)*cos(phi);
        trial_spin[3] = sin(theta)*sin(phi);
        trial_spin[4] = cos(theta);
        
        return;
}


double trial_spin_energy_at_site(int site, double ** si, double ** neighbor_list, int num_neighbors, double * J, double * trial_spin ){
    
    int i,j,k;
    double e_trial = 0.0;
    double spinX, spinY,spinZ;
    
    
    double * spinSum = new double [3];
    
    for (i =0; i<3; i++){
        spinSum[i] = 0.0;
    }
    
    
    //lookup neighbors and sum
    for (int j = 0; j < num_neighbors; j++)
    {
        //sum ( Sj_k * J ) , k=x,y,z
        spinSum[0] = spinSum[0] + si[(int)neighbor_list[site][j]][2]*J[j]; //J * Sj_x
        spinSum[1] = spinSum[1] + si[(int)neighbor_list[site][j]][3]*J[j]; //J * Sj_y
        spinSum[2] = spinSum[2] + si[(int)neighbor_list[site][j]][4]*J[j]; //J * Sj_z
    }
    
    spinX = trial_spin[2];
    spinY = trial_spin[3];
    spinZ = trial_spin[4];
    
    //( Si_k * Sj_k * J )
    double dotProduct_J = spinX*spinSum[0] + spinY*spinSum[1] + spinZ*spinSum[2];
    
    //Energy summation calculation: sum(-Jj * Si_k * Sj_k) over j 
    e_trial = - dotProduct_J;
    
    free (spinSum);
    
    return e_trial;
}


bool trial_flip( double &deltaE, double current_site_energy, double trial_energy, double beta ){
    
    deltaE = trial_energy - current_site_energy;
    
    if (deltaE <= 0)
        {
           return true;
        }
        
        else //deltaE > 0
        {
            double prob = exp(-beta*deltaE);
            double randNum = dist(engine);
            
            if (randNum < prob)
            {
                return true;
            }
        }
                
        return false;
}



void mc_step(int N, double ** si, int num_neighbors, double ** neighbor_list, double * J,  double temperature){
    
    uniform_int_distribution<int> nsites_gen(0, N - 1);
    
    double beta = 1.0/temperature;
    
    double * trial_spin = new double [5];
    
    double deltaE;
    double current_site_energy;
    double trial_spin_energy;
    
    bool keep_trial_spin = false;
    
    //Sweeps
    for (int sweep = 0; sweep < N; sweep++)
    {
       //random spin site
        int i = nsites_gen(engine);
        
        
        for (int j = 0; j < 5; j++){
            trial_spin[j] = 0.0;
        }
        
        get_trial_spin_for_site(trial_spin);

        current_site_energy = site_energy( i, si, neighbor_list, num_neighbors, J );
        trial_spin_energy = trial_spin_energy_at_site( i, si, neighbor_list, num_neighbors,J, trial_spin );
    
        deltaE = 0.0;
        keep_trial_spin = trial_flip (deltaE, current_site_energy, trial_spin_energy, beta);
        
        //update is trial spin is accepted
        if (keep_trial_spin){
            
            si[i][0] = trial_spin[0]; //theta
            si[i][1] = trial_spin[1]; //phi
            si[i][2] = trial_spin[2]; //sin(theta)*cos(phi);
            si[i][3] = trial_spin[3]; //sin(theta)*sin(phi);
            si[i][4] = trial_spin[4]; //cos(theta);
        }
        
    }//end sweep
    
    
    free (trial_spin);
    
    return;
    
}


void ave_var(int cc, double &Avg, double &Var, double U)
{
    
    double old_Avg = Avg;
    
    Avg = cc == 0 ? U : Avg + ((U - Avg) / cc);
    
    Var = cc == 0 ? 0 : Var * cc + (U - old_Avg) * (U - Avg);
    Var /= (cc + 1);
}

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
    double X = (beta/MN)*( M2Avg - (MAvg*MAvg));
    
    //Specific Heat
    double Cv = (beta2*MN)*(H2Avg - (HAvg*HAvg));
    
    double b_par = 1.0 - ( (1.0/3.0) * (M4Avg)/ (M2Avg*M2Avg) );
    //write to file
    hfile << T << "\t" << beta << "\t"
    << MAvg << "\t" << MVar << "\t" << M2Avg << "\t" << M2Var << "\t"
    << HAvg << "\t" << HVar << "\t" << H2Avg << "\t" << H2Var << "\t"
    << X    << "\t" << Cv   << "\t" 
    << M4Avg << "\t" << M4Var << "\t" << b_par << "\n" ;
    
    
    hfile.close();
    
    return;
}


void temp_step( int N, double ** si, int num_neighbors, double ** neighbor_list, double *J,  
                double temperature, int MCSteps, int sampling_freq){
    
    
    //total energy magnetization
    double H = 0.0 , M = 0.0;
    
    double kB = 8.617330350E-5; //eVK^-1
    double kT = 1000.0*kB*temperature; //meV
    //double kT = temperature;
    
    int step;
    
    //equilibration
    for ( step = 0; step < 50000; step++){
        
        mc_step( N, si, num_neighbors, neighbor_list, J,  kT);
    }
    
    
    //holder for final calculated quantities
    double * data = new double  [14];
    for (int i =0; i < 14; i++) data[i] = 0.0;
    
    //Average and Variance
    double MAvg=0.0, HAvg=0.0, MVar=0.0, HVar=0.0;
    double M2Avg=0.0, H2Avg=0.0, M2Var=0.0, H2Var=0.0;
    double M4Avg=0.0, H4Avg=0.0, M4Var=0.0, H4Var=0.0;
    
    
    //production steps
    int cc =0;
    for ( step = 0; step < MCSteps; step++){
        

        mc_step( N, si, num_neighbors, neighbor_list, J,  kT);
        
         if (step % sampling_freq == 0 )
        { 
           
            H = total_energy( N, si, neighbor_list, num_neighbors, J);
            M = magnetization (N, si);
            
            ave_var(cc, MAvg, MVar, M);
            ave_var(cc, HAvg, HVar, H);
            ave_var(cc, M2Avg, M2Var, M*M);
            ave_var(cc, H2Avg, H2Var, H*H);
            ave_var(cc, M4Avg, M4Var, M*M*M*M);
            //cout << cc <<"\t" <<  H << "\t" << HAvg  << "\n";
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
    
    
    //print to file
    save_data( N, (int) ((double)step/(double)sampling_freq),data);
    
    
    free(data);
    
    
    return;    
}

void save_snapshot(int N, double ** si, double temperature ){

    /*open file for writing*/
    ofstream hfile;
    std::ostringstream hf;
    hf <<"snapshots/config_"<<temperature<<".txt" ;
    std::string hmf = hf.str();
    hfile.open(hmf.c_str(), ios_base::app);
    
    
    double theta, phi, spinX, spinY, spinZ;
     
    for (int i = 0; i < N; i++ ){
      theta = si[i][0];
      phi   = si[i][1];
      
      spinX = sin(theta)*cos(phi);
      spinY = sin(theta)*sin(phi);
      spinZ = cos(theta);
      
      //write to file
      hfile << i  <<"\t" << theta << "\t" << phi << "\t" << spinX << "\t" << spinY << "\t" << spinZ << "\n";
          
    }
      
    hfile.close();   
    
}

int main(int argc, char *argv[]){

    
    int N;
    int num_neighbors;
    double ** neighbor_list;
    double * J;
    
    int MCSteps = 100000;
    int sampling_freq = 50;
    
    //double J1 = 1.0;
    //double J2 = 0.0;
    double J1 = 13.1; //meV
    double J2 = 13.7; //meV
    J1 = J1*2.0;
    J2 = J2*2.0;
    
    //read in tables and build lists
    int num_1st_nn,num_2nd_nn;
    ifstream infile2 ("neighbor_lists.txt");
    
    if (infile2.is_open())
    {
        
        infile2>>N;
        infile2>>num_1st_nn>>num_2nd_nn;
        
        num_neighbors = num_1st_nn + num_2nd_nn;
        
        //populate neighbor lists
        neighbor_list =  new double * [N];
        
        for (int i = 0; i < N; i++)
        {
            neighbor_list[i] = new double [num_neighbors];
        }  
        
        for (int i = 0; i < N; i++){
            for (int j = 0; j < num_neighbors; j++){
                infile2 >> neighbor_list[i][j]; 
            }
        }
        
        infile2.close();
    }
    else{ 
        cout << "'neighbor_lists.txt' not found! \n";
        return 0; 
    }
    
    
    J = new double [num_neighbors];
    
    for (int i = 0; i < num_neighbors; i++){
        if (i<num_1st_nn) J[i] = J1;
        else if (i >= num_1st_nn || i < num_2nd_nn)
            J[i] = J2;
        
    }
    
    double **si = new double* [N];
    for (int i = 0; i < N; i++){
        si[i] = new double [5];
    }
    
   
    //random initialization
    initialize( N,  si);

   
    
    //temperature steps
    //1-initialize list of temperature values:
    double * temps = new double [15];
    temp_init(temps);
    
    
    //2-loop over temperature values
    for (int i = 3; i < 7; i++){
        
        double temperature = temps[i];
        temp_step( N, si, num_neighbors, neighbor_list, J,temperature, MCSteps, sampling_freq);
        
        //save configuration
        save_snapshot(N, si, temperature );
    }
    
    
    
       free(si);
       free (temps);
       free(neighbor_list);
       free(J);
        

    return 0;
}
