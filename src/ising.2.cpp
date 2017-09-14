//ising 2.0
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>


using namespace std;

double rnd (){
    
    /*http://stackoverflow.com/a/22935167*/
    random_device rd{};    
    mt19937 engine{rd()};
    uniform_real_distribution<double> dist{0.0, 1.0};
    double x = dist(engine);
    return x;
}

void initialize(int N, int *** si){
    
    for (int x = 0; x < N; x++){
        si[x] = new int * [N];
        
        for (int y = 0; y < N; y++){
            si[x][y] = new int [N];
            
            for (int z = 0; z < N; z++){
                int p  = (int) rnd();
                if (p > 0.5) si[x][y][z] = 1;
                else si[x][y][z] = -1;
            }
        }
    }
    
    return;
}

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


double Magnetization(int N, int *** si){
    
    double M = 0.0;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                M = M+(double)si[i][j][k];
            
            return M;
        
}

double Energy( int N, int *** si, double h, double J){
    
    double H = 0.0;
    int * point = new int [3];
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                point[0] = i; point[1] = j; point[2]=k;
                
                int sj = neighbors(N, si, point);
                H = H - ((J*si[i][j][k]*(sj)) + h*si[i][j][k]);
            }
            
            free (point);
        
        return H / 4.;
}

double deltaE (int N, int *** si, double h, double J, int * point){
    
    int sj = neighbors(N, si , point);
    double dE = -2.0 *(J*si[point[0]][point[1]][point[2]]*sj 
                                + h*si[point[0]][point[1]][point[2]]);
    
    return dE;    
}


void MC_move(int N, int *** si, double h, double J, double &M, double &H, double T){
    
    double beta = 1.0/T;
    int * point = new int [3];
    int i = (int) rnd()%N; int j = (int)rnd()%N; int k = (int) rnd()%N;
    
    point[0] = i; point[1] = j; point[2]=k;
    
    
    double dE = deltaE(N, si, h, J, point);
    double dm=0.0,dh=0.0;
    double p = T*log(rnd());
    if (dE > p)
    {
        
        si[i][j][k] *= -1;
        dm=2.0*si[i][j][k];
        dh=dE;
    }
    
    M = M+dm;
    H = H+dh;
    
    free (point);
    return;
}

void collect_data(double * data, double M, double H){
    data[0] = data[0]+M;
    data[1] = data[1]+H;
    data[2] = data[2]+M*M;
    data[3] = data[3]+H*H;
    data[4] = data[4]+fabs(M);
}

void configuration(int N, int *** si, string filename){
    
    ofstream cfile;
    cfile.open(filename.c_str());
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                cfile<<i <<"\t" << j << "\t" << k << "\t" << si[i][j][k] << "\n";
            }
            
            cfile.close();
        return;
    
}

int main(){
   
    
  int N = 10;
  int MCSteps = 300*N*N*N;
  double h = 1.0;
  double J = 2.0;
  
  int *** si = new int ** [N];
  initialize(N, si);
  
  double T;
  double M = Magnetization(N, si);
  double H = Energy(N, si, h, J);
  
  // do{
        cout << "Running h = " << h << "\n";
        ofstream tfile;
        std::ostringstream tf;
        tf << h <<"_mg.txt" ;
        std::string tmf = tf.str();
        tfile.open(tmf.c_str());
        
        
        for (T=1.0; T <= 10.; T+=1.){
            
            
            double * data;
            data = new double [5];
            for (int i=0; i < 5; i++) data[i] = 0.0;
            
            
            //Equilibration Steps
            for (int step = 0; step < MCSteps; step++)
                MC_move(N, si, h, J, M, H, T);
            
            //Production
            for (int step = 0; step < 0.5*MCSteps; step++){
                
                MC_move(N, si, h, J, M, H, T);
                
                //recording and updating the M and H values 
                //Frequency can be varied depending on length/size of simulation
                collect_data(data, M, H);
            }
            
            //Normalization
            double norm   = 1.0/(MCSteps*N*N*N);
            double M_avg  = data[0]*norm;
            double H_avg  = data[1]*norm;
            double M_2avg = data[2]*norm;
            double H_2avg = data[3]*norm;
            double M_abs  = data[4]*norm;
            
            
            //Quantities
            
            double beta = 1.0/T;
            double beta2 = beta*beta;
            
            //Susceptibility
            double X = beta*( M_2avg - (M_avg*M_avg));
            //Specific Heat
            double Cv = beta2*(H_2avg - (H_avg*H_avg));
            
            tfile << T << "\t" << M_avg << "\t" << H_avg << "\t" << X << "\t" << Cv <<"\n";
            
           /* 
            //save configuration after every T sweep. 
            //Frequency can be varied for longer/larger simulations
            std::ostringstream sfile;
            sfile <<"config_"<<h<<"_"<<T <<".txt" ;
            std::string sfilename = sfile.str();
            configuration(N, si, sfilename.c_str());
           */
           
            free(data);
            cout << "Done h = " << h << ", T = " << T <<"\n";

        }
        
        tfile.close();
        cout << "Done h = " << h << "\n";
 //        h=h+2.;

 //   }while(h<=3.0);
    
    
    free (si);
  
}
