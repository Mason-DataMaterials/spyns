#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int N = 10; //lattice size
int *** si; //spin configuration

double H; //Energy
double M; //Magnetization

double h; //external field
double T; //temperature values

int MCSteps = 200000;

void initialize(){
    
    si =new int ** [N];
    
    for (int x = 0; x < N; x++){
        si[x] = new int * [N];
        
        for (int y = 0; y < N; y++){
            si[x][y] = new int [N];
            
            for (int z = 0; z < N; z++){
                
                double r = (double) rand()/RAND_MAX;
                si[x][y][z] = 1;
            }
        }
    }
    
    return;
}

int neighbors(int * point ){
    
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

double Magnetization(){
    
    double M = 0.0;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                M = M+ (double)si[i][j][k];
    
    return M;
    
}


double Energy(){
    
    double H = 0.0;
    int * point = new int [3];
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                point[0] = i; point[1] = j; point[2]=k;
                
                int sj = neighbors(point);
                H = H + ( - si[i][j][k]*(h+sj));
            }
    
    free (point);
    
    return H / 2.;
}

double deltaE (int * point){
    
    int sj = neighbors(point);
    double dE = -2. * si[point[0]][point[1]][point[2]]*(h+sj);
    
    return dE;
    
}

void MC_move(){
    
    double beta = 1./T;
    int * point = new int [3];
    int i = rand()%N; int j = rand()%N; int k = rand()%N;
    
    point[0] = i; point[1] = j; point[2]=k;
    
    
    double dE = deltaE(point);
    double dm=0.0,dh=0.0;
    double p = T*log((double)rand()/RAND_MAX);
    if (dE > p)
    {
        
        si[i][j][k] *= -1;
        dm=2.*si[i][j][k];
        dh=dE;
    }
    
    M = M+dm;
    H = H+dh;
    
    free (point);
    return;
}

void collect_data(double * data){
    data[0] = data[0]+M;
    data[1] = data[1]+H;
    data[2] = data[2]+M*M;
    data[3] = data[3]+H*H;
    data[4] = data[4]+fabs(M);
}

void configuration(){
    
    ofstream cfile;
    std::ostringstream cf;
    cf <<"h="<<h<<"/configurations/t="<<int(T)<<".txt" ;
    std::string cnf = cf.str();
    
    cfile.open(cnf.c_str());
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                cfile<<i <<"\t" << j << "\t" << k << "\t" << si[i][j][k] << "\n";
            }
    
    return;
    
}


int main(){
    h = 1.;
    
    do{
        ofstream tfile;
        std::ostringstream tf;
        tf <<"h="<<h<<"/mg.txt" ;
        std::string tmf = tf.str();
        
        tfile.open(tmf.c_str());

        
        for (T=1.0; T <= 12.; T++){
            
            initialize();
            M = Magnetization();
            H = Energy();
            
            double * data;
            data = new double [5];
            for (int i=0; i < 5; i++) data[i] = 0.0;
            
            for (int step = 0; step < MCSteps; step++){
                
                MC_move();
                collect_data(data);
            }
            
            double norm   = 1.0/(MCSteps*N);
            double M_avg  = data[0]*norm;
            double H_avg  = data[1]*norm;
            double M_2avg = data[2]*norm;
            double H_2avg = data[3]*norm;
            double M_abs  = data[4]*norm;
            
            tfile << T << "\t" << M_avg << "\t" << H_avg << "\n";
            configuration();
        }
        
        
        
        h=h+1.;
        
    }while(h<=10.0);
    
    
    return 0;
}
