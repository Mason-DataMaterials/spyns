//
//  Heisenberg-Model
//  heisenberg.cpp
//

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

int N = 20; //lattice size
double **** si; //spin configuration

double H; //Energy
double M; //Magnetization
double * m;

double * h; //external field
double T; //temperature values

double theta;
double phi;

int MCSteps = 2525000;

void initialize(){
    
    si =new double *** [N];
    
    for (int x = 0; x < N; x++){
        si[x] = new double ** [N];
        
        for (int y = 0; y < N; y++){
            si[x][y] = new double * [N];
            
            for (int z = 0; z < N; z++){
                
                si[x][y][z] = new double [3];
                //(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
                theta = 1.0/cos(1.0-2.0 * (double) rand()/RAND_MAX);
                phi = M_PI*2.0* (double)rand()/RAND_MAX;
                //Sx
                si[x][y][z][0] = sin(theta)*cos(phi);
                //Sy
                si[x][y][z][1] = sin(theta)*sin(phi);
                //Sz
                si[x][y][z][2] = cos(phi);
                
            }
        }
    }
    
    M = 0.0;
    m = new double [N*N*N];
    for (int i =0; i <N*N*N; i++) m[i] = 0.0;
    
    
    return;
}

void configuration(){
    
    ofstream cfile;
    std::ostringstream cf;
    cf <<"h="<<int(h[0])<<"/configurations/t="<<int(T)<<".txt" ;
    //cf << int(T) << ".txt";
    std::string cnf = cf.str();
    
    cfile.open(cnf.c_str());
    
    int n = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                cfile<< i <<"\t" << j << "\t" << k << "\t" << m[n] <<"\t"
                << si[i][j][k][0] << "\t"<< si[i][j][k][1] << "\t"<< si[i][j][k][2] << "\n";
                
                n++;
            }
    
    return;
    
}

double Magnetization(){
    
    int n = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                for(int p = 0; p <3; p++)
                {
                    m[n] += si[i][j][k][p]*si[i][j][k][p];
                }
                m[n] = sqrt(m[n]);
                M = M+m[n];
                n++;
            }
    
    
    return M;
    
}

int neighbors(int * point, int i){
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    int xp = (x+N-1) % N;
    int xm = (x+1) % N;
    int yp = (y+N-1) % N;
    int ym = (y+1) % N;
    int zp = (z+N-1) % N;
    int zm = (z+1) % N;
    
    
    int sj = si[xp][y][z][i] + si[xm][y][z][i] + si[x][yp][z][i] + si[x][ym][z][i] + si[x][y][zp][i] + si[x][y][zm][i];
    
    return sj;
}

double Energy(){
    
    double * Hi = new double [3];
    for (int i =0; i <3; i++) Hi[i] = 0.0;
    
    int * point = new int [3];
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){
                point[0] = i; point[1] = j; point[2]=k;
                
                int sj_x = neighbors(point, 0);
                int sj_y = neighbors(point, 1);
                int sj_z = neighbors(point, 2);
                
                
                Hi[0] = Hi[0] + ( - si[i][j][k][0]*(h[0]+sj_x));
                Hi[1] = Hi[1] + ( - si[i][j][k][1]*(h[1]+sj_y));
                Hi[2] = Hi[2] + ( - si[i][j][k][2]*(h[2]+sj_z));
            }
    
    double Ht = (Hi[0] + Hi[1] + Hi[2]);
    free (point);
    return Ht / 2.;
}

double deltaE (int * point){
    
    int sj_x = neighbors(point, 0);
    int sj_y = neighbors(point, 1);
    int sj_z = neighbors(point, 2);
    
    double dE_x = -2. * si[point[0]][point[1]][point[2]][0]*(h[0]+sj_x);
    double dE_y = -2. * si[point[0]][point[1]][point[2]][1]*(h[1]+sj_y);
    double dE_z = -2. * si[point[0]][point[1]][point[2]][2]*(h[2]+sj_z);
    
    double dE = dE_x + dE_y + dE_z;
    
    return dE;
    
}

double delta_M(int * point){
    double dm[3] = {0.0};
    
    int x = point[0];
    int y = point[1];
    int z = point[2];
    
    dm[0] = si[x][y][z][0];
    dm[1] = si[x][y][z][1];
    dm[2] = si[x][y][z][2];
    
    double dM = 2. * sqrt(dm[0]*dm[0] + dm[1]*dm[1] + dm[2]*dm[2] );
    
    
    return dM;
}


void MC_move(){
    
    int * point = new int [3];
    int i = rand()%N; int j = rand()%N; int k = rand()%N;
    
    point[0] = i; point[1] = j; point[2]=k;
    
    //theta update, phi constant
    double dE = deltaE(point);
    
    double p = T*log((double)rand()/RAND_MAX);
    if (dE > p)
    {
        theta = 1.0/cos(1.0-2.0 * (double) rand()/RAND_MAX);
        
        //Sx
        si[i][j][k][0] = sin(theta)*cos(phi);
        //Sy
        si[i][j][k][1] = sin(theta)*sin(phi);
        //Sz
        si[i][j][k][2] = cos(phi);
        
        M = M+delta_M(point);
        H = H+dE;
    }
    
    //phi update, theta constant
    dE = deltaE(point);
    
    p = T*log((double)rand()/RAND_MAX);
    if (dE > p)
    {
        phi = M_PI*2.0* (double)rand()/RAND_MAX;
        //Sx
        si[i][j][k][0] = sin(theta)*cos(phi);
        //Sy
        si[i][j][k][1] = sin(theta)*sin(phi);
        //Sz
        si[i][j][k][2] = cos(phi);
        
        M = M+delta_M(point);
        H = H+dE;
    }
    
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

int main(int argc, const char * argv[]) {
    
    double ht = 0.;
    
    for (; ht <=12; ht++)
    {
        h = new double [3];
        h[0] = ht; h[1] = ht; h[2] = ht;
        
        ofstream tfile;
        std::ostringstream tf;
        tf <<"h="<<int(ht)<<"/mg.txt" ;
        std::string tmf = tf.str();
        
        tfile.open(tmf.c_str());

        

        for (T =0.; T <=12; T++){
            
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
            double norm   = 1.0/(MCSteps*N*N*N);
            double M_avg  = data[0]*norm;
            double H_avg  = data[1]*norm;
            double M_2avg = data[2]*norm;
            double H_2avg = data[3]*norm;
            double M_abs  = data[4]*norm;
            
            tfile << T << "\t" << M_avg << "\t" << H_avg << "\n";
            configuration();
            
            free(data);
        }
        
        tfile.close();
        
    }
    return 0;
}
