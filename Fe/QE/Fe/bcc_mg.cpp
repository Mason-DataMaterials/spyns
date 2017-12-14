//
//  bcc_mg.cpp
//  cs-mc
//
//  Created by Swabir Silayi on 11/25/17.
//  Copyright © 2017 Swabir Silayi. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <time.h>


using namespace std;

double u = 1.0;
double d = -1.0;


const int N = 9826;
double ** r;
double ** r_old;
double * dist = new double [N];

double a = 5.3970578; //a.u.
double L;
double rCutOff;

double *Ui;



void initPositions(int key) {
    
    r = new double * [N];
    r_old = new double * [N];
    
    Ui = new double [N];
    
    for (int i = 0; i < N; i++){
        r[i] = new double [3];
        r_old[i] = new double [3];
        Ui[i] = 0.0;
        
    }
    
    
    
    // find M large enough to fit N atoms on an fcc lattice
    int M = 1;
    while (2 * M * M * M < N)
        ++M;
    L = a*M;           // lattice constant of conventional cell
    rCutOff = 0.5*L;
    
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
                        r[n][0] = (x + dx[p]) * a;
                        r[n][1] = (y + dy[p]) * a;
                        r[n][2] = (z + dz[p]) * a;
                        
                        if (key == 0)//nm
                        {
                            Ui[n] = 0.0;
                        }

                        
                        if (key == 1)//fm
                        {
                                Ui[n] = u;
                        }

                        
                        if (key == 2)//afm_1
                        {
                            if (n%2 == 0 )
                                Ui[n] = u;
                            else if (n%2 == 1 )
                                Ui[n] = d;
                        }

                        if (key == 3)//afm_2
                        {
                            if(n==0) Ui[n] = u;
                            else if(n==1) Ui[n] = d;
                            
                            else if (n%2 == 0 && n != 0)
                                Ui[n] = -Ui[n-2];
                            else if (n%2 == 1 && n > 1)
                                Ui[n] = -Ui[n-2];
                        }
                        ++n;
                    }
                }
        }
    }
    
    double rCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            rCM[i] += r[n][i];
    for (int i = 0; i < 3; i++)
        rCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            r[n][i] -= rCM[i];
    
    
    
    double center[3] = {0.0, 0.0, 0.0};
    
    for (int i = 0; i < N; i++){
        
        
        double rSqd;
        double dr[3] = {0.0};
        
        
        dr[0] = -r[i][0] + center[0];
        dr[1] = -r[i][1] + center[1];
        dr[2] = -r[i][2] + center[2];
        
        
        rSqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        dist[i] = sqrt(rSqd);
        
    }
    
    
    int min = dist[0];
    int index = 0;
    for (int i = 1; i < N; i++)
        if (dist[i]<min) {min = dist[i]; index = i;}
    
    double shift_x = r[index][0];
    double shift_y = r[index][1];
    double shift_z = r[index][2];
    
    for (int i = 0; i < N; i++){
        
        r[i][0] =  (r[i][0]-shift_x);
        r[i][1] =  (r[i][1]-shift_y);
        r[i][2] =  (r[i][2]-shift_z);
        
        double dx = center[0] - r[i][0];
        double dy = center[1] - r[i][1];
        double dz = center[2] - r[i][2];
        double rSqd = dx*dx + dy*dy + dz*dz;
        
        dist[i] = sqrt(rSqd);
        
    }
    
    
    
    return;
    
}


double energy(int key){
    
    
    double H = 0.0;
    
    double nn[4] = {0.0, sqrt(3.0)/2.0, 1.0, sqrt(2.0)};
    for (int i = 0; i < 4; i++) nn[i] = nn[i] * a;
    
    double h0 = 0.0;
    double h1 = 0.0;
    double h2  = 0.0;
    
    
    ofstream file;
    if (key == 0) file.open("nm.xyz");
    else if (key == 1) file.open("fm.xyz");
    else if (key == 2) file.open("afm_1.xyz");
    else if (key == 3) file.open("afm_2.xyz");
    
    file <<"15 \n";
    file << "Lattice=0 0 0 0.5 0.5 0.5 Properties=species:S:1:pos:R:3 Time=0 \n";
    
    
    
    
    
    for (int i = 0; i < N; i++){
        
        
        if ( round (dist[i] - nn[0]) == 0.0){
            if(Ui[i] == 0)
                file << "Fe" <<"\t";
            if(Ui[i] == u)
                file << "Fe1" <<"\t";
            if(Ui[i] == d)
                file << "Fe2" <<"\t";
            
            for (int j = 0; j < 3; j++)
                file << r[i][j] <<setw(12)<< "\t";
            
            file << "\n";
            
            h0 += Ui[i];
            
        }
    }
    
    for (int i = 0; i < N; i++){
        
        if (round (dist[i] - nn[1]) == 0.0){
            if(Ui[i] == 0)
                file << "Fe" <<"\t";
            if(Ui[i] == u)
                file << "Fe1" <<"\t";
            if(Ui[i] == d)
                file << "Fe2" <<"\t";
            
            for (int j = 0; j < 3; j++)
                file << r[i][j] <<setw(12)<< "\t";
            
            file << "\n";
            
            h1 += Ui[i];
            
        }
        
    }
    
    
    for (int i = 0; i < N; i++){
        
        
        if ( round(dist[i]-nn[2]) == 0.0){
            
            if(Ui[i] == 0)
                file << "Fe" <<"\t";
            if(Ui[i] == u)
                file << "Fe1" <<"\t";
            if(Ui[i] == d)
                file << "Fe2" <<"\t";
            
            for (int j = 0; j < 3; j++)
                file << r[i][j] <<setw(12)<< "\t";
            
            file << "\n";
            
            
            
            h2 += Ui[i];
        }
        
    }
    
    
    if (key == 0) {    cout << "MG    \t Si \t Sj1 \t Sj2 \n";
                       cout << "NM    \t";
    }
    else if (key == 1) cout << "FM    \t";
    else if (key == 2) cout << "AFM_1 \t";
    else if (key == 3) cout << "AFM_2 \t";
    
    cout << h0/2.0 <<"\t" << h0*h1/2.0 << "\t" << h0*h2/2.0 << "\n";
    
    
    return H;
    
}

void record(string filename){
    ofstream file; file.open(filename.c_str());
    
    
    file << N<<"\n";
    file << "Lattice=0 0 0 0.5 0.5 0.5 Properties=species:S:1:pos:R:3 Time=0 \n";
    for (int i = 0; i < N; i++)
    {
        if(Ui[i] == 0)
            file << "Fe" <<"\t";
        if(Ui[i] == 0.8)
            file << "Fe1" <<"\t";
        if(Ui[i] == -0.1)
            file << "Fe2" <<"\t";
        for (int j = 0; j < 3; j++)
            file << r[i][j] << setw(12)<< "\t";
        
        file << "\n";
        
    }
    
    file.close();
    
    return;
    
}

int main(){
    
    int nm = 0;
    int fm = 1;
    int afm1 = 2;
    int afm2  = 3;
    
    
    
    initPositions(nm);
    record("init_nm.xyz");
    double H = energy(nm);
    
    initPositions(fm);
    record("init_fm.xyz");
    H = energy(fm);
    
    initPositions(afm1);
    record("init_afm1.xyz");
    H = energy(afm1);
    
    initPositions(afm2);
    record("init_afm2.xyz");
    H = energy(afm2);

    
    return 0;
    
}
