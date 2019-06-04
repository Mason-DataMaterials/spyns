#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <omp.h>
#include <algorithm>
#include<bits/stdc++.h> 


using namespace std;


random_device rd;
mt19937 engine(rd());
uniform_real_distribution<double> dist(0.0, 1.0);


void positions( int key, int N, double ** r ); 
void initialize_angles(int N, double ** si);
void unique_z( int N, double * z, double * z_u, int* count); 
void sort_positions(int N, double ** r, int& n_layers);
double theta( double * u, double * v );

void xy_per_layer(int N, int n_layers, double ** si);

int main(){
    
    
    int N,i,j,k, n_layers;
    
    N = 1024;
    
    
    double **si = new double* [N];
    for (int i = 0; i < N; i++){
        si[i] = new double [8];
    }
    
    
    //position coordinates and layer indexing to be done in python from input structure file
    positions( 1, N, si );
    initialize_angles(N, si);
    sort_positions(N, si, n_layers);
    
    /*
     *     cout << N << "\n";
     *     cout << n_layers << " z layers in supercell \n"; 
     * 
     *     for (i = 0; i < N; i++ ){
     * 
     *         cout << "Mn" << "\t " << si[i][0] << "\t" << si[i][1] << "\t" << si[i][2] <<"\t"<< si[i][3] << "\t" << si[i][4] << "\t";
     *	       cout                  << si[i][5] <<"\t" << si[i][6] << "\t" << si[i][7] << "\t" << si[i][8] <<"\n" ; 
            }
     */
    
    
    //calculate the relative angles between layers
    xy_per_layer(N, n_layers,  si);
    
    return 0;
    
}


void xy_per_layer(int N, int n_layers, double ** si){
    
    
    int i,j,k;
    double ** theta_rel = new double * [n_layers];
    
    double ** xy_vector = new double * [n_layers];
    for (i = 0; i < n_layers; i++){
        
        xy_vector[i] = new double [2];
        theta_rel[i] = new double [2];
    }
    
    
    for (int i = 0; i < n_layers; i++){
        
        int l_count = 0;
        for (int j = 0; j<N; j++){
            if (si[j][8] == i ){
                xy_vector[i][0] += si[j][2]; 
                xy_vector[i][1] += si[j][3];
                l_count +=1; 
            }
            
        }
        
        xy_vector[i][0] /= l_count;
        xy_vector[i][1] /= l_count;
        
    }
    
    
    for (i = 0; i < n_layers; i++){
        
        
        theta_rel[i][0] = 0.0 ;
        theta_rel[i][1] = 0.0 ;
        
        //angles relative to the bottom layer
        theta_rel[i][0] = theta(xy_vector[0], xy_vector[i]) ;
        
        //angles relative to the last layer
        if (i > 0){
            theta_rel[i][1] = theta(xy_vector[i-1], xy_vector[i]) ;
        }
        
        cout << i << "\t" << xy_vector[i][0] << "\t" << xy_vector[i][1] << "\t" << theta_rel[i][0] << "\t" << theta_rel[i][1]<<"\n" ;
    }
    
    
    return;
    
}


double theta( double * u, double * v )
{
    
    double dot_uv = u[0]*v[0] + u[1]*v[1];
    
    double u_ = sqrt(u[0]*u[0] + u[1]*u[1]);
    double v_ = sqrt(v[0]*v[0] + v[1]*v[1]);
    
    double cos_theta = dot_uv / ( u_ * v_ );
    
    cos_theta = ( cos_theta < -1.0 ? -1.0 : ( cos_theta > 1.0 ? 1.0 : cos_theta ) );
    
    double theta = acos ( cos_theta );
    
    
    return theta * ( 180.00 / M_PI); 
    
}


void unique_z( int N, double * z, double * z_u, int* count) 
{ 
    unordered_set<double> s; 
    
    int j = 0;
    
    for (int i=0; i<N; i++) 
    { 
        if (s.find(z[i])==s.end()) 
        { 
            s.insert(z[i]); 
            z_u[j] = z[i];
            j++; 
        } 
    } 
    
    *count = j;
    
    return;
} 

void sort_positions(int N, double ** r, int& n_layers){
    
    int i,j,k;
    double z_val[N], z_list[N];
    int count;
    
    for (i = 0; i < N; i++){
        z_val[i] = r[i][7];
    }
    
    unique_z( N, z_val, z_list, &count); 
    
    
    for (j =0; j < count; j++ ){
        
        for (i = 0; i < N; i++){
            if(r[i][7] == z_list[j]){
                r[i][8] = j;
            }
            
        }
    }
    
    n_layers = count;
    
    return;
}

void positions( int key, int N, double ** r ) {
    
    double lc = 1.0;
    
    int nbcc = 2;
    int nfcc = 4;
    
    int nc; 
    double * dx; double *dy; double *dz;
    if (key == 1) {
        nc = nbcc;
        
        dx = new double [nc];
        dy = new double [nc];
        dz = new double [nc];
        
        // positions in fcc unit cell
        dx[0] = 0.0; dx[1] = 0.5;
        dy[0] = 0.0; dy[1] = 0.5;
        dz[0] = 0.0; dz[1] = 0.5;
        
    }
    
    else if (key == 2){
        
        nc = nfcc;   
        dx = new double [nc];
        dy = new double [nc];
        dz = new double [nc];
        
        // positions in fcc unit cell
        dx[0] = 0.0; dx[1] = 0.5;dx[2] = 0.5; dx[3] = 0.0;
        dy[0] = 0.0; dy[1] = 0.5;dy[2] = 0.0; dy[3] = 0.5;
        dz[0] = 0.0; dz[1] = 0.0;dz[2] = 0.5; dz[3] = 0.5;
        
    }
    
    // find M large enough to fit N atoms on a lattice
    int M = 0.0;
    while (nc*M*M*M < N)
        M = M+1;
    
    double L = M*lc;
    
    int n = 0;
    for (int x = 0; x < M; x++){
        for (int y = 0; y < M; y++){
            for (int z = 0; z < M; z++)
                for (int p = 0; p < nc; p++){
                    if (n < N){
                        r[n][5] = (x + dx[p]) * lc;
                        r[n][6] = (y + dy[p]) * lc;
                        r[n][7] = (z + dz[p]) * lc;
                        ++n;
                    }
                }
        }
    }
    
    
    double rCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++){
        rCM[0] += r[n][5];
        rCM[1] += r[n][6];
        rCM[2] += r[n][7];
    }
    
    
    for (int i = 0; i < 3; i++)
        rCM[i] /= N;
    
    for (int n = 0; n < N; n++){
        r[n][5] -= rCM[0];
        r[n][6] -= rCM[1];
        r[n][7] -= rCM[2];
    }
    
    
    
    return;
}

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
void initialize_angles(int N, double ** si)
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


