//
//  main.cpp
//  Heisenberg-Model
//
//  Created by Swabir Silayi on 9/3/17.
//  Copyright © 2017 Swabir Silayi. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;


int main(int argc, const char * argv[]) {
    
    int n = 10;
    
    double theta = 0.0;
    double phi = 0.0;
    
    double **S;
    
    S = new double *[n];
    for (int i = 0; i < n; i++)
    {
        S[i] = new double [3];
    }
    
    
    //initialize random orientation
    for (int i = 0; i< n; i++){
        
        theta = ((double)rand() / (RAND_MAX + 1.0))* M_PI;
        phi = ((double)rand() / (RAND_MAX + 1.0))* 2.* M_PI;

        S[i][0] = sin(theta)* cos(phi);
        S[i][1] = sin(theta)*sin(phi);
        S[i][2] = cos(theta);
    }
    
    
    
    double *H;
    H = new double [3];
    H[0] = H[1] = H[2] = 1.0;
    
    double T = 1.0;
    
    for (int step = 0; step < 10000; step++){
        
        double d_theta = ((double)rand() / (RAND_MAX + 1.0))* M_PI;
        double d_phi = ((double)rand() / (RAND_MAX + 1.0))* 2.* M_PI;
        
        
        
        
        
    }
    
    
    return 0;
}
