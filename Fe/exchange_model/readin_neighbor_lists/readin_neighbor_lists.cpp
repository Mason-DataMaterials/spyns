#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main () {
    
    

      
    
  string line;
  ifstream infile ("neighbor_lists.txt");
  
  int N;
  int num_1st_nn,num_2nd_nn;
  double * num_nn;
  
  num_nn = new double * [2];
  
  if (infile.is_open())
  {
   
      infile>>N;
      infile>>num_1st_nn>>num_2nd_nn;
      
      num_nn[0] = num_1st_nn;
      num_nn[1] = num_2nd_nn;
      
      //populate neighbor lists
      double **first_nn = new double * [N];
      double **second_nn = new double * [N];
      for (int i = 0; i < N; i++)
      {
          first_nn[i] = new double [num_1st_nn];
          second_nn[i] = new double [num_2nd_nn];
      }  
      
      int nn_per_site = num_1st_nn + num_2nd_nn;
      
      for (int i = 0; i < N; i++){
          cout << i << "\t"; 
          
          for (int j = 0; j < num_1st_nn; j++){
              infile >> first_nn[i][j]; 
              
              cout << first_nn[i][j] << "\t";
          }
          
          for (int j = 0; j < num_2nd_nn; j++){
              infile >> second_nn[i][j]; 
              cout << second_nn[i][j] << "\t";
          }
          
          
          cout << "\n";
      }
        
       
      
      infile.close();
  }

  else cout << "Unable to open file"; 

  return 0;
}
