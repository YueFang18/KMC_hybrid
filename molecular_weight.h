#include<iomanip>
#include<list>
#include<stdlib.h>
#include<vector>
#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include "my_rand.h"
using namespace std;

double M_total, M_total_square;
double num_M_chain[2] = {0.0};  //number of monomers in each chain

void molecular_weight(double *mol_wt,
                      double &Mn,
                      double &Mw,
                      double *num_M,
                      vector<int> poly_M1,
                      vector<int> poly_M2){
  M_total = 0;
  M_total_square = 0;
  for (int i = 0; i < poly_M1.size(); i++) {
    num_M_chain[0] = poly_M1[i];
    num_M_chain[1] = poly_M2[i];
    num_M[0] += num_M_chain[0];
    num_M[1] += num_M_chain[1];
    M_total += num_M_chain[0] * mol_wt[0] + num_M_chain[1] * mol_wt[1];
    M_total_square += pow((num_M_chain[0] * mol_wt[0] + num_M_chain[1] * mol_wt[1]), 2);
  }
  Mn = M_total / poly_M1.size();
  Mw = M_total_square / M_total;
}

