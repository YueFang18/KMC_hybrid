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
struct chain {
  vector<int> M1, M2;
  int dyad[4];  //AA,AB,BA,BB
  int triad[8]; //AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
};

struct poly_chain {
  poly_chain(){};
  vector<int> chain1, chain2;
  vector<int> num_M1, num_M2;
//	vector<int> AA,AB,BA,BB;
} poly;

std::vector<int> M1M1;
std::vector<int> M1M2;
std::vector<int> M2M1;
std::vector<int> M2M2;

typedef vector<chain> chain_pool;
chain_pool allchains;          //completed segments and dead chains
chain new_chain;
int r1, r2, rmax;                       //random numbers for chain selection

vector<int>::iterator p_chain;
int c_tab1, c_tab2;               //2chain labels

void efficient_explicit_sequence_record_number(int reaction_index,
                                               double *num_monomer,
                                               double *species,
                                               int &num_chain,
                                               vector<int> &M1M1,
                                               vector<int> &M1M2,
                                               vector<int> &M2M1,
                                               vector<int> &M2M2,
                                               chain_pool &allchains,
                                               poly_chain &poly,
                                               double *tot_dyad,
                                               double *tot_triad) {
  int M1_counter = 0;
  int M2_counter = 0;
/* if((reaction_index==3) || (reaction_index==4) || (reaction_index==5) || (reaction_index==6) )
	{
          num_chain=num_chain+1;
          allchains.push_back(new_chain);
    }
*/
  switch (reaction_index) {
    /************** Initiator Decomposition **************/
    case 50:              //I2=2I*.
      species[0] -= 1;
      species[1] += 2;
      break;
      /**************** Chain Initiation*************/
    case 51 :             //I*+M1=I-M1*
      species[1] -= 1;
      species[2] -= 1;
      species[4] += 1;
      num_monomer[0] += 1;
      break;
    case 52 :       //I*+M2=I-M2*
      species[1] -= 1;
      species[3] -= 1;
      species[5] += 1;
      num_monomer[1] += 1;
      break;
      /**************** Chain Propagation*************/
    case 53 :             //IM1*+M1=IM1M1*
      species[2] -= 1;
      species[4] -= 1;
      species[6] += 1;
      num_monomer[0] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(2);            //add 2 M1's
      allchains[allchains.size() - 1].dyad[0] = 1;                //add 1 AA
      allchains[allchains.size() - 1].dyad[1] = 0;
      allchains[allchains.size() - 1].dyad[2] = 0;
      allchains[allchains.size() - 1].dyad[3] = 0;
      for (int i = 0; i < 8; i++) {
        allchains[allchains.size() - 1].triad[i] = 0;
      }
      M1M1.push_back(num_chain);                                //Type 1 Chain
      break;
    case 54 :                  //IM1*+M2=IM1M2*
      species[3] -= 1;
      species[4] -= 1;
      species[7] += 1;
      num_monomer[1] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);

      allchains[allchains.size() - 1].M1.push_back(1);            //add 1 M1
      allchains[allchains.size() - 1].M2.push_back(1);            //add 1 M2
      allchains[allchains.size() - 1].dyad[0] = 0;                //add 1 AA
      allchains[allchains.size() - 1].dyad[1] = 1;
      allchains[allchains.size() - 1].dyad[2] = 0;
      allchains[allchains.size() - 1].dyad[3] = 0;
      for (int i = 0; i < 8; i++) {
        allchains[allchains.size() - 1].triad[i] = 0;
      }
      M1M2.push_back(num_chain);                                //Type2 Chain
      break;
    case 55 :              //IM2*+M1---->IM2M1*
      species[2] -= 1;
      species[5] -= 1;
      species[8] += 1;
      num_monomer[0] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);

      allchains[allchains.size() - 1].M1.push_back(1);            //add 1 M1
      allchains[allchains.size() - 1].M2.push_back(1);            //add 1 M2
      allchains[allchains.size() - 1].dyad[0] = 0;                //add 1 AA
      allchains[allchains.size() - 1].dyad[1] = 0;
      allchains[allchains.size() - 1].dyad[2] = 1;
      allchains[allchains.size() - 1].dyad[3] = 0;
      for (int i = 0; i < 8; i++) {
        allchains[allchains.size() - 1].triad[i] = 0;
      }
      M2M1.push_back(num_chain);                                //Type3 Chain
      break;

    case 56 :            //IM2*+M2---->IM2M2*
      species[3] -= 1;
      species[5] -= 1;
      species[9] += 1;
      num_monomer[1] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);

      allchains[allchains.size() - 1].M2.push_back(2);            //add 2 M2
      allchains[allchains.size() - 1].dyad[0] = 0;                //add 1 AA
      allchains[allchains.size() - 1].dyad[1] = 0;
      allchains[allchains.size() - 1].dyad[2] = 0;
      allchains[allchains.size() - 1].dyad[3] = 1;
      for (int i = 0; i < 8; i++) {
        allchains[allchains.size() - 1].triad[i] = 0;
      }
      M2M2.push_back(num_chain);                                //Type4 Chain
      break;
      /****************Propagation with Penultimate effect*************/
    case 0 :                        //IM1M1*+M1=IM1M1*
      species[2] -= 1;
      species[6] -= 1;
      species[6] += 1;
      num_monomer[0] += 1;
      r1 = my_rand(M1M1.size());
      c_tab1 = M1M1[r1];
      allchains[c_tab1].M1.back()++;
      allchains[c_tab1].dyad[0]++;                    //+AA
      allchains[c_tab1].triad[0]++;                    // (AAA),AAB,ABA,ABB,BAA,BAB,BBA,BBB
      break;
    case 1 :                          //IM1M2*+M1=IM2M1*
      species[2] -= 1;
      species[7] -= 1;
      species[8] += 1;
      num_monomer[0] += 1;
      r1 = my_rand(M1M2.size());
      c_tab1 = M1M2[r1];
      allchains[c_tab1].M1.push_back(1);
      allchains[c_tab1].dyad[2]++;                    //+BA
      allchains[c_tab1].triad[2]++;                    // AAA,AAB,(ABA),ABB,BAA,BAB,BBA,BBB
      M2M1.push_back(c_tab1);
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);

      break;
    case 2 :                          //IM2M1*+M1=IM1M1*
      species[2] -= 1;
      species[8] -= 1;
      species[6] += 1;
      num_monomer[0] += 1;
      r1 = my_rand(M2M1.size());
      c_tab1 = M2M1[r1];
      allchains[c_tab1].M1.back()++;
      allchains[c_tab1].dyad[0]++;                    //+AA
      allchains[c_tab1].triad[4]++;                    // AAA,AAB,ABA,ABB,(BAA),BAB,BBA,BBB
      M1M1.push_back(c_tab1);
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);

      break;
    case 3 :                         //IM2M2*+M1=IM2M1*
      species[2] -= 1;
      species[9] -= 1;
      species[8] += 1;
      num_monomer[0] += 1;
      r1 = my_rand(M2M2.size());
      c_tab1 = M2M2[r1];
      allchains[c_tab1].M1.push_back(1);
      allchains[c_tab1].dyad[2]++;                    //+BA
      allchains[c_tab1].triad[6]++;                    // AAA,AAB,ABA,ABB,BAA,BAB,(BBA),BBB
      M2M1.push_back(c_tab1);
      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      break;
    case 4 :                        //IM1M1*+M2=IM1M2*
      species[3] -= 1;
      species[6] -= 1;
      species[7] += 1;
      num_monomer[1] += 1;
      r1 = my_rand(M1M1.size());
      c_tab1 = M1M1[r1];
      allchains[c_tab1].M2.push_back(1);
      allchains[c_tab1].dyad[1]++;                    //+AB
      allchains[c_tab1].triad[1]++;                    // AAA,(AAB),ABA,ABB,BAA,BAB,BBA,BBB
      M1M2.push_back(c_tab1);
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      break;
    case 5 :                        //IM1M2*+M2=IM2M2*
      species[3] -= 1;
      species[7] -= 1;
      species[9] += 1;
      num_monomer[1] += 1;
      r1 = my_rand(M1M2.size());
      c_tab1 = M1M2[r1];
      allchains[c_tab1].M2.back()++;
      allchains[c_tab1].dyad[3]++;                    //+BB
      allchains[c_tab1].triad[3]++;                    // AAA,AAB,ABA,(ABB),BAA,BAB,BBA,BBB
      M2M2.push_back(c_tab1);
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      break;
    case 6 :                         //IM2M1*+M2=IM1M2*
      species[3] -= 1;
      species[8] -= 1;
      species[7] += 1;
      num_monomer[1] += 1;
      r1 = my_rand(M2M1.size());
      c_tab1 = M2M1[r1];
      allchains[c_tab1].M2.push_back(1);
      allchains[c_tab1].dyad[1]++;                    //+AB
      allchains[c_tab1].triad[5]++;                    // AAA,AAB,ABA,ABB,BAA,(BAB),BBA,BBB
      M1M2.push_back(c_tab1);
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      break;
    case 7 :                          //IM2M2*+M2=IM2M2*
      species[3] -= 1;
      species[9] -= 1;
      species[9] += 1;
      num_monomer[1] += 1;
      r1 = my_rand(M2M2.size());
      c_tab1 = M2M2[r1];
      allchains[c_tab1].M2.back()++;
      allchains[c_tab1].dyad[3]++;                    //+BB
      allchains[c_tab1].triad[7]++;                    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,(BBB)
      break;
      /****************ReCombination of type I-MM.+I-MM. --> dead end polymer*************/
    case 8 :                           //IM1M1*+IM1M1*=POLY
      species[6] -= 1;
      species[6] -= 1;
      species[10] += 1;
//  pick two random chains
      do {
        r1 = my_rand(M1M1.size());
        r2 = my_rand(M1M1.size());
      } while (r1 == r2);

//  mark the two chains with poly category
      c_tab1 = M1M1[r1];
      c_tab2 = M1M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0] + 1;        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
//  triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0] + 2;
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M1M1.begin() + r2;
      M1M1.erase(p_chain);

//  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;

    case 9 :                          //IM1M1*+IM1M2*=poly;
      species[6] -= 1;
      species[7] -= 1;
      species[10] += 1;

      r1 = my_rand(M1M1.size());
      r2 = my_rand(M1M2.size());

      c_tab1 = M1M1[r1];
      c_tab2 = M1M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2] + 1;        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4] + 1;
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2] + 1;
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M1M2.begin() + r2;
      M1M2.erase(p_chain);
      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);

      break;
    case 10 :                          //IM1M1*+IM2M1*=POLY
      species[6] -= 1;
      species[8] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M1.size());
      r2 = my_rand(M2M1.size());

      c_tab1 = M1M1[r1];
      c_tab2 = M2M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0] + 1;        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
//triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0] + 1;
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4] + 1;
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M2M1.begin() + r2;
      M2M1.erase(p_chain);
//  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);

      break;
    case 11 :                       //IM1M1*+IM2M2*=POLY;
      species[6] -= 1;
      species[9] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M1.size());
      r2 = my_rand(M2M2.size());

      c_tab1 = M1M1[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2] + 1;        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4] + 1;
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6] + 1;
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      //delete active segments r1 and r2
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);

      break;

    case 12 :                          //IM1M2*+IM1M2*=POLY
      species[7] -= 1;
      species[7] -= 1;
      species[10] += 1;
      do {
        r1 = my_rand(M1M2.size());
        r2 = my_rand(M1M2.size());
      } while (r1 == r2);

      c_tab1 = M1M2[r1];
      c_tab2 = M1M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3] + 1;        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6] + 1;
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3] + 1;
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      p_chain = M1M2.begin() + r2;
      M1M2.erase(p_chain);
      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);

      break;
    case 13 :                          //IM1M2*+IM2M1*=POLY
      species[7] -= 1;
      species[8] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M2.size());
      r2 = my_rand(M2M1.size());

      c_tab1 = M1M2[r1];
      c_tab2 = M2M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1] + 1;        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2] + 1;
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5] + 1;
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      p_chain = M2M1.begin() + r2;
      M2M1.erase(p_chain);

      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 14 :                          //IM1M2*+IM2M2*=POLY
      species[7] -= 1;
      species[9] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M2.size());
      r2 = my_rand(M2M2.size());

      c_tab1 = M1M2[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3] + 1;        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1] + 1;
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7] + 1;

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);

      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 15 :                          //IM2M1*+IM2M1*=POLY
      species[8] -= 1;
      species[8] -= 1;
      species[10] += 1;
      do {
        r1 = my_rand(M2M1.size());
        r2 = my_rand(M2M1.size());
      } while (r1 == r2);

      c_tab1 = M2M1[r1];
      c_tab2 = M2M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0] + 1;        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4] + 1;
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1] + 1;
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      p_chain = M2M1.begin() + r2;
      M2M1.erase(p_chain);
      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 16 :                          //IM2M1*+IM2M2*=POLY
      species[8] -= 1;
      species[9] -= 1;
      species[10] += 1;
      r1 = my_rand(M2M1.size());
      r2 = my_rand(M2M2.size());

      c_tab1 = M2M1[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2] + 1;        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6] + 1;
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5] + 1;
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 17 :                          //IM2M2*+IM2M2*=POLY
      species[9] -= 1;
      species[9] -= 1;
      species[10] += 1;
      do {
        r1 = my_rand(M2M2.size());
        r2 = my_rand(M2M2.size());
      } while (r1 == r2);

      c_tab1 = M2M2[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(c_tab2);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[2];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[1];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3] + 1;        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[4];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[6];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[1];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[3];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7] + 2;

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      //  update the numbers of M1 and M2 in poly chains
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
      /****************ReCombination of type I-M.+I-M. --> dead end polymer*************/
    case 18 :                         //IM1*+IM1*=POLY
      species[4] -= 1;
      species[4] -= 1;
      species[10] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      for (int i = 0; i < 3; i++) { allchains[allchains.size() - 1].dyad[i] = 0; }
      for (int i = 0; i < 8; i++) { allchains[allchains.size() - 1].triad[i] = 0; }

      allchains[allchains.size() - 1].M1.push_back(2);            //add 2 M1
      allchains[allchains.size() - 1].dyad[0]++;
      tot_dyad[0]++;
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(2);
      poly.num_M2.push_back(0);
      break;
    case 19 :                               //IM1*+IM2*=POLY
      species[4] -= 1;
      species[5] -= 1;
      species[10] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      for (int i = 0; i < 3; i++) { allchains[allchains.size() - 1].dyad[i] = 0; }
      for (int i = 0; i < 8; i++) { allchains[allchains.size() - 1].triad[i] = 0; }
      allchains[allchains.size() - 1].M1.push_back(1);
      allchains[allchains.size() - 1].M2.push_back(1);
//	 allchains[allchains.size()-1].dyad[1]++;
//	 tot_dyad[1]++;
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(1);
      break;
    case 20 :                               //IM2*+IM2*=POLY
      species[5] -= 1;
      species[5] -= 1;
      species[10] += 1;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      for (int i = 0; i < 3; i++) { allchains[allchains.size() - 1].dyad[i] = 0; }
      for (int i = 0; i < 8; i++) { allchains[allchains.size() - 1].triad[i] = 0; }
      allchains[allchains.size() - 1].M2.push_back(2);            //add 2 M2
      allchains[allchains.size() - 1].dyad[3]++;
      tot_dyad[3]++;
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(2);
      break;
      /****************ReCombination of type I-M.+I-MM. --> dead end polymer*************/
    case 21 :                                //IM1*+IM1M1*=POLY
      species[4] -= 1;
      species[6] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M1.size());
      c_tab1 = M1M1[r1];
      allchains[c_tab1].M1.back()++;
      allchains[c_tab1].dyad[0]++;
      allchains[c_tab1].triad[0]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_dyad[0]++;
      tot_triad[0]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter + 1);
      poly.num_M2.push_back(M2_counter);

      break;
    case 22 :                              //IM1*+IM1M2*=POLY
      species[4] -= 1;
      species[7] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M2.size());
      c_tab1 = M1M2[r1];
      allchains[c_tab1].M1.push_back(1);
      allchains[c_tab1].dyad[2]++;

      tot_dyad[2]++;
      allchains[c_tab1].triad[2]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[2]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);

      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter + 1);
      poly.num_M2.push_back(M2_counter);
      break;
    case 23 :                              //IM1*+IM2M1*=POLY
      species[4] -= 1;
      species[8] -= 1;
      species[10] += 1;
      r1 = my_rand(M2M1.size());
      c_tab1 = M2M1[r1];
      allchains[c_tab1].M1.back()++;
      allchains[c_tab1].dyad[0]++;
      tot_dyad[0]++;
      allchains[c_tab1].triad[4]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[4]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);

      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter + 1);
      poly.num_M2.push_back(M2_counter);

      break;
    case 24 :                              //IM2M2*+IM1*=POLY
      species[4] -= 1;
      species[9] -= 1;
      species[10] += 1;
      r1 = my_rand(M2M2.size());
      c_tab1 = M2M2[r1];
      allchains[c_tab1].M1.push_back(1);
      allchains[c_tab1].dyad[2]++;
      tot_dyad[2]++;
      allchains[c_tab1].triad[6]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[6]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter + 1);
      poly.num_M2.push_back(M2_counter);
      break;
    case 25 :                              //IM1M1*+M2=POLY
      species[5] -= 1;
      species[6] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M1.size());
      c_tab1 = M1M1[r1];
      allchains[c_tab1].M2.push_back(1);
      allchains[c_tab1].dyad[1]++;
      tot_dyad[1]++;
      allchains[c_tab1].triad[1]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter + 1);
      break;
    case 26 :                              //IM1M2*+M2=POLY
      species[5] -= 1;
      species[7] -= 1;
      species[10] += 1;
      r1 = my_rand(M1M2.size());
      c_tab1 = M1M2[r1];
      allchains[c_tab1].M2.back()++;
      allchains[c_tab1].dyad[3]++;
      tot_dyad[3]++;
      allchains[c_tab1].triad[3]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[3]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter + 1);
      break;
    case 27 :                              //IM2M1*+M2=POLY
      species[5] -= 1;
      species[8] -= 1;
      species[10] += 1;
      r1 = my_rand(M2M1.size());
      c_tab1 = M2M1[r1];
      allchains[c_tab1].M2.push_back(1);
      allchains[c_tab1].dyad[1]++;
      tot_dyad[1]++;
      allchains[c_tab1].triad[5]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[5]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter + 1);
      break;
    case 28 :                              //IM2M2*+M2=POLY
      species[5] -= 1;
      species[9] -= 1;
      species[10] += 1;
      r1 = my_rand(M2M2.size());
      c_tab1 = M2M2[r1];
      allchains[c_tab1].M2.back()++;
      allchains[c_tab1].dyad[3]++;
      tot_dyad[3]++;
      allchains[c_tab1].triad[7]++; // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[7]++;
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter + 1);
      break;
      /****************Disproportination of type I-MM.+I-MM. --> dead end polymer (1 to 2 radical transfer)*************/
    case 29 :                                //IM1M1*+IM1M1*=POLY+M
      species[6] -= 1;
      species[6] -= 1;
      species[10] += 2;
      do {
        r1 = my_rand(M1M1.size());
        r2 = my_rand(M1M1.size());
      } while (r1 == r2);

      c_tab1 = M1M1[r1];
      c_tab2 = M1M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);

//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M1M1.begin() + r2;
      M1M1.erase(p_chain);

      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);

      break;
    case 30 :                          //IM1M1*+IM1M2*=poly;
      species[6] -= 1;
      species[7] -= 1;
      species[10] += 2;

      r1 = my_rand(M1M1.size());
      r2 = my_rand(M1M2.size());

      c_tab1 = M1M1[r1];
      c_tab2 = M1M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);
//  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M1M2.begin() + r2;
      M1M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 31 :                          //IM1M1*+IM2M1*=POLY
      species[6] -= 1;
      species[8] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M1.size());
      r2 = my_rand(M2M1.size());

      c_tab1 = M1M1[r1];
      c_tab2 = M2M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);
      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      //delete active segments r1 and r2
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M2M1.begin() + r2;
      M2M1.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 32 :                       //IM1M1*+IM2M2*=POLY;
      species[6] -= 1;
      species[9] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M1.size());
      r2 = my_rand(M2M2.size());

      c_tab1 = M1M1[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 33 :                          //IM1M2*+IM1M2*=POLY
      species[7] -= 1;
      species[7] -= 1;
      species[10] += 2;
      do {
        r1 = my_rand(M1M2.size());
        r2 = my_rand(M1M2.size());
      } while (r1 == r2);

      c_tab1 = M1M2[r1];
      c_tab2 = M1M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      p_chain = M1M2.begin() + r2;
      M1M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 34 :                          //IM1M2*+IM2M1*=POLY
      species[7] -= 1;
      species[8] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M2.size());
      r2 = my_rand(M2M1.size());

      c_tab1 = M1M2[r1];
      c_tab2 = M2M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      p_chain = M2M1.begin() + r2;
      M2M1.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 35 :                          //IM1M2*+IM2M2*=POLY
      species[7] -= 1;
      species[9] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M2.size());
      r2 = my_rand(M2M2.size());

      c_tab1 = M1M2[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 36 :                          //IM2M1*+IM2M1*=POLY
      species[8] -= 1;
      species[8] -= 1;
      species[10] += 2;
      do {
        r1 = my_rand(M2M1.size());
        r2 = my_rand(M2M1.size());
      } while (r1 == r2);

      c_tab1 = M2M1[r1];
      c_tab2 = M2M1[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      p_chain = M2M1.begin() + r2;
      M2M1.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 37 :                          //IM2M1*+IM2M2*=POLY
      species[8] -= 1;
      species[9] -= 1;
      species[10] += 2;
      r1 = my_rand(M2M1.size());
      r2 = my_rand(M2M2.size());

      c_tab1 = M2M1[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
    case 38 :                          //IM2M2*+IM2M2*=POLY
      species[9] -= 1;
      species[9] -= 1;
      species[10] += 2;
      do {
        r1 = my_rand(M2M2.size());
        r2 = my_rand(M2M2.size());
      } while (r1 == r2);

      c_tab1 = M2M2[r1];
      c_tab2 = M2M2[r2];
      poly.chain1.push_back(c_tab1);
      poly.chain1.push_back(c_tab2);
      poly.chain2.push_back(-1);
      poly.chain2.push_back(-1);

      //  add the dyads to total dyads
      tot_dyad[0] += allchains[c_tab1].dyad[0] + allchains[c_tab2].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1] + allchains[c_tab2].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2] + allchains[c_tab2].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3] + allchains[c_tab2].dyad[3];        //BB
      //triad
      tot_triad[0] += allchains[c_tab1].triad[0] + allchains[c_tab2].triad[0];
      tot_triad[1] += allchains[c_tab1].triad[1] + allchains[c_tab2].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2] + allchains[c_tab2].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3] + allchains[c_tab2].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4] + allchains[c_tab2].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5] + allchains[c_tab2].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6] + allchains[c_tab2].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7] + allchains[c_tab2].triad[7];

      // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
//delete active segments r1 and r2
      if (r1 < r2) {
        rmax = r2;
        r2 = r1;
        r1 = rmax;
      }

      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      p_chain = M2M2.begin() + r2;
      M2M2.erase(p_chain);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      M1_counter = 0;
      poly.num_M2.push_back(M2_counter);
      M2_counter = 0;
      for (int w = 0; w < allchains[c_tab2].M1.size(); w++) { M1_counter += allchains[c_tab2].M1[w]; };
      for (int w = 0; w < allchains[c_tab2].M2.size(); w++) { M2_counter += allchains[c_tab2].M2[w]; };
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      break;
      /****************Disproportination of type I-M.+I-M. --> dead end polymer*************/
    case 39 :                         //IM1*+IM1*=POLY
      species[4] -= 1;
      species[4] -= 1;
      species[10] += 2;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);            //add 1 M1
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);
      num_chain = num_chain + 1;                                    //add a second chain
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);            //add 1 M1
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);
      break;
    case 40 :                               //IM1*+IM2*=POLY
      species[4] -= 1;
      species[5] -= 1;
      species[10] += 2;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M2.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);
      break;
    case 41 :                               //IM2*+IM2*=POLY
      species[5] -= 1;
      species[5] -= 1;
      species[10] += 2;
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M2.push_back(1);            //add 2 M2
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M2.push_back(1);            //add 2 M2
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);

      break;
      /****************Disprop of type I-M.+I-MM. --> dead end polymer*************/
    case 42 :                                //IM1*+IM1M1*=POLY
      species[4] -= 1;
      species[6] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M1.size());
      c_tab1 = M1M1[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];

      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);

      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);

      break;
    case 43 :                              //IM1*+IM1M2*=POLY
      species[4] -= 1;
      species[7] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M2.size());
      c_tab1 = M1M2[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      break;
    case 44 :                              //IM1*+IM2M1*=POLY
      species[4] -= 1;
      species[8] -= 1;
      species[10] += 2;
      r1 = my_rand(M2M1.size());
      c_tab1 = M2M1[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      break;
    case 45 :                              //IM2M2*+IM1*=POLY
      species[4] -= 1;
      species[9] -= 1;
      species[10] += 2;
      r1 = my_rand(M2M2.size());
      c_tab1 = M2M2[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(1);
      poly.num_M2.push_back(0);
      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      break;
    case 46 :                              //IM1M1*+IM2*=POLY
      species[5] -= 1;
      species[6] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M1.size());
      c_tab1 = M1M1[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);
      p_chain = M1M1.begin() + r1;
      M1M1.erase(p_chain);
      break;
    case 47 :                              //IM1M2*+IM2*=POLY
      species[5] -= 1;
      species[7] -= 1;
      species[10] += 2;
      r1 = my_rand(M1M2.size());
      c_tab1 = M1M2[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);
      p_chain = M1M2.begin() + r1;
      M1M2.erase(p_chain);
      break;
    case 48 :                              //IM2M1*+IM2*=POLY
      species[5] -= 1;
      species[8] -= 1;
      species[10] += 2;
      r1 = my_rand(M2M1.size());
      c_tab1 = M2M1[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);
      p_chain = M2M1.begin() + r1;
      M2M1.erase(p_chain);
      break;
    case 49 :                              //IM2M2*+IM2*=POLY
      species[5] -= 1;
      species[9] -= 1;
      species[10] += 2;
      r1 = my_rand(M2M2.size());
      c_tab1 = M2M2[r1];
      poly.chain1.push_back(c_tab1);
      poly.chain2.push_back(-1);
      for (int w = 0; w < allchains[c_tab1].M1.size(); w++) { M1_counter += allchains[c_tab1].M1[w]; };
      for (int w = 0; w < allchains[c_tab1].M2.size(); w++) { M2_counter += allchains[c_tab1].M2[w]; };
      tot_dyad[0] += allchains[c_tab1].dyad[0];        //AA
      tot_dyad[1] += allchains[c_tab1].dyad[1];        //AB
      tot_dyad[2] += allchains[c_tab1].dyad[2];        //BA
      tot_dyad[3] += allchains[c_tab1].dyad[3];        //BB
      tot_triad[0] += allchains[c_tab1].triad[0];    // AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB
      tot_triad[1] += allchains[c_tab1].triad[1];
      tot_triad[2] += allchains[c_tab1].triad[2];
      tot_triad[3] += allchains[c_tab1].triad[3];
      tot_triad[4] += allchains[c_tab1].triad[4];
      tot_triad[5] += allchains[c_tab1].triad[5];
      tot_triad[6] += allchains[c_tab1].triad[6];
      tot_triad[7] += allchains[c_tab1].triad[7];
      poly.num_M1.push_back(M1_counter);
      poly.num_M2.push_back(M2_counter);
      num_chain = num_chain + 1;
      allchains.push_back(new_chain);
      allchains[allchains.size() - 1].M1.push_back(1);
      poly.chain1.push_back(num_chain);
      poly.chain2.push_back(-1);
      poly.num_M1.push_back(0);
      poly.num_M2.push_back(1);
      p_chain = M2M2.begin() + r1;
      M2M2.erase(p_chain);
      break;
  }
}
