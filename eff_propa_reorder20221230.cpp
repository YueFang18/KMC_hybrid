/*
 * combine KMC and ODE to calculate scale_factor
 */

#include<cstdio>
#include<iomanip>
#include<stdlib.h>
#include<vector>
#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include "my_rand.h"    // header file for random number generation
#include "efficient_explicit_sequence_record_number.h"  // header file for species and sequence update
#include "molecular_weight.h"
#include "ode.h"
#include <utility>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>

using namespace std;


int main( ) {

    srand(time(NULL));

    clock_t start,end;
    start=clock();

    std::string data_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/GAO_project/P1_KMC_Acceleration_Scaling/code/effPropaRecorder0719/src/";
    std::string out_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/GAO_project/P1_KMC_Acceleration_Scaling/code/effPropaRecorder0719/output/";

    //............. Declaration of varibles begins.....//
    const double Na = 6.022e23;
    const double R = 0.001986;    //kcal/mol/K
    double mol_den[2];
    double mol_wt[2];
    double mol_frac[2], sol_mf = 0.0, mf_denm = 0.0, mf_monom = 0.0;
    double temp;        //Tempetaure in K
    double init_conc;   //mol/lit
    double solvent_conc;
    double time_limit;
    double M0;
    int reaction_index;
    int num_chain = -1;     //chain counter, every time a new chain added, make sure to assign an end unit to it
    double tot_dyad[4] = {0};
    double tot_dyad_tot = 0;
    double dyad_frac[4] = {0.0};
    double tot_triad[8] = {0};
    double tot_triad_tot = 0;
    double triad_frac[8] = {0.0};
    double reaction_counter = 0.0, tau = 0.0;
    double total_rate = 0.0, reaction = 0.0, time = 0.0, init_vol = 0.0, out_time = 0.0,d_time = 0.0,D_time = 0.0,Aver_d_time = 0.0;
    double rand_11, rand_22;
    double species[12] = {0};
    double all_species =0.0;
    double add_spe3=0.0;
    double species_all = 0;
    double Area = 0;
    double all_scale_fac = 0;
    double Nscale = 0;
    char str_in;
    int c_tab11, c_tab22=0;
    double aver_spe3=0;

    double rate[57] = {0.0};
    double A[57] = {0.0}, E[57] = {0.0}, c[57] = {0.0};
    double conversion[2] = {0.0}, over_x = 0.0, over_den = 0.0, over_num = 0.0;
    double volume = 0.0, species_sum = 0.0, M_initial[2] = {0.0}, num_monomer[2] = {0.0};
    double Mn = 0.0, Mw = 0.0, pdi = 0.0, Dpn = 0.0, M_total = 0.0, M_total_square = 0.0;
    double M_start[2] = {0.0}, M_present[2] = {0.0}, f_inst[2] = {0.0}, F_inst[2] = {0.0}, M_sum = 0.0, M_diff = 0.0;
    double num_M[2] = {0.0}, f[2] = {0.0}, F[2] = {0.0}, num_M_sum = 0;
    const int ccd_size = 1000;
    double ccd[ccd_size][ccd_size];
    int rxn_matrix[57][11] = {0}; //row column
    int dep_graph[57][57] = {0};
    double dyad_frac_int[4] = {0};
    double over_x_int = 0;
    double Mn_int = 0;
    double sum_sq_025 = 0;
    double min_sum_sq_025 = 100;
    double scale_fac = 1.0;
    double radical_ode = 0.0,Aver_add_spe = 0.0;
    double x[11]={0.0};
    double start_t=0.0, end_t=5000.0, dt_t=5.0;
    double radical[11] = {0.0};
//    vector_type r(11);
//    r[0]=0.005;
//    r[1]=0.0;
//    r[2]=1.50499;
//    r[3]=0.644795;
//    r[4]=0;

    vector_type r(11);
    r[0]=0.005;
    r[2]=1.64617;
    r[3]=1.34687;
    vector_type Rad(11);
    Rad[0]=0.005;
    Rad[2]=1.64617;
    Rad[3]=1.34687;

    stiff_system stiff_sys;
    stiff_system_jacobi system_jacobi;

// declaration of flags
    int fl_conv = 0, fl_temp = 0, fl_init = 0, fl_mono = 0,fl_aver=0;
    int fl_mwd = 0, fl_ads = 0, fl_eq = 0;
    double conv_in, temp_in, init_in, mono1_in, mono2_in;
    double solv_mol = 0;

    double all_species_counter=0;
    double species_counter=0;
    for (int i = 0; i < ccd_size; i++) {
        for (int j = 0; j < ccd_size; j++) {
            ccd[i][j] = 0;
        }
    }
    //............. Declaration of varibles ends......//
    //............. Assigning Input  data to the varibles from input.txt, Begins ......//


    ifstream input0(data_dir + "input.txt");
    input0 >> scale_fac;
    input0 >> temp;
    for (int i = 0; i < 2; i++)
        input0 >> mol_den[i]; //7.01  9.4
    for (int i = 0; i < 2; i++)
        input0 >> mol_wt[i];  //128  100.0
    for (int i = 0; i < 2; i++)
        input0 >> mol_frac[i];  //0.165  0.135
    input0 >> init_conc;        //0.005
    input0 >> solvent_conc;    //11.23
    input0 >> time_limit;      //30000
    input0 >> M0;              //1e10
    for (int i = 0; i < 57; i++) {
        input0 >> A[i] >> E[i];
        c[i] = A[i] * exp(-E[i] / (R * temp));
    }

    input0.close();

    //..............defining dependency graph......................//
    ifstream file(data_dir + "dependency_graph.txt");
    for (int i = 0; i < 57; i++)
        for (int j = 0; j <= 10; j++) {
            file >> rxn_matrix[i][j];
        }
    file.close();

    // to figure out what reaction was happened
    for (int i = 0; i < 57; i++)
        for (int j = 0; j <= 10; j++)
            for (int k = 0; k < 57; k++) {
                if ((rxn_matrix[i][j] == 1 || rxn_matrix[i][j] == -1) && (rxn_matrix[k][j] == -1 || rxn_matrix[k][j] == 2)) {
                    dep_graph[i][k] = 1;
                }
            }


//..............defining dependency graph..end....................//
//............. Assigning Input  data to the varibles from input.txt, Ends ......//
//..............reading flags from user input.....................//
//..............reading flags from user input...end..................//

    ofstream conv(out_dir + "Hy_conversion.out");
    ofstream molwt(out_dir+"Hy_molecularweight.out");
    ofstream fract(out_dir+"Hy_fractions.out");
    ofstream speciesout(out_dir+"Hy_species.out");
    ofstream species_con(out_dir+"Hy_spe_con.out");
    ofstream tim(out_dir+"Hy_time.out");
    ofstream sf(out_dir+"Hy_SF.out");
    ofstream m1m1(out_dir+"Hy_M1M1.out");
    ofstream polyinfo(out_dir+"Hy_polychain.out");
    ofstream polyinfomw(out_dir+"Hy_polychain_mw.out");
    ofstream aver(out_dir + "Hy_aver_spe3.out");

    ofstream ads_dis(out_dir+"Hy_ads_distribution.out");
    ofstream triad_dis(out_dir+"Hy_triad_distribution.out");
    ofstream eq_dyad(out_dir+"Hy_equal_dyad_point.out");
    ofstream ccd_count1(out_dir + "Hy_ccd_count1.out");
    ofstream sequence(out_dir + "Hy_chianinfo.out");
    ccd_count1 << "\t" << "ccd" << endl;
    conv << "time" << "\t" << "conversion(s)" << "\t" << "overall_conversion" << endl;
    molwt << "time" << "\t" << "\t" << "Mn" << "\t" << "\t" << "Mw" << "\t" << "\t" << "pdi" << "\t" << "\t" << "Dpn"
          << endl;
    fract << "time" << "\t" << "\t" << "f" << "\t" << "\t" << "\t" << "\t" << "F" << "\t" << "\t" << "\t" << "\t"
          << "f_inst" << "\t" << "\t" << "\t" << "\t" << "F_inst" << endl;

    for (int i = 0; i < 2; i++) {
        init_vol += mol_frac[i] * M0 / (mol_den[i] * Na);  //Volume in liters v1+v2
        mf_denm += mol_den[i];
        mf_monom += mol_frac[i];
    }

    sol_mf = 1 - mf_monom;
    init_vol += sol_mf * M0 / (solvent_conc * Na); //V=v1+v2+vsol

    volume = init_vol;                       //total volume
    cout << volume << endl;

    //........................number of molecules of the species  ......//

    species[0] = round(init_conc * init_vol * Na + 0.5); //add0.5??

    species[2] = round(M0 * mol_frac[0]);
    species[3] = round(M0 * mol_frac[1]);
    species[11] = round(M0 * sol_mf);


    cout << species[0]/(volume * Na) <<"\t"<<species[2]/(volume * Na) <<"\t"<<species[3]/(volume * Na)<<"\t"<<species[11]/(volume * Na)<< endl;
    cout << species[0] <<"\t"<<species[2] <<"\t"<<species[3]<<"\t"<<species[11]<<endl;

    M_initial[0] = species[2];
    M_initial[1] = species[3];
    M_start[0] = species[2];
    M_start[1] = species[3];

    //initialize the propensity functions//
    rate[ 50 ]=c[ 0  ]*species[0];
    rate[ 51 ]=1/scale_fac*c[ 1  ]/(Na*volume)*species[1]*species[ 2 ];
    rate[ 52 ]=1/scale_fac*c[ 2  ]/(Na*volume)*species[1]*species[ 3 ];
    rate[ 53 ]=1/scale_fac*c[ 3  ]/(Na*volume)*species[ 4 ]*species[ 2 ];
    rate[ 54 ]=1/scale_fac*c[ 4  ]/(Na*volume)*species[ 4 ]*species[ 3 ];
    rate[ 55 ]=1/scale_fac*c[ 5  ]/(Na*volume)*species[ 5 ]*species[ 2 ];
    rate[ 56 ]=1/scale_fac*c[ 6  ]/(Na*volume)*species[ 5 ]*species[ 3 ];
    rate[ 0  ]=	1/scale_fac*	c[ 7  ]/(Na*volume)*species[ 2 ]*species[ 6 ];
    rate[ 1  ]=	1/scale_fac*	c[ 8  ]/(Na*volume)*species[ 2 ]*species[ 7 ];
    rate[ 2  ]=	1/scale_fac*	c[ 9  ]/(Na*volume)*species[ 2 ]*species[ 8 ];
    rate[ 3  ]=	1/scale_fac*	c[ 10 ]/(Na*volume)*species[ 2 ]*species[ 9 ];
    rate[ 4  ]=	1/scale_fac*	c[ 11 ]/(Na*volume)*species[ 3 ]*species[ 6 ];
    rate[ 5  ]=	1/scale_fac*	c[ 12 ]/(Na*volume)*species[ 3 ]*species[ 7 ];
    rate[ 6  ]=	1/scale_fac*	c[ 13 ]/(Na*volume)*species[ 3 ]*species[ 8 ];
    rate[ 7  ]=	1/scale_fac*	c[ 14 ]/(Na*volume)*species[ 3 ]*species[ 9 ];
    rate[ 8  ]=	1/scale_fac/scale_fac	*	c[ 15 ]*2/(Na*volume)*species[ 6 ]*(species[ 6 ]-1)/2;
    rate[ 9  ]=	1/scale_fac/scale_fac	*	c[ 16 ]/(Na*volume)*species[ 6 ]*species[ 7 ];
    rate[ 10 ]=	1/scale_fac/scale_fac	*	c[ 17 ]/(Na*volume)*species[ 6 ]*species[ 8 ];
    rate[ 11 ]=	1/scale_fac/scale_fac	*	c[ 18 ]/(Na*volume)*species[ 6 ]*species[ 9 ];
    rate[ 12 ]=	1/scale_fac/scale_fac	*	c[ 19 ]*2/(Na*volume)*species[ 7 ]*(species[ 7 ]-1)/2;
    rate[ 13 ]=	1/scale_fac/scale_fac	*	c[ 20 ]/(Na*volume)*species[ 7 ]*species[ 8 ];
    rate[ 14 ]=	1/scale_fac/scale_fac	*	c[ 21 ]/(Na*volume)*species[ 7 ]*species[ 9 ];
    rate[ 15 ]=	1/scale_fac/scale_fac	*	c[ 22 ]*2/(Na*volume)*species[ 8 ]*(species[ 8 ]-1)/2;
    rate[ 16 ]=	1/scale_fac/scale_fac	*	c[ 23 ]/(Na*volume)*species[ 8 ]*species[ 9 ];
    rate[ 17 ]=	1/scale_fac/scale_fac	*	c[ 24 ]*2/(Na*volume)*species[ 9 ]*(species[ 9 ]-1)/2;
    rate[ 18 ]=	1/scale_fac/scale_fac	*	c[ 25 ]*2/(Na*volume)*species[ 4 ]*(species[ 4 ]-1)/2;
    rate[ 19 ]=	1/scale_fac/scale_fac	*	c[ 26 ]/(Na*volume)*species[ 4 ]*species[ 5 ];
    rate[ 20 ]=	1/scale_fac/scale_fac	*	c[ 27 ]*2/(Na*volume)*species[ 5 ]*(species[ 5 ]-1)/2;
    rate[ 21 ]=	1/scale_fac/scale_fac	*	c[ 28 ]/(Na*volume)*species[ 4 ]*species[ 6 ];
    rate[ 22 ]=	1/scale_fac/scale_fac	*	c[ 29 ]/(Na*volume)*species[ 4 ]*species[ 7 ];
    rate[ 23 ]=	1/scale_fac/scale_fac	*	c[ 30 ]/(Na*volume)*species[ 4 ]*species[ 8 ];
    rate[ 24 ]=	1/scale_fac/scale_fac	*	c[ 31 ]/(Na*volume)*species[ 4 ]*species[ 9 ];
    rate[ 25 ]=	1/scale_fac/scale_fac	*	c[ 32 ]/(Na*volume)*species[ 5 ]*species[ 6 ];
    rate[ 26 ]=	1/scale_fac/scale_fac	*	c[ 33 ]/(Na*volume)*species[ 5 ]*species[ 7 ];
    rate[ 27 ]=	1/scale_fac/scale_fac	*	c[ 34 ]/(Na*volume)*species[ 5 ]*species[ 8 ];
    rate[ 28 ]=	1/scale_fac/scale_fac	*	c[ 35 ]/(Na*volume)*species[ 5 ]*species[ 9 ];
    rate[ 29 ]=	1/scale_fac/scale_fac	*	c[ 36 ]*2/(Na*volume)*species[ 6 ]*(species[ 6 ]-1)/2;
    rate[ 30 ]=	1/scale_fac/scale_fac	*	c[ 37 ]/(Na*volume)*species[ 6 ]*species[ 7 ];
    rate[ 31 ]=	1/scale_fac/scale_fac	*	c[ 38 ]/(Na*volume)*species[ 6 ]*species[ 8 ];
    rate[ 32 ]=	1/scale_fac/scale_fac	*	c[ 39 ]/(Na*volume)*species[ 6 ]*species[ 9 ];
    rate[ 33 ]=	1/scale_fac/scale_fac	*	c[ 40 ]*2/(Na*volume)*species[ 7 ]*(species[ 7 ]-1)/2;
    rate[ 34 ]=	1/scale_fac/scale_fac	*	c[ 41 ]/(Na*volume)*species[ 7 ]*species[ 8 ];
    rate[ 35 ]=	1/scale_fac/scale_fac	*	c[ 42 ]/(Na*volume)*species[ 7 ]*species[ 9 ];
    rate[ 36 ]=	1/scale_fac/scale_fac	*	c[ 43 ]*2/(Na*volume)*species[ 8 ]*(species[ 8 ]-1)/2;
    rate[ 37 ]=	1/scale_fac/scale_fac	*	c[ 44 ]/(Na*volume)*species[ 8 ]*species[ 9 ];
    rate[ 38 ]=	1/scale_fac/scale_fac	*	c[ 45 ]*2/(Na*volume)*species[ 9 ]*(species[ 9 ]-1)/2;
    rate[ 39 ]=	1/scale_fac/scale_fac	*	c[ 46 ]*2/(Na*volume)*species[ 4 ]*(species[ 4 ]-1)/2;
    rate[ 40 ]=	1/scale_fac/scale_fac	*	c[ 47 ]/(Na*volume)*species[ 4 ]*species[ 5 ];
    rate[ 41 ]=	1/scale_fac/scale_fac	*	c[ 48 ]*2/(Na*volume)*species[ 5 ]*(species[ 5 ]-1)/2;
    rate[ 42 ]=	1/scale_fac/scale_fac	*	c[ 49 ]/(Na*volume)*species[ 4 ]*species[ 6 ];
    rate[ 43 ]=	1/scale_fac/scale_fac	*	c[ 50 ]/(Na*volume)*species[ 4 ]*species[ 7 ];
    rate[ 44 ]=	1/scale_fac/scale_fac	*	c[ 51 ]/(Na*volume)*species[ 4 ]*species[ 8 ];
    rate[ 45 ]=	1/scale_fac/scale_fac	*	c[ 52 ]/(Na*volume)*species[ 4 ]*species[ 9 ];
    rate[ 46 ]=	1/scale_fac/scale_fac	*	c[ 53 ]/(Na*volume)*species[ 5 ]*species[ 6 ];
    rate[ 47 ]=	1/scale_fac/scale_fac	*	c[ 54 ]/(Na*volume)*species[ 5 ]*species[ 7 ];
    rate[ 48 ]=	1/scale_fac/scale_fac	*	c[ 55 ]/(Na*volume)*species[ 5 ]*species[ 8 ];
    rate[ 49 ]=	1/scale_fac/scale_fac	*	c[ 56 ]/(Na*volume)*species[ 5 ]*species[ 9 ];




    //................................... KMC Loop Begins .............//

    int flag1=0, Count_n=0,Count_n1=0;

    while (time <= time_limit) {

        rand_11 =(float) rand() / RAND_MAX;
        rand_22 = (float) rand() / RAND_MAX;

//        rand_11 = 1.0 * (my_rand(RAND_MAX - 1) + 1.0) / RAND_MAX;      //random number for next rxn time step
//        rand_22 = 1.0 * (my_rand(RAND_MAX - 1) + 1.0) / RAND_MAX;      //random number for event excution


        total_rate = 0.0;

        for (int i = 0; i < 57; i++) {
            total_rate = total_rate + rate[i];
        }

        tau = (1 / total_rate) * log(1/rand_11);      // reaction time step

        reaction = rand_22 * total_rate;
        reaction_counter = rate[0];

        reaction_index = 0;

        while (reaction_counter < reaction)       // condition for which reaction to be excuted
        {
            reaction_index++;
            reaction_counter += rate[reaction_index];
        }

        over_den = 0.0;
        over_num = 0.0;

        time += tau;       // reaction time, sec
        out_time += tau;
        d_time += tau;
        D_time += tau;

        all_species_counter += (species[1]+species[4]+species[5]+species[6]+species[7]+species[8]+species[9])*tau;
        species_counter++;

        do{
            all_species = species[1] + species[4] + species[5] + species[6] + species[7] + species[8] + species[9];
            Aver_add_spe += all_species * tau;
            Aver_d_time += tau;
        }while (D_time >80 && Aver_add_spe/Aver_d_time < 2);


//        do{
//            all_species = species[1] + species[4] + species[5] + species[6] + species[7] + species[8] + species[9];
//            Aver_add_spe += all_species * tau;
//            Aver_d_time += tau;
//            if (Aver_add_spe/Aver_d_time>2){break;}
//        } while(D_time > 95);


//calculate the SF at beginning with 2
        if ( scale_fac == 1 && d_time > 0.05 && flag1 == 0) {

            ode(stiff_sys, system_jacobi,r,Rad,d_time);
            radical_ode = r[1]+r[4]+r[5]+r[6]+r[7]+r[8]+r[9];
            scale_fac = 2 / (Na * volume * radical_ode);
            flag1 =1;
        }


//setting for re-update with ODE
        if ( D_time > 100 &&  Aver_add_spe !=0 && Aver_d_time !=0 ){

            cout<< Aver_add_spe<<"  "<<Aver_d_time<<endl;

            aver_spe3 = Aver_add_spe/Aver_d_time;
            aver<<aver_spe3<<endl;

            ode(stiff_sys, system_jacobi,r,Rad,D_time);
            radical_ode = r[1]+r[4]+r[5]+r[6]+r[7]+r[8]+r[9];
            scale_fac = aver_spe3 / (Na * volume * radical_ode);
            sf<<"   "<< time <<"   "<< Aver_add_spe / Aver_d_time <<"   "<< scale_fac <<"   "<< radical_ode <<"    "<<r[2]<<"    "<<r[3]<<"    "<<"222   222"<<"    "<< rand_11 <<"    "<< rand_22 <<endl;

            D_time = 0.0;
            Aver_add_spe = 0;
            Aver_d_time = 0;
            cout<<"111          111"<<"    "<<rand_11<<"    "<<rand_22<<"    "<<tau<<endl;
        }


        if (out_time >= 100.0) {
//            sf<<Aver_add_spe/Aver_d_time<<"   "<<scale_fac<<"   "<<all_species/(Na * volume) <<"   "<<radical_ode<<endl;

            speciesout << time;
            species_con << time;
            out_time = 0.0;

            for (int i = 0; i < 11; i++) {
                speciesout << "    " << species[i];
                species_con << " " << species[i] / (Na * volume);
            }
            speciesout << endl;
            species_con << endl;


            for (int i = 0; i < 2; i++) {
                conversion[i] = num_monomer[i] / M_initial[i];
                M_sum += M_start[i] + M_present[i];
                M_diff += M_start[i] - M_present[i];
                over_num += num_monomer[i];
                over_den += M_initial[i];
            }

            over_x = over_num / over_den; //overall conversion
            conv << time << "\t" << conversion[0] << "\t" << conversion[1] << "\t" << over_x << endl;
            f_inst[0] = (M_present[0] + M_start[0]) / (M_sum);
            F_inst[0] = (M_start[0] - M_present[0]) / (M_diff);
            f_inst[1] = (M_present[1] + M_start[1]) / (M_sum);
            F_inst[1] = (M_start[1] - M_present[1]) / (M_diff);
            M_start[0] = M_present[0];
            M_start[1] = M_present[1];
            num_M[0] = 0;
            num_M[1] = 0;
            num_M_sum = 0.0;
            species_sum = 0.0;


            if (fl_mwd == 0) {
//         molecular_weight(time, mol_wt, Mn, Mw, num_M, poly); //function for calculating molecular weight
                molecular_weight(mol_wt,
                                 Mn,
                                 Mw,
                                 num_M,
                                 species,
                                 scale_fac,
                                 poly); //function for calculating molecular weight

                pdi = Mw / Mn;
                for (int i = 0; i < 2; i++) {
                    species_sum += species[2 + i];
                    num_M_sum += num_M[i];
                }
                for (int i = 0; i < 2; i++) {
                    f[i] = species[2 + i] / species_sum;   //avg compostion of monomers in the feed
                    F[i] = num_M[i] / (num_M_sum);       //avg compostion of monomers in the dead end polymer chains
                }
                Dpn = num_M_sum / (poly.chain1.size());    // degree of polymerization
                molwt << time << "\t" << Mn << "\t" << "\t" << Mw << "\t" << "\t" << pdi << "\t" << "\t" << Dpn << endl;
                fract << time << "\t" << f[0] << "\t" << f[1] << "\t" << "\t" << F[0] << "\t" << F[1] << "\t" << "\t"
                      << f_inst[0] << "\t" << f_inst[1] << "\t" << "\t" << F_inst[0] << "\t" << F_inst[1] << endl;
            }

            if (fl_ads == 1) {
                tot_dyad_tot = tot_dyad[0] + tot_dyad[1] + tot_dyad[2] + tot_dyad[3];
                ads_dis << over_x << "\t" << "\t";
                sum_sq_025 = 0;
                for (int i = 0; i < 4; i++) {
                    dyad_frac[i] = tot_dyad[i] / tot_dyad_tot;
                    ads_dis << dyad_frac[i] << "\t" << "\t";
                    if (fl_eq == 0) {
                        sum_sq_025 += (dyad_frac[i] - 0.25) * (dyad_frac[i] - 0.25);
                    }
                }
                ads_dis << endl;
                if (fl_eq == 0) {
                    if (sum_sq_025 < min_sum_sq_025) {
                        min_sum_sq_025 = sum_sq_025;
                        for (int i = 0; i < 4; i++) { dyad_frac_int[i] = dyad_frac[i]; }
                        over_x_int = over_x;
                        Mn_int = Mn;
                    }
                }
                triad_dis << over_x << "\t";
                tot_triad_tot = 0;
                for (int i = 0; i < 8; i++) { tot_triad_tot += tot_triad[i]; }
                for (int i = 0; i < 8; i++) {
                    triad_frac[i] = tot_triad[i] / tot_triad_tot;
                    triad_dis << triad_frac[i] << "\t";
                }
                triad_dis << endl;
            }
        }


        efficient_explicit_sequence_record_number(reaction_index,
                                                  num_monomer,
                                                  species,
                                                  num_chain,
                                                  M1M1,
                                                  M1M2,
                                                  M2M1,
                                                  M2M2,
                                                  allchains,
                                                  poly,
                                                  tot_dyad,
                                                  tot_triad);  //function for updating species info when event occurs



        all_scale_fac += scale_fac * tau;

        // update rates
        if (dep_graph[reaction_index][50]==1) rate[50]=c[0]*species[0];
        if (dep_graph[reaction_index][51]==1) rate[51]=1/scale_fac*c[1]/(Na*volume)*species[1]*species[2];
        if (dep_graph[reaction_index][52]==1) rate[52]=1/scale_fac*c[2]/(Na*volume)*species[1]*species[3];
        if (dep_graph[reaction_index][53]==1) rate[53]=1/scale_fac*c[3]/(Na*volume)*species[4]*species[2];
        if (dep_graph[reaction_index][54]==1) rate[54]=1/scale_fac*c[4]/(Na*volume)*species[4]*species[3];
        if (dep_graph[reaction_index][55]==1) rate[55]=1/scale_fac*c[5]/(Na*volume)*species[5]*species[2];
        if (dep_graph[reaction_index][56]==1) rate[56]=1/scale_fac*c[6]/(Na*volume)*species[5]*species[3];
        if (dep_graph[reaction_index][0]==1) rate[0]=1/scale_fac*c[7]/(Na*volume)*species[2]*species[6];
        if (dep_graph[reaction_index][1]==1) rate[1]=1/scale_fac*c[8]/(Na*volume)*species[2]*species[7];
        if (dep_graph[reaction_index][2]==1) rate[2]=1/scale_fac*c[9]/(Na*volume)*species[2]*species[8];
        if (dep_graph[reaction_index][3]==1) rate[3]=1/scale_fac*c[10]/(Na*volume)*species[2]*species[9];
        if (dep_graph[reaction_index][4]==1) rate[4]=1/scale_fac*c[11]/(Na*volume)*species[3]*species[6];
        if (dep_graph[reaction_index][5]==1) rate[5]=1/scale_fac*c[12]/(Na*volume)*species[3]*species[7];
        if (dep_graph[reaction_index][6]==1) rate[6]=1/scale_fac*c[13]/(Na*volume)*species[3]*species[8];
        if (dep_graph[reaction_index][7]==1) rate[7]=1/scale_fac*c[14]/(Na*volume)*species[3]*species[9];
        if (dep_graph[reaction_index][8]==1) rate[8]=1/scale_fac/scale_fac*c[15]*2/(Na*volume)*species[6]*(species[6]-1)/2;
        if (dep_graph[reaction_index][9]==1) rate[9]=1/scale_fac/scale_fac*c[16]/(Na*volume)*species[6]*species[7];
        if (dep_graph[reaction_index][10]==1) rate[10]=1/scale_fac/scale_fac*c[17]/(Na*volume)*species[6]*species[8];
        if (dep_graph[reaction_index][11]==1) rate[11]=1/scale_fac/scale_fac*c[18]/(Na*volume)*species[6]*species[9];
        if (dep_graph[reaction_index][12]==1) rate[12]=1/scale_fac/scale_fac*c[19]*2/(Na*volume)*species[7]*(species[7]-1)/2;
        if (dep_graph[reaction_index][13]==1) rate[13]=1/scale_fac/scale_fac*c[20]/(Na*volume)*species[7]*species[8];
        if (dep_graph[reaction_index][14]==1) rate[14]=1/scale_fac/scale_fac*c[21]/(Na*volume)*species[7]*species[9];
        if (dep_graph[reaction_index][15]==1) rate[15]=1/scale_fac/scale_fac*c[22]*2/(Na*volume)*species[8]*(species[8]-1)/2;
        if (dep_graph[reaction_index][16]==1) rate[16]=1/scale_fac/scale_fac*c[23]/(Na*volume)*species[8]*species[9];
        if (dep_graph[reaction_index][17]==1) rate[17]=1/scale_fac/scale_fac*c[24]*2/(Na*volume)*species[9]*(species[9]-1)/2;
        if (dep_graph[reaction_index][18]==1) rate[18]=1/scale_fac/scale_fac*c[25]*2/(Na*volume)*species[4]*(species[4]-1)/2;
        if (dep_graph[reaction_index][19]==1) rate[19]=1/scale_fac/scale_fac*c[26]/(Na*volume)*species[4]*species[5];
        if (dep_graph[reaction_index][20]==1) rate[20]=1/scale_fac/scale_fac*c[27]*2/(Na*volume)*species[5]*(species[5]-1)/2;
        if (dep_graph[reaction_index][21]==1) rate[21]=1/scale_fac/scale_fac*c[28]/(Na*volume)*species[4]*species[6];
        if (dep_graph[reaction_index][22]==1) rate[22]=1/scale_fac/scale_fac*c[29]/(Na*volume)*species[4]*species[7];
        if (dep_graph[reaction_index][23]==1) rate[23]=1/scale_fac/scale_fac*c[30]/(Na*volume)*species[4]*species[8];
        if (dep_graph[reaction_index][24]==1) rate[24]=1/scale_fac/scale_fac*c[31]/(Na*volume)*species[4]*species[9];
        if (dep_graph[reaction_index][25]==1) rate[25]=1/scale_fac/scale_fac*c[32]/(Na*volume)*species[5]*species[6];
        if (dep_graph[reaction_index][26]==1) rate[26]=1/scale_fac/scale_fac*c[33]/(Na*volume)*species[5]*species[7];
        if (dep_graph[reaction_index][27]==1) rate[27]=1/scale_fac/scale_fac*c[34]/(Na*volume)*species[5]*species[8];
        if (dep_graph[reaction_index][28]==1) rate[28]=1/scale_fac/scale_fac*c[35]/(Na*volume)*species[5]*species[9];
        if (dep_graph[reaction_index][29]==1) rate[29]=1/scale_fac/scale_fac*c[36]*2/(Na*volume)*species[6]*(species[6]-1)/2;
        if (dep_graph[reaction_index][30]==1) rate[30]=1/scale_fac/scale_fac*c[37]/(Na*volume)*species[6]*species[7];
        if (dep_graph[reaction_index][31]==1) rate[31]=1/scale_fac/scale_fac*c[38]/(Na*volume)*species[6]*species[8];
        if (dep_graph[reaction_index][32]==1) rate[32]=1/scale_fac/scale_fac*c[39]/(Na*volume)*species[6]*species[9];
        if (dep_graph[reaction_index][33]==1) rate[33]=1/scale_fac/scale_fac*c[40]*2/(Na*volume)*species[7]*(species[7]-1)/2;
        if (dep_graph[reaction_index][34]==1) rate[34]=1/scale_fac/scale_fac*c[41]/(Na*volume)*species[7]*species[8];
        if (dep_graph[reaction_index][35]==1) rate[35]=1/scale_fac/scale_fac*c[42]/(Na*volume)*species[7]*species[9];
        if (dep_graph[reaction_index][36]==1) rate[36]=1/scale_fac/scale_fac*c[43]*2/(Na*volume)*species[8]*(species[8]-1)/2;
        if (dep_graph[reaction_index][37]==1) rate[37]=1/scale_fac/scale_fac*c[44]/(Na*volume)*species[8]*species[9];
        if (dep_graph[reaction_index][38]==1) rate[38]=1/scale_fac/scale_fac*c[45]*2/(Na*volume)*species[9]*(species[9]-1)/2;
        if (dep_graph[reaction_index][39]==1) rate[39]=1/scale_fac/scale_fac*c[46]*2/(Na*volume)*species[4]*(species[4]-1)/2;
        if (dep_graph[reaction_index][40]==1) rate[40]=1/scale_fac/scale_fac*c[47]/(Na*volume)*species[4]*species[5];
        if (dep_graph[reaction_index][41]==1) rate[41]=1/scale_fac/scale_fac*c[48]*2/(Na*volume)*species[5]*(species[5]-1)/2;
        if (dep_graph[reaction_index][42]==1) rate[42]=1/scale_fac/scale_fac*c[49]/(Na*volume)*species[4]*species[6];
        if (dep_graph[reaction_index][43]==1) rate[43]=1/scale_fac/scale_fac*c[50]/(Na*volume)*species[4]*species[7];
        if (dep_graph[reaction_index][44]==1) rate[44]=1/scale_fac/scale_fac*c[51]/(Na*volume)*species[4]*species[8];
        if (dep_graph[reaction_index][45]==1) rate[45]=1/scale_fac/scale_fac*c[52]/(Na*volume)*species[4]*species[9];
        if (dep_graph[reaction_index][46]==1) rate[46]=1/scale_fac/scale_fac*c[53]/(Na*volume)*species[5]*species[6];
        if (dep_graph[reaction_index][47]==1) rate[47]=1/scale_fac/scale_fac*c[54]/(Na*volume)*species[5]*species[7];
        if (dep_graph[reaction_index][48]==1) rate[48]=1/scale_fac/scale_fac*c[55]/(Na*volume)*species[5]*species[8];
        if (dep_graph[reaction_index][49]==1) rate[49]=1/scale_fac/scale_fac*c[56]/(Na*volume)*species[5]*species[9];



        M_present[0] = species[2];
        M_present[1] = species[3];
        M_sum = 0.0;
        M_diff = 0.0;

        // Data analuysis (for conversion and composition etc) loop begins...............//

    }

    /*................................ poly_chain info..........................*/
    for ( int i = 0; i < M1M1.size(); i++){
        m1m1 << M1M1[i] << "; ";
    }
    m1m1<< endl;
    for ( int i = 0; i < M1M2.size(); i++){
        m1m1 << M1M2[i] << "; ";
    }
    m1m1<< endl;
    for ( int i = 0; i < M2M1.size(); i++){
        m1m1 << M2M1[i] << "; ";
    }
    m1m1<< endl;
    for ( int i = 0; i < M2M2.size(); i++){
        m1m1 << M2M2[i] << "; ";
    }
    m1m1<< endl;

    double all_mw=0.0;
    for ( int i = 0; i < poly.chain1.size(); i++){
        polyinfo << poly.chain1[i] << " "<< poly.chain2[i] << " "<< poly.num_M1[i] << " "<<poly.num_M2[i] << " "<<poly.num_M1[i]*mol_wt[0]+poly.num_M2[i]*mol_wt[1]<< endl;
        all_mw += poly.num_M1[i]*mol_wt[0]+poly.num_M2[i]*mol_wt[1];
    }




    /*................................ KMC Loop ends here..........................*/
    if (fl_mwd==1)
    {
//        molecular_weight(time,mol_wt,Mn,Mw,num_M, poly); // function for calculating molecular weight
        molecular_weight(mol_wt,
                         Mn,
                         Mw,
                         num_M,
                         species,
                         scale_fac,
                         poly); //function for calculating molecular weight

        pdi=Mw/Mn;
        for(int i=0;i< 2 ;i++)
        {
            species_sum+=species[2+i];
            num_M_sum+=num_M[i];
        }
        for(int i=0;i< 2 ;i++)
        {
            f[i]=species[2+i]/species_sum;   //avg compostion of monomers in the feed
            F[i]=num_M[i]/(num_M_sum);       //avg compostion of monomers in the dead end polymer chains
        }
        Dpn=num_M_sum/(poly.chain1.size());    // degree of polymerization
        molwt<<time<<"\t"<<Mn<<"\t"<<"\t"<<Mw<<"\t"<<"\t"<<pdi<<"\t"<<"\t"<<Dpn<<endl;
        fract<<time<<"\t"<< f[0]<<"\t"<<f[1] <<"\t"<<"\t"<< F[0]<<"\t"<<F[1] <<"\t"<<"\t"<< f_inst[0]<<"\t"<<f_inst[1] <<"\t"<<"\t"<< F_inst[0]<<"\t"<<F_inst[1] <<endl;
    }

    if (fl_ads==0)
    {
        tot_dyad_tot=tot_dyad[0]+tot_dyad[1]+tot_dyad[2]+tot_dyad[3];
        ads_dis<<over_x<<"\t"<<"\t";
        sum_sq_025=0;
        for(int i=0;i<4;i++)
        {	dyad_frac[i]=tot_dyad[i]/tot_dyad_tot;
            ads_dis<<dyad_frac[i]<<"\t"<<"\t";
            if (fl_eq==0)
            {
                sum_sq_025+=(dyad_frac[i]-0.25)*(dyad_frac[i]-0.25);
            }
        }
        ads_dis<<endl;
        if(fl_eq==1)
        {
            if (sum_sq_025<min_sum_sq_025)
            {	min_sum_sq_025=sum_sq_025;
                for(int i=0;i<4;i++) {dyad_frac_int[i]=dyad_frac[i];}
                over_x_int=over_x;
                Mn_int=Mn;
            }
        }
        triad_dis<<over_x<<"\t";
        tot_triad_tot=0;
        for (int i=0;i<8;i++){tot_triad_tot+=tot_triad[i];}
        for (int i=0;i<8;i++)
        {	triad_frac[i]=tot_triad[i]/tot_triad_tot;
            triad_dis<<triad_frac[i]<<"\t";
        }
        triad_dis<<endl;
    }

    if (fl_eq==0)
    {
        eq_dyad<<dyad_frac_int[0]<<"\t"<<dyad_frac_int[1]<<"\t"<<dyad_frac_int[2]<<"\t"<<dyad_frac_int[3]<<endl;
        eq_dyad<<over_x_int<<"\t"<<Mn_int<<endl;
    }





    /*................................ Recording sequence starts..........................*/

    end=clock();
    tim<<"average scaling_factor= "<< (all_scale_fac/time) <<endl;
    tim<<"runnning time = "<<(double)(end-start)/CLOCKS_PER_SEC<<"s"<<endl;
    tim<<"all_species_counter= "<<all_species_counter<<endl;
    tim<<"species_times_counter= "<<species_counter<<endl;
    tim<<"average_radical_species= "<<all_species_counter/time<<endl;


    return 0;
}