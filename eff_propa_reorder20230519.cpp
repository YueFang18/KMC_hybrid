/*
 * combine KMC and ODE to calculate scale_factor
 */
#include <random>
#include <cstdio>
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
#include <random>

using namespace std;
//int int_DX;


int main( ) {

    srand(time(NULL));

    clock_t start,end;
    start=clock();

    std::string data_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P1_KMC_Acceleration_Scaling/code/effATRP_new/src/";
    std::string out_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P1_KMC_Acceleration_Scaling/code/effATRP_new/output/data20230519/kd2.4e7/";

    //............. Declaration of varibles begins.....//
    const double Na = 6.022e23;
    const double R = 0.001986;    //kcal/mol/K
    double mol_den[2];
    double mol_wt[2];
    double mol_frac[2], sol_mf = 0.0, mf_denm = 0.0, mf_monom = 0.0;
    double temp;        //Tempetaure in K
    double init_conc[2];   //mol/lit
    double solvent_conc;
    double time_limit;
    double M0;
    int reaction_index,index0=0,index1=0,index2=0,index3=0;

    int Nscale=0;
    int num_chain = -1;     //chain counter, every time a new chain added, make sure to assign an end unit to it

    double reaction_counter = 0.0, tau = 0.0;
    double total_rate = 0.0, reaction = 0.0, time = 0.0, init_vol = 0.0, D_time = 0.0, d_time_0 = 0.0,d_time_1 = 0.0,pre_d_time = 0.0,pre_d_time_spe3 = 0.0,pre_d_time_spe2 = 0.0,aver_spe3=0.0,aver_radical=0.0;
    double rand_11,rand_22;
    double species[6] = {0};
    long double s3=0.0;
    double s5=0.0;
    double species_all = 0;
    double Area = 0;
    double all_scale_fac = 0;

    double rate[4] = {0.0};
    double A[4] = {0.0}, E[4] = {0.0}, c[4] = {0.0};
    double conversion[2] = {0.0}, over_x = 0.0, over_den = 0.0, over_num = 0.0;
    double volume = 0.0, species_sum = 0.0, M_initial[2] = {0.0}, num_monomer[2] = {0.0};
    double Mn = 0.0, Mw = 0.0, pdi = 0.0, Dpn = 0.0, M_total = 0.0, M_total_square = 0.0;
    double M_start[2] = {0.0}, M_present[2] = {0.0}, f_inst[2] = {0.0}, F_inst[2] = {0.0}, M_sum = 0.0, M_diff = 0.0;
    double num_M[2] = {0}, f[2] = {0.0}, F[2] = {0.0}, num_M_sum = 0;
    double num_Dr[1] = {0};
    double num_Pr[1] = {0};
    double num_P[1] = {0};
    const int ccd_size = 1000;
    double ccd[ccd_size][ccd_size];

    double over_x_int = 0;
    double Mn_int = 0;
    double sum_sq_025 = 0;
    double min_sum_sq_025 = 100;
    double scale_fac = 1.0;
    double scale_fac_s2 = 1.0;
    double radical_ode0 = 0.0,radical_ode1 = 0.0,radical_ode2 = 0.0,radical_ode3 = 0.0,radical_ode4 = 0.0,radical_ode5 = 0.0;
    double radical_kmc = 0.0;
    double x[6]={0.0};
    double start_t=0.0, end_t=5000.0, dt_t=5.0;
    double radical[6] = {0.0};
    double pre_add_spe3 = 0.0,add_spe2 = 0.0,pre_add_spe2 = 0.0;
    int N[319] = {0};
    vector_type r(6);
    r[0]=0.023;
    r[1]=0.0115;
    r[4]=4.67;
    vector_type Rad(6);
    Rad[0]=0.023;
    Rad[1]=0.0115;
    Rad[4]=4.67;
    poly_chain poly_test;
    stiff_system stiff_sys;
    stiff_system_jacobi system_jacobi;

// declaration of flags
    int fl_conv = 0, fl_temp = 0, fl_init = 0, fl_mono = 0,fl_aver = 0;
    int fl_mwd = 0, fl_ads = 0, fl_eq = 0;
    int flag1 = 0,flag0 = 0;
    double conv_in, temp_in, init_in, mono1_in, mono2_in;
    double solv_mol = 0;

//    int iter_counter=0;
    for (int i = 0; i < ccd_size; i++) {
        for (int j = 0; j < ccd_size; j++) {
            ccd[i][j] = 0;
        }
    }
    //............. Declaration of varibles ends......//
    //............. Assigning Input  data to the varibles from input.txt, Begins ......//

    ifstream input0(data_dir + "input0.txt");
    input0 >> scale_fac;
    input0 >> temp;
    for (int i = 0; i < 2; i++)
        input0 >> mol_den[i]; //7.01  9.4
    for (int i = 0; i < 2; i++)
        input0 >> mol_wt[i];  //128  100.0
    for (int i = 0; i < 2; i++)
        input0 >> mol_frac[i];  //0.58 0.0
    for (int i = 0; i < 2; i++)
        input0 >> init_conc[i];

    input0 >> solvent_conc;
    input0 >> time_limit;
    input0 >> M0;
    for (int i = 0; i < 4; i++) {
        input0 >> c[i];     //43.8; 2.4e7;1.65e3;2.82e7
        //c[i] = A[i] * exp(-E[i] / (R * temp));
    }
    input0.close();

    //..............defining dependency graph......................//



//............. Assigning Input  data to the varibles from input.txt, Ends ......//


    ofstream conv(out_dir + "Hy_conversion.out");
    ofstream molwt(out_dir+"Hy_molecularweight.out");
    ofstream speciesout(out_dir+"Hy_species.out");
    ofstream species_con(out_dir+"Hy_spe_con.out");
    ofstream wf(out_dir + "Hy_weight_frication.out");
    ofstream rat(out_dir + "Hy_rate.out");
    ofstream ode3(out_dir + "Hy_ode3.out");
    ofstream index(out_dir + "Hy_reaction_index.out");
    ofstream sf(out_dir + "Hy_scaling_factor.out");
    ofstream aver(out_dir + "Hy_aver_spe3.out");


    conv << "time" << "\t" << "conversion(s)" << "\t" << "overall_conversion" << endl;
    molwt << "time" << "\t" <<"conversion"<< "\t" << "Mn" << "\t" << "\t" << "Mw" << "\t"  << "pdi" << "\t"  << "Dpn"<< endl;
    ofstream fract(out_dir+"fractions.out");
    fract << "time" << "\t" << "\t" << "f" << "\t" << "\t" << "\t" << "\t" << "F" << "\t" << "\t" << "\t" << "\t"
          << "f_inst" << "\t" << "\t" << "\t" << "\t" << "F_inst" << endl;

    speciesout <<"time" << "\t"<<"species[0]" <<"\t"  <<"species[1] "<<"\t" <<"species[2]"<< "\t" <<"species[3] "<< "\t"<<" species[4]"<<  "\t"<<"species[5] " << endl;
    speciesout <<"time" << "\t"<<"Dr" <<"\t"  <<"C "<<"\t" <<"CX"<< "\t" <<"Pr* "<< "\t"<<" M"<<  "\t"<<"P" << endl;


    init_vol=(M0*mol_frac[0])/(Na*4.67); //4.67=species[4]/ (Na * volume);

    volume = init_vol;                   //total volume
    cout << volume << endl;

    //........................number of molecules of the species  ......//
    species[0] = round(init_conc[0] * init_vol * Na + 0.5); //add0.5??
    species[1] = round(init_conc[1]* init_vol * Na + 0.5);
    species[4] = round(4.67 * init_vol * Na + 0.5);

    vector<int> Dr(species[0],0);

    cout << species[0]/ (Na * volume)<< "\t" << species[1]/ (Na * volume)<< "\t" << species[2] / (Na * volume)<< "\t" << species[3]/ (Na * volume)<< "\t"<<species[4]/ (Na * volume) << "\t" << species[5]/ (Na * volume)<< "\t"<< endl;
    cout << species[0]<< "\t" << species[1]<< "\t" << species[2]<< "\t" << species[3]<< "\t"<<species[4] << "\t" << species[5]<< "\t"<< endl;

    M_initial[0] = species[4];
    M_start[0] = species[4];

    //initialize the propensity functions//
    rate[0] = c[0] * species[0] * species[1]/ (Na * volume);
    rate[1] = 1 / scale_fac /scale_fac_s2 *  c[1] / (Na * volume) * species[2] * species[3];
    rate[2] = 1 / scale_fac *  c[2] / (Na * volume) * species[3] * species[4];
    rate[3] = 1 / scale_fac /scale_fac * c[3] * 2 / (Na * volume) * species[3] * (species[3]-1) / 2;
//    rate[3] = 1  * c[3] * 2 / (Na * volume) * species[3] * (species[3]-1) / 2;

    //................................... KMC Loop Begins ......//
    int n=0;
    while (time <= time_limit) {

        rand_11 =(double) rand() / RAND_MAX;
        rand_22 = (double) rand() / RAND_MAX;

//        rand_11 = 1.0 * (my_rand(RAND_MAX - 1) + 1.0) / RAND_MAX;      //random number for next rxn time step
//        rand_22 = 1.0 * (my_rand(RAND_MAX - 1) + 1.0) / RAND_MAX;      //random number for event excution

        total_rate = 0.0;
        for (int i = 0; i < 4; i++) total_rate = total_rate + rate[i];
        tau = (1 / total_rate) * log (1 / rand_11);      // reaction time step

        reaction = rand_22 * total_rate;
        reaction_counter = rate[0];

        reaction_index = 0;
        while (reaction_counter < reaction)       // condition for which reaction to be excuted
        {
            reaction_index++;
            reaction_counter += rate[reaction_index];
        }

        if(reaction_index==0){index0++;}
        if(reaction_index==1){index1++;}
        if(reaction_index==2){index2++;}
        if(reaction_index==3){index3++;}

        over_den = 0.0;
        over_num = 0.0;

        time += tau;       // reaction time, sec
        D_time += tau;
//        d_time += tau;
        d_time_0 += tau;
        d_time_1 += tau;

        pre_add_spe3 += species[3]*tau;
        pre_d_time_spe3 += tau;

        if (species[2]!=0){
            pre_add_spe2 += species[2]*tau;
            pre_d_time_spe2 += tau;
        }

//calculate the SF at beginning with 1
        if (d_time_1 > 0.0001 && flag0 == 0 && pre_add_spe2!=0){

//ODE
            ode(stiff_sys, system_jacobi,r,Rad, d_time_1 );
            radical_ode3 = r[3];
            radical_ode2 = r[2];
            scale_fac =1/ (Na * volume) / radical_ode3;
//            scale_fac_s2 = (pre_add_spe2/pre_d_time_spe2/(Na * volume))/radical_ode2;

//            pre_add_spe2 = 0;
//            pre_d_time_spe2 = 0;

            flag0 = 1;
            sf<<time<<"  "<<scale_fac<<"  "<<scale_fac_s2<<"    "<<s3<<" "<< pre_add_spe3/pre_d_time_spe3<<endl;
            cout<<scale_fac<<" "<< pre_add_spe3/pre_d_time_spe3<<endl;
//            d_time_1 = 0.0;

//            cout<<"000         000"<<endl;
        }


// calculate the SF within 100s (average)
        if (d_time_1>1 && pre_d_time_spe3!=0.0 && fl_aver==0){

                aver_spe3 = pre_add_spe3/pre_d_time_spe3;
//              cout<<aver_spe3<<endl;

                aver<<aver_spe3<<endl;

                ode(stiff_sys, system_jacobi,r,Rad,d_time_1);
                radical_ode3=r[3];
                radical_ode2 = r[2];
//                if ((aver_spe3/(Na * volume *radical_ode3)) > 1){
                    scale_fac = (aver_spe3/(Na * volume))/radical_ode3;
 //                   scale_fac_s2 = (pre_add_spe2/pre_d_time_spe2/(Na * volume))/radical_ode2;
                    fl_aver=1;
//                    cout<<scale_fac<<" "<<aver_spe3<<endl;
//                 }



                sf<<" "<<scale_fac<<"    "<<scale_fac_s2<<"    "<<aver_spe3<<"   "<< aver_spe3/(Na * volume)<<"   "<<radical_ode3<<endl;

//                cout<<"111         000"<<endl;
        }

//setting for re-update with 1
            if ( d_time_1 >1){

            d_time_1 = 0.0;
            fl_aver=0;
//            pre_add_spe3 = 0;
//            pre_d_time_spe3 = 0;
//            pre_add_spe2 = 0;
//            pre_d_time_spe2 = 0;
//            cout<<"222         000"<<endl;
        }

            if (D_time >= 10) {

            D_time = 0;
            speciesout << time;
            species_con << time;
            index << time;
            ode3<<time;


            speciesout<<" " << species[0]<<" "<< species[1]<<" "<< species[2]<<" "<< species[3]<<" "<< species[4]<<" " << species[5]<<" "<< endl;
            species_con<<" " << species[0]/ (Na * volume)<<" "<< species[1]/ (Na * volume)<<" "<< species[2]/ (Na * volume)<<" "<< species[3]/ (Na * volume)<<" "<< 1-species[4]/ (4.67*Na * volume)<<" " << species[5]/ (Na * volume)<<" "<< endl;
            index << " " << index0 << " " << index1 << " " << index2 << " " << index3 <<endl;
            index << " " << rate[0] << " " << rate[1]  << " " << rate[2]  << " " << rate[3]  <<endl;

            ode3<<"   "<<scale_fac<<"    "<<pre_add_spe3/pre_d_time_spe3/(Na * volume)<<"    "<<radical_ode3<<"    "<<s3 <<endl;

            rat <<time<<"  " <<rate[0]<<" " <<rate[1]<<" " <<rate[2]<<" " <<rate[3]<<endl;

            conversion[0] = num_monomer[0] / M_initial[0];
            M_sum += M_start[0] + M_present[0];
            M_diff += M_start[0] - M_present[0];
            over_num += num_monomer[0];
            over_den += M_initial[0];

            M_sum += M_start[0] + M_present[0];
            M_diff += M_start[0] - M_present[0];
            over_num += num_monomer[0];
            over_den += M_initial[0];

            over_x = over_num / over_den;
            conv << time << "\t" << conversion[0] <<  "\t" << over_x << endl;

            f_inst[0] = (M_present[0] + M_start[0]) / (M_sum);
            F_inst[0] = (M_start[0] - M_present[0]) / (M_diff);
            //          f_inst[1] = (M_present[1] + M_start[1]) / (M_sum);
            //          F_inst[1] = (M_start[1] - M_present[1]) / (M_diff);
            M_start[0] = M_present[0];
            //          M_start[1] = M_present[1];
            num_M[0] = 0;
            //          num_M[1] = 0;
            num_M_sum = 0.0;
            species_sum = 0.0;
            num_M_sum = 0.0;
            species_sum = 0.0;

                if (fl_mwd == 1) {
                    molecular_weight(mol_wt,
                                     Mn,
                                     Mw,
                                     num_M,
                                     num_Dr,
                                     num_Pr,
                                     num_P,
                                     Dr,
                                     Pr,
                                     N,
                                     scale_fac,
                                     poly_test.chain1); // function for calculating molecular weight
                    pdi = Mw / Mn;
                    species_sum += species[4];
                    num_M_sum += num_M[0];


                    f[0] = species[4] / species_sum;   //avg compostion of monomers in the feed
                    F[0] = num_M[0] / (num_M_sum);       //avg compostion of monomers in the dead end polymer chains

                    Dpn = num_M_sum / (poly_test.chain1.size()+ Dr.size()+ Pr.size()/scale_fac);    // degree of polymerization
                    molwt << time << "\t" << Mn << "\t" << "\t" << Mw << "\t" << "\t" << pdi << "\t" << "\t" << Dpn << endl;
                    fract << time << "\t" << f[0] << "\t" << F[0] <<  "\t" << f_inst[0]
                          << "\t" << F_inst[0]  << endl;

                }

        }           // Data analuysis loop ends.............//


        efficient_explicit_sequence_record_number(reaction_index,
                                                  num_monomer,
                                                  species,
                                                  num_chain,
                                                  Dr,
                                                  Pr,
                                                  allchains,
                                                  poly_test
        );  //function for updating species info when event occurs


        all_scale_fac += scale_fac*tau;

        // update rates

        rate[0] = c[0] * species[0]* species[1]/ (Na * volume);
        rate[1] =1 / scale_fac /scale_fac_s2  *  c[1] / (Na * volume) * species[2] * species[3];
        rate[2] =1 / scale_fac *  c[2] / (Na * volume) * species[3] * species[4];
        rate[3] = 1 / scale_fac / scale_fac * c[3] * 2 / (Na * volume) * species[3] * (species[3] - 1) / 2;
//        rate[3] = 1  * c[3] * 2 / (Na * volume) * species[3] * (species[3]-1) / 2;

        M_present[0] = species[4];
        //      M_present[1] = species[3];
        M_sum = 0.0;
        M_diff = 0.0;

        // Data analuysis (for conversion and composition etc) loop begins...............//
    }



    /*................................ KMC Loop ends here..........................*/
    if (fl_mwd == 0) {
        molecular_weight(mol_wt,
                         Mn,
                         Mw,
                         num_M,
                         num_Dr,
                         num_Pr,
                         num_P,
                         Dr,
                         Pr,
                         N,
                         scale_fac,
                         poly_test.chain1); // function for calculating molecular weight
        pdi = Mw / Mn;
        species_sum += species[4];
        num_M_sum += num_M[0];

        for(int i=0;i<319;i++){
            wf<<i+1<<" "<<(i+1)*N[i]/num_M_sum <<endl;
        }


        f[0] = species[4] / species_sum;   //avg compostion of monomers in the feed
        F[0] = num_M[0] / (num_M_sum);       //avg compostion of monomers in the dead end polymer chains

        Dpn = num_M_sum / (poly_test.chain1.size()+ Dr.size()+ Pr.size()/scale_fac);    // degree of polymerization
        molwt << time << "\t" << Mn << "\t" << "\t" << Mw << "\t" << "\t" << pdi << "\t" << "\t" << Dpn << endl;
        fract << time << "\t" << f[0] << "\t" << F[0] <<  "\t" << f_inst[0]
              << "\t" << F_inst[0]  << endl;

    }


    end=clock();
    ode3<<"runnning time = "<<(double)(end-start)/CLOCKS_PER_SEC<<"s"<<endl;
    return 0;
}
