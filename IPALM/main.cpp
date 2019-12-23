//using namespace std;
#include <iostream>
#include <string>
#include "Basis_pursuit3.h"
#include "Basis_pursuit4.h"
#include "Basis_pursuit6.h"
#include "L_1_Lasso3.h"
#include "L_1_Lasso4.h"
#include "L_1_Lasso5.h"
#include "Fused_Lasso.h"
#include "Fused_lasso3.h"
#include "Fused_lasso4.h"
#include "Fused_lasso5.h"
#include "SMSVM.h"
#include "SMSVM2.h"
#include "SMSVM3.h" 
#include "SMSVM4.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>






    int main(int argc,char * argv[])
    {
    	if(argc<8)  std::cerr << "try:  ./main 4 b w4a 1 1 100000 100 1 1 50 " << std::endl;
        int i= atoi(argv[1]);
        char j= argv[2][0];
        double lambda1= atof(argv[4]);
        double lambda2= 0.0;
        double lambda3= atof(argv[5]);
        long n;
        long m;
        long p_N_1= 1;
        long p_N_2= 100;
		std::vector<double> x0;
        //for (int j = 0; j < n; j++)  {readx0_ >>x0[j]; }
        std::vector<double> y0;
        double beta_0= 1;
        double epsilon_0= 0;
        double eta;
        double rho;
		long val_tau=1;
        stringstream  tau_convert;
        string tau_str;
     	string filename, filename1, filename2, filenameMatrix, filenameMatrix_M;
     	filename= argv[3];
     	filenameMatrix="../datas/matrix_"+filename;
     	long max_nb_outer=atoi(argv[6]);
     	double time= atof(argv[7]);
     	if(argc== 9){
     		beta_0= atof(argv[8]);
		 }
		 else if(argc== 10){
		 	beta_0= atof(argv[8]);
		 	p_N_2= atoi(argv[9]);
		 }
		 else{
		 	beta_0= atof(argv[8]);
		 	epsilon_0= atof(argv[9]);
		 	p_N_2= atoi(argv[10]);
		 }
        switch(i){
        	case 1: {
            switch(j){
            	case 'a':{
			    Basis_pursuit4<long,double> bp4(filenameMatrix.c_str(),lambda1,lambda2);
                n=bp4.get_n();
                m=bp4.get_m(); 
                x0.resize(n,0);
                y0.resize(m,0);
                val_tau= 1;
                tau_convert<<val_tau;
				tau_str=tau_convert.str();
			    filename1="BP_outer_"+filename+"tau_"+tau_str;
                filename2="BP_inner_"+filename+"tau_"+tau_str;
                rho= 0.9;
                eta= 0.95;
                bp4.ALM3_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time);
                break;
                }
                case 'b': {
				cout<<"Numbber of features is much larger than samples. Katyusha is inefficent."<< endl;
                break;
                }
                case 'c':{
			    Basis_pursuit6<long,double> bp6(filenameMatrix.c_str(),lambda1,lambda2);
                n=bp6.get_n();
                m=bp6.get_m(); 
                x0.resize(n,0);
                y0.resize(m,0);
                val_tau= 1;
                tau_convert<<val_tau;
				tau_str=tau_convert.str();
			    filename1="BP_outer_"+filename+"tau_"+tau_str;
                bp6.SMART_CD_solver(beta_0,x0,y0,val_tau,max_nb_outer,p_N_2,filename1,time);
                break;
                }
                case 'd':{
			    Basis_pursuit3<long,double> bp3(filenameMatrix.c_str(),lambda1,lambda2);
                n=bp3.get_n();
                m=bp3.get_m(); 
                x0.resize(n,0);
                y0.resize(m,0);
                val_tau= n;
                tau_convert<<val_tau;
				tau_str=tau_convert.str();
			    filename1="BP_outer_"+filename+"tau_"+tau_str;
			    filename2="BP_inner_"+filename+"tau_"+tau_str;
                rho= 0.9;
                eta= 0.95;
                bp3.ALM2_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time);
                break;
                }
				default: {
                cout<< "no inner solver input"<< endl;	
                break;
                }
            }
            break;
            }
            case 2: {
            switch(j){
            	case 'a':{
				L_1_Lasso4<long,double> l1l4(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=l1l4.get_n();
            	m=l1l4.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= 1;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="LAD_outer_"+filename+"tau_"+tau_str;
            	filename2="LAD_inner_"+filename+"tau_"+tau_str;
            	rho= 0.95;
                eta= 0.95;
            	l1l4.DLRCSGR3_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time);
            	break;
                }
            	case 'b': {
				cout<<"Numbber of features is much larger than samples. Katyusha is inefficent."<< endl;
				break;
			    }
			    case 'c':{
				L_1_Lasso5<long,double> l1l5(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=l1l5.get_n();
            	m=l1l5.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= 1;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="LAD_outer_"+filename+"tau_"+tau_str;
            	l1l5.SMART_CD_solver(beta_0,x0,y0,val_tau,max_nb_outer,p_N_2, filename1, time);
            	break;
                }
                case 'd':{
				L_1_Lasso3<long,double> l1l3(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=l1l3.get_n();
            	m=l1l3.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= n;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="LAD_outer_"+filename+"tau_"+tau_str;
            	filename2="LAD_inner_"+filename+"tau_"+tau_str;
                eta= 1.0/1.2;
            	l1l3.DLRCSGR2_solver(beta_0,epsilon_0,1.0/eta,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time);
            	break;
                }
				default: {
					cout<< "no inner solver input"<< endl;
					break;
				}	
			}
			break;
		    }
            case 3: {
            switch(j){
            	case 'a':{
                Fused_lasso4<long,double> fl4(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=fl4.get_n();
            	m=fl4.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= 1;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="FL_outer_"+filename+"tau_"+tau_str;
            	filename2="FL_inner_"+filename+"tau_"+tau_str;
            	rho= 0.95;
                eta= 0.95;
            	fl4.DLRCSGR3_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time); 
            	break;
                }
            	case 'b':  {
				Fused_lasso<long,double> fl(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=fl.get_n();
            	m=fl.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= floor(sqrt(m));;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="FL_outer_"+filename+"tau_"+tau_str;
            	filename2="FL_inner_"+filename+"tau_"+tau_str;
            	rho= 0.9;
                eta= 0.95;
            	fl.Katyusha_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time); 
            	break;
                }
                case 'c':{
                Fused_lasso5<long,double> fl5(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=fl5.get_n();
            	m=fl5.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= 1;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="FL_outer_"+filename+"tau_"+tau_str;
            	fl5.SMART_CD_solver(beta_0,x0,y0,val_tau,max_nb_outer,p_N_2,filename1,time); 
            	break;
                }
                case 'd':{
                Fused_lasso3<long,double> fl3(filenameMatrix.c_str(),lambda1,lambda2,lambda3);
            	n=fl3.get_n();
            	m=fl3.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= n;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="FL_outer_"+filename+"tau_"+tau_str;
            	filename2="FL_inner_"+filename+"tau_"+tau_str;
                eta= 1.0/1.2;
            	fl3.DLRCSGR2_solver(beta_0,epsilon_0,1.0/eta,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time); 
            	break;
                }
				default: {
					cout<< "no inner solver input"<< endl;
					break;
				}	
			}
			break;
		    }
			case 4: {
            switch(j){
            	case 'a': {
			    SMSVM3<long,double> svm3(filenameMatrix.c_str(),lambda1,lambda2);
            	n=svm3.get_n();
            	m=svm3.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= 1;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="SVM_outer_"+filename+"tau_"+tau_str;
            	filename2="SVM_inner_"+filename+"tau_"+tau_str;
            	rho= 0.9;
                eta= 0.95;
            	svm3.DLRCSGR3_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time); 
            	break;
                }
            	case 'b': {
				SMSVM<long,double> svm(filenameMatrix.c_str(),lambda1,lambda2);
            	n=svm.get_n();
            	m=svm.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= floor(sqrt(m));;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="SVM_outer_"+filename+"tau_"+tau_str;
            	filename2="SVM_inner_"+filename+"tau_"+tau_str;
            	rho= 0.9;
                eta= 0.95;
            	svm.Katyusha_solver(beta_0,epsilon_0,eta,rho,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2, time); 
            	break;
                }
                case 'c': {
			    SMSVM4<long,double> svm4(filenameMatrix.c_str(),lambda1,lambda2);
            	n=svm4.get_n();
            	m=svm4.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= 1;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="SVM_outer_"+filename+"tau_"+tau_str;
            	svm4.SMART_CD_solver(beta_0,x0,y0,val_tau,max_nb_outer,p_N_2,filename1,time); 
            	break;
                }
                case 'd': {
			    SMSVM2<long,double> svm2(filenameMatrix.c_str(),lambda1,lambda2);
            	n=svm2.get_n();
            	m=svm2.get_m(); 
            	x0.resize(n,0);
                y0.resize(m,0);
            	val_tau= n;
            	tau_convert<<val_tau;
				tau_str=tau_convert.str();
				filename1="SVM_outer_"+filename+"tau_"+tau_str;
            	filename2="SVM_inner_"+filename+"tau_"+tau_str;
                eta= 0.95;
            	svm2.DLRCSGR2_solver(beta_0,epsilon_0,1.0/eta,x0,y0,val_tau,max_nb_outer,p_N_1,p_N_2,filename1,filename2,time); 
            	break;
                }
				default: {
					cout<< "no inner solver input"<< endl;
					break;
				}	
			}
            break;
            }
            default:{
            	cout<< "No such type of problem"<< endl;
				break;
			}
		}
		 

        
 
    }
