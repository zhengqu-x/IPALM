#ifndef QCQP_B_H
#define QCQP_B_H



#include "Matrix.h"
#include "ALM_L_Katyusha.h"
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>      /* printf */
#include <time.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <math.h>

/*This class solves problem of the form
min 1/2 x^\top Q_0 x+b_0^\top x
s.t. 1/2 x^\top Q_i x+b_i^\top x -1\leq 0,\enspace \forall i=1,\dots,m
-b\leq x\leq b

using IPALM with loopless Katyusha as inner solver.
*/



template<typename L, typename D>
class QCQP_b: public ALM_L_Katyusha<L, D>
{
private:

  vector<D> Lambda; // maximal eigenvalues of Q0, Q1,...Qm.

  vector<D> norm_c; //||[b0;b1;...;bm] ||_2

  vector<D> norm_1_c; //||[b0;b1;...;bm] ||_1

  D b;  // absolute bound on each variable

  D nb_q_c;  //number of quadratic constraints

  D nb_v; //number of variables

  Matrix<L,D> data_Q; // Q=[Q0;Q1;...;Qm]; b=[b0;b1;...;bm].



protected:

public:

  QCQP_b(const char* Matrix_file,D val_b)
  :ALM_L_Katyusha<L,D>(), data_Q(Matrix_file)
  {
    b= val_b;
    nb_v=data_Q.nfeatures;
    nb_q_c=data_Q.nsamples/data_Q.nfeatures- 1;
    this->mu_g=0;
  }

  inline D get_n(){
    return data_Q.nfeatures;
  }

  inline D get_m(){
    return data_Q.nsamples/data_Q.nfeatures- 1;
  }

 void set_dimension(){
    this->d= nb_v;
    this->m= nb_q_c;
  }

  void set_is_feasibility_constraint(){
    for(L i=0;i<nb_q_c;i++)
    this->is_feasiblity_constraint[i]=1;
  }


  inline D value_of_P(vector<D> & x){
    D res=0;
    for(L i=0;i<this->d;i++)
    {
      if (fabs(x[i])> b){
        cout<< "error: the variable x_"<< i<<"="<<x[i]<<". It should be in [-"<<b<<", "<<b<<" ]"<< endl;
        system("pause");
        return std::numeric_limits<double>::max();
      }
    }
    return res;
  }


  inline void prox_of_P(D alpha,vector<D> & x, vector<D> &  y){
    for(L i=0;i<this->d;i++)
    {
      D x2=x[i];
      if (x2>= b){
        y[i]=b;
      }
      else if(x2<=-b){
        y[i]=-b;
      }
      else{
        y[i]=x2;
      }
    }
  }


  D prox_of_h_star_j(D x,D beta, L j){
    if(x>=0) return x;
    else return 0;
  }

  D value_of_h_star_j(D x,L j){
    if(x>=0) return 0;
    else {
      cout<< "Error! x="<<x<<" should be nonnegative."<<endl;
      system("pause");
      return std::numeric_limits<double>::max();
    }
  }

  void compute_gradient_f0(vector<D> & x, vector<D>& g){
    for(L j=0;j<nb_v;j++){
      g[j]=0;
      for(L k = data_Q.ptr[j]; k < data_Q.ptr[j + 1];k++){
        L kj=data_Q.row_idx[k];
        g[j]+=x[kj]*data_Q.A[k];
      }
      g[j]+=data_Q.b[j];
    }
  }

  D compute_gradient_and_value_f_i(vector<D> & x ,vector<D> & g, L i){
    D fi=0;
    for(L j=0;j<nb_v;j++){
      g[j]=0;
      L j2=i*nb_v+j;
      for(L k = data_Q.ptr[j2]; k < data_Q.ptr[j2 + 1];k++){
        L kj=data_Q.row_idx[k];
        g[j]+=x[kj]*data_Q.A[k];
      }
      D bj=data_Q.b[j2];
      fi+=x[j]*g[j]/2+x[j]*bj;
      g[j]+=bj;
    }
    return fi-1;
  }




  D value_of_f0(vector<D> & x){
    return my_value_of_f_i( x, 0);
  }

  D value_of_f_i(vector<D> & x, L i){
    return my_value_of_f_i( x, i)-1;
  }

  D my_value_of_f_i(vector<D> & x, L i){
    D fi=0;
    for(L j=0;j<nb_v;j++){
      D Qjx=0;
      L j2=i*nb_v+j;
      for(L k = data_Q.ptr[j2]; k < data_Q.ptr[j2 + 1];k++){
        L kj=data_Q.row_idx[k];
        Qjx+=x[kj]*data_Q.A[k];
      }
      fi+=x[j]*Qjx/2+x[j]*data_Q.b[j2];
    }
    return fi;
  }


  D value_of_h_j(D x, L i){
    cout<< "Error! The "<<i<<"th function should be a feasibility constraint"<< endl;
    system("pause");
    return std::numeric_limits<double>::max();
  }

  D distance_to_domain_of_h_j(D x, L i){
    if(x<=0)
    return 0;
    else
    return x;
  }


  D compute_Qi_x(vector<D> & x ,vector<D> & g, L i){
    D normg=0;
    for(L j=0;j<nb_v;j++){
      g[j]=0;
      L j2=i*nb_v+j;
      for(L k = data_Q.ptr[j2]; k < data_Q.ptr[j2 + 1];k++){
        L kj=data_Q.row_idx[k];
        g[j]+=x[kj]*data_Q.A[k];
      }
      normg+=g[j]*g[j];
    }
    return sqrt(normg);
  }


  void compute_lambda(){
    D maxv=0;
    D minv=std::numeric_limits<double>::max();
    Lambda.resize(nb_q_c+1);
    norm_c.resize(nb_q_c+1,0);
    norm_1_c.resize(nb_q_c+1,0);
    std::vector<D> xk(nb_v,1);
    std::vector<D> yk(nb_v);
    D normk;
    D tmp;
    for (L i= 0; i<nb_q_c+1 ; i++){
      for(L kk=0;kk<20;kk++){
        normk= compute_Qi_x(xk , yk ,i);
        for (L j=0;j<nb_v;j++)
        xk[j]=yk[j]/normk;
      }
      D res=0;
      normk= compute_Qi_x(xk, yk ,i);
      for (L j=0;j<nb_v;j++)
      res+=xk[j]*yk[j];
      Lambda[i]= res;
      if (res> maxv){
        maxv= res;
      }
      if (res< minv){
        minv= res;
      }
      D res2= 0;
      D res3= 0;
      for (L j=0; j< nb_v; j++){
        res2+= data_Q.b[i*nb_v+j]*data_Q.b[i*nb_v+j];
        res3+= fabs(data_Q.b[i*nb_v+j]);
      }
      norm_c[i]= sqrt(res2);
      norm_1_c[i]=res3;
    }
    cout<< "max Lambda= "<< maxv<< " min Lambda= "<< minv<< endl;
  }


  inline void set_Li_Lf(){
    this->Li.clear();
    this->Li.resize(this->nsamples,0);
    this->Li[0]=this->nsamples*Lambda[0];
    this->Lf=Lambda[0];
    D maxLi=this->Li[0];
    D minLi=this->Li[0];
    for (L i=0; i< nb_q_c; i++){
      D M_p_j= nb_v*Lambda[i+ 1]*b*b+ b*norm_1_c[i+1]+ 1;
      D M_dp_j= sqrt(nb_v)*Lambda[i+ 1]*b+ norm_c[i+1];
      D L_dp_j= Lambda[i+1];
      D d_s_j= M_p_j+ this->beta_s*this->lambda_s[i];
      D tmp=L_dp_j*d_s_j/this->beta_s+ M_dp_j*M_dp_j/this->beta_s;
      this->Li[i+1]= this->nsamples*tmp;
      this->Lf+=tmp;
      if (this->Li[i+1]> maxLi){
        maxLi= this->Li[i+ 1];
      }
      if (this->Li[i+ 1]< minLi){
        minLi= this->Li[i+ 1];
      }
    }
    this->sumLi=this->Lf*this->nsamples;
    cout<<"  max of Li: "<<maxLi<<"  min of Li: "<<minLi<<" sumLi="<<this->sumLi<<" Lf= "<<this->Lf<<endl;
  }




  void ALM_QCQP_solver(D beta_0, D epsilon_0,  D eta, D rho,vector<D> & x0,vector<D> & y0,L val_tau, L max_nb_outer, L p_N_1, L p_N_2,string filename1, string filename2, D time){
    compute_lambda();
    this->ALM_solve_with_L_Katyusha(beta_0, epsilon_0,  eta, rho,x0,y0, val_tau, max_nb_outer,  p_N_1,  p_N_2,  filename1,  filename2, time);
  }



};

#endif
