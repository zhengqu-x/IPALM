#ifndef ALM_L_KATYUSHA_H
#define ALM_L_KATYUSHA_H

#include "L_Katyusha.h"


#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>      /* printf */
#include <time.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <sstream>
//#include "cmd_line.h"


//This class implements the method IPALM_L_Katyusha

/*
The optimization problem to solve is:

min f_0(x)+\sum_{i=1}^m h_i(f_i(x))+P(x)

x is a d-dimensional vector.

// Each subproblem solves problem of the form
f_0(x)+ \sum_{i=1}^m h_{\beta_s}^i(f_i(x);\lambda_s^i) +g(x)+ \beta_s/2\|x-x_s\|^2
 by L_Katyuhsa
*/

template<typename L, typename D>
class ALM_L_Katyusha: public L_Katyusha<L, D>
{
private:



protected:


  std::vector<D> x_s;

  std::vector<D> old_x_s;

  std::vector<D> lambda_s;

  std::vector<D> old_lambda_s;

  std::vector<D> x_tmp;

  std::vector<D> g_tmp;

  std::vector<D> g_tmp_x;
  std::vector<D> g_tmp_w;
  std::vector<D> valuefi;




  std::vector<D> is_feasiblity_constraint;

  D beta_s;

  D epsilon_s;

  D m_s;

  D m_0;

  L m;

  L d;

  D eta;

  D rho;

  D tau_s;

  std::vector<D> L_phi;

  D function_value;

  D infeasibility;

  L print_every_N_ALM;

  D running_time_ALM;

  L nb_outer_iters;

  ofstream samp_ALM;


public:

  D mu_g;

  D lambda1;

  D lambda2;




  virtual inline void set_Li_Lf(){}

  virtual inline D value_of_P(vector<D> &){return D(NULL);}
  virtual inline void prox_of_P(D, vector<D> &, vector<D> &){} //prox_of_P(L,x,y) computes y=argmin{P(u)+L/2||u-x||^2}
  virtual inline D prox_of_h_star_j(D,D, L){return D(NULL);} //prox_of_h_star_j(x,beta,j) computes argmin{h*(y)+beta/2(y-x)^2}
  virtual inline D value_of_h_star_j(D, L){return D(NULL);}
  virtual inline void compute_gradient_f0(vector<D> &,vector<D> &){} //compute_gradient_f0(x,g) computes g=nabla f0(x)
  virtual inline D compute_gradient_and_value_f_i(vector<D> &,vector<D> &, L){return D(NULL);}//compute_gradient_and_value_f_i(i,x,g) computes g=nabla f_i(x) and return f_i(x)


  virtual inline D value_of_h_j(D, L){return D(NULL);}
  virtual inline D distance_to_domain_of_h_j(D,L){return D(NULL);}

  virtual inline D value_of_f0(vector<D> &){return D(NULL);} //f0(x)
  virtual inline D value_of_f_i(vector<D> &,L){return D(NULL);} // f_i(x)
  virtual inline void set_dimension(){} // set the value of m and d
  virtual inline void set_is_feasibility_constraint(){}

  ALM_L_Katyusha()
  : L_Katyusha<L,D>()
  {

  }



  D value_of_phi_i(L i){
    D res= 0;
    if (i==0)
    res=value_of_f0(this->x);
    else
    {
      D fi=value_of_f_i(this->x,i);
      i=i-1;
      D lambdasi=lambda_s[i];
      D tmp= prox_of_h_star_j(fi/beta_s+ lambdasi,beta_s,i);
      res=fi*tmp- value_of_h_star_j(tmp,i)- beta_s/2*(tmp- lambdasi)*(tmp- lambdasi);
    }
    return this->nsamples*res;
  }

  void compute_full_gradient(vector<D> & x, vector<D> & gx){
    compute_gradient_f0(x,gx);
    for(L i=0;i<m;i++){
          D fi=compute_gradient_and_value_f_i(x,g_tmp,i+1);
          D tmp= prox_of_h_star_j(fi/beta_s+ lambda_s[i],beta_s,i);
          for(L j=0;j<d;j++)
            gx[j]+=g_tmp[j]*tmp;
    }
  }

  void compute_batch_delta_gradient(){
    for (L j=0; j< this->batch_size; j++){
      L s= this->batch_i[j];
      if(s==0){
        compute_gradient_f0(this->x,g_tmp_x);
        compute_gradient_f0(this->w,g_tmp_w);
        add_up(1.,1.,this->theta_S[s]);
      }
      else{
        L i=s-1;
        D fi_x=compute_gradient_and_value_f_i(this->x,g_tmp_x,s);
        D tmp_x= prox_of_h_star_j(fi_x/beta_s+ lambda_s[i],beta_s,i);
        D fi_w=compute_gradient_and_value_f_i(this->w,g_tmp_w,s);
        D tmp_w= prox_of_h_star_j(fi_w/beta_s+ lambda_s[i],beta_s,i);
        add_up(tmp_x,tmp_w,this->theta_S[s]);
      }
    }
  }

  void  add_up(D fx, D fw, D c){
      for(L i=0;i<d;i++)
        this->batch_delta_gradient[i]+=(fx*g_tmp_x[i]-fw*g_tmp_w[i])*c;
    }



  inline D value_of_g(){
    D res= value_of_P(this->x);
    for(L i=0;i<d;i++)
      res+= 0.5*beta_s*(this->x[i]- x_s[i])*(this->x[i]- x_s[i]);
    return res;
  }

  void prox_of_g(D stepsz, vector<D> & x, vector<D> & g, vector<D> & nextx){
    D tmp0=beta_s+stepsz;
    for(L i=0;i<d;i++)
     x_tmp[i]=(beta_s*x_s[i]+stepsz*x[i]-g[i])/tmp0;
    prox_of_P(tmp0, x_tmp, nextx);
  }



  void compute_x(){
    for (L i=0; i< d; i++){
      x_s[i]= this->x[i];
    }
    for (L i= 0; i< m; i++){
      valuefi[i]=value_of_f_i(x_s,i+1);
    }
  }

  void update_lambda(){
    for(L j=0;j<m;j++)
    {
      D fi= valuefi[j];
      lambda_s[j]= prox_of_h_star_j(fi/beta_s+ lambda_s[j],beta_s,j);
    }
  }



  void compute_function_value(){
    function_value=value_of_f0(x_s);
    D res= 0;
    for (L i= 0; i< m; i++){
      D fi=value_of_f_i(x_s,i+1);
      if(is_feasiblity_constraint[i]==1){
        D tmp=distance_to_domain_of_h_j(fi,i);
        res+=tmp*tmp;
      }
      else
       function_value+=value_of_h_j(fi,i);
    }
    function_value+=value_of_P(x_s);
    infeasibility=sqrt(res);
  }





inline void compute_m0(D beta0, D val_eta, D val_rho, L val_tau){
  cout<<"beta0="<<beta0<<"val_eta="<<val_eta<<"; val_rho="<<val_rho<<"; val_tau="<<val_tau<<endl;
  m_s= this->nsamples/val_tau*1e+6;

  //D tmp7=1-val_tau/this->n*sqrt(beta0/(max_Lf_s+max_M_s/beta0+beta0));
  //cout<<"tmp7: "<<tmp7<<"; "<<val_tau/(this->n+0.)*sqrt(beta0/(max_Lf_s+max_M_s/beta0+beta0))<<endl;
  //m_0=(2*log(val_rho)+log(val_eta)+log(2))/log(tmp7);

  cout<<" m_s="<<m_s<<endl;
}



void Initialize(D beta_0, D epsilon_0,  D val_eta, D val_rho,L val_tau, vector<D> & x0,vector<D> & y0){
  cout<<"start initializing ALM"<<endl;
  set_dimension();
  this->nsamples=m+1;
  this->nfeatures=d;
  this->tau=val_tau;


  beta_s=beta_0;
  tau_s= 1;
  epsilon_s=epsilon_0;
  compute_m0(beta_s,val_eta,val_rho,val_tau);

  eta=val_eta;
  rho=val_rho;

  x_s.resize(d,0);
  old_x_s.resize(d,0);
  lambda_s.resize(m,0);
  old_lambda_s.resize(m,0);
  for(L i=0;i<d;i++){
    x_s[i]=x0[i];
    old_x_s[i]= x0[i];
  }
  for(L j=0;j<m;j++){
    lambda_s[j]=y0[j];
    old_lambda_s[j]= y0[j];
  }
  g_tmp.clear();
  g_tmp.resize(d,0);
  g_tmp_x.clear();
  g_tmp_x.resize(d,0);
  g_tmp_w.clear();
  g_tmp_w.resize(d,0);
  valuefi.clear();
  valuefi.resize(m,0);
  x_tmp.clear();
  x_tmp.resize(d,0);
  is_feasiblity_constraint.resize(m,0);

  set_is_feasibility_constraint();
  cout<<"Initialization ALM finished!"<<endl;

}

void reset_everything(){
  epsilon_s*=rho;
  beta_s*=eta;
}

void update_m_s(L val_tau){
  D tmp= 0;
  D tmpx= 0;
  D tmpy= 0;
  D tmpy2= 0;
  for(L j=0;j<m;j++)
  {
    D tmp_p= prox_of_h_star_j(valuefi[j]/beta_s/eta+ lambda_s[j], beta_s/eta, j);
    tmpy+=(tmp_p- lambda_s[j])*(tmp_p- lambda_s[j]);
    tmp+= (lambda_s[j]- old_lambda_s[j])*(lambda_s[j]- old_lambda_s[j]);
    tmpy2+= (beta_s*eta*lambda_s[j]- beta_s*old_lambda_s[j])*(beta_s*eta*lambda_s[j]- beta_s*old_lambda_s[j]);
  }
  for (L i= 0; i< d;i++){
    tmpx= (x_s[i]- old_x_s[i])*(x_s[i]- old_x_s[i]);
  }
  D tmp4= 2*epsilon_s+ beta_s*tmp;
  D tmp5= (1- eta)*beta_s/2*tmpy;
  D tmp6= 0;
  D tmp7= sqrt(tmp)*sqrt(tmpy2);
  D tmp1= log(tmp4+ tmp5+ tmp6+ tmp7+ 2/(2*eta- 1)*tau_s*beta_s/2*tmpx);
  D tmp2= log(epsilon_s)+ log(rho)- log(2);
  D tmp3= sqrt(beta_s/this->sumLi);
  m_s= ceil((tmp1- tmp2)/log(2)*4*max(1.0+this->nsamples,1/tmp3)/val_tau);
  cout<<"here here="<<tmp4<<" tmp5="<<tmp5<<" tmp6="<<tmp6<<" tmp7="<<tmp7<<" tmp1="<<tmp1<<" tmp2="<<tmp2<<" tmp3="<<tmp3<<endl;
}

inline void compute_and_record_res(){
  if(nb_outer_iters%print_every_N_ALM==0){
    compute_function_value();
    cout<<setprecision(9)<<"Iteration: "<<nb_outer_iters<<"; time="<<running_time_ALM<< "; function value="<<function_value<<"; infeasibility= "<< infeasibility<<endl;
    samp_ALM<<setprecision(9)<<nb_outer_iters<<" "<<running_time_ALM<<" "<< function_value<<" "<< infeasibility<<endl;
  }
}

void ALM_solve_with_L_Katyusha(D beta_0, D epsilon_0,  D eta, D rho,vector<D> & x0,vector<D> & y0, L val_tau, L max_nb_outer, L p_N_1, L p_N_2, string filename1, string filename2, D time){
  Initialize(beta_0, epsilon_0, eta, rho,val_tau, x0, y0);
  nb_outer_iters=0;
  filename1= "results/ALM_"+filename1;
  samp_ALM.open(filename1.c_str());
  running_time_ALM=0;
  print_every_N_ALM=p_N_1;
  this->set_print_every_N(p_N_2);
  compute_and_record_res();

 while(nb_outer_iters<max_nb_outer){
  D start= std::clock();
  cout<<"m_s= "<< ceil(m_s/this->nsamples*val_tau)<<"; beta_s="<<beta_s<<"; epsilon_s="<<epsilon_s<<endl;
  this->loopless2(x_s, filename2, mu_g+ tau_s*beta_s,  ceil(m_s/this->nsamples*val_tau), epsilon_s, val_tau, 1, 1, 0, 1);
  for(L i=0;i<d;i++){
    old_x_s[i]=x_s[i];
  }
  for(L i=0;i<m;i++){
    old_lambda_s[i]=lambda_s[i];
  }
  compute_x();
  update_lambda();
  update_m_s(val_tau);
  nb_outer_iters++;
  running_time_ALM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  compute_and_record_res();
  start = std::clock();
  reset_everything();
  running_time_ALM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    if (running_time_ALM> time){
      break;
    }
  }
}

};

#endif /* MIN_SMOOTH_CONVEX_H */
