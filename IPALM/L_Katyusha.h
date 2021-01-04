#ifndef L_Katyusha_H
#define L_Katyusha_H


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
#include <math.h>


//This class implements the loopless variance reduced type methods with arbitrary sampling

/*
The optimization problem to solve is:

F(x):=\sum_{i=1}^n 1/n* phi_i(x)+ g(x)
*/

template<typename L, typename D>
class L_Katyusha
{


protected:

  // involved variables







  std::vector<D> proba_vector;

  std::vector<D>  tilde_proba_vector;


  std::vector<D> group_C;

  std::vector<D> index_group_C;

  std::vector<D> maxp_group_C;

  std::vector<D> isolated_I;

  std::vector<D> sump_group_C;
  std::vector<D> theta_S;   //theta_S in the paper




  std::vector<D> Li;  // the Lipchitz constant of phi_i




  std::vector<D> gk;     // the vector g^k in the paper

  std::vector<D> gradient_f_w;  // gradient of f at w
  std::vector<D> gradient_f_x;  // gradient of f at x
  std::vector<D> batch_delta_gradient; //  1/n(G(x^k)-G(w^k))theta_{S_k} I_{S_k} e in the paper

  L c;

  L noverc;





  // auxiliary variables

  L nb_of_iters_per_loop;

  L max_nb_loops;

  D max_p;


  L nsamples;
  L nfeatures;


  D running_time;

  L nb_loops;



  L print_every_N;

  vector<L> batch_i;
  vector<L> my_batch;

  L batch_size;

  L nb_groups;    //number of groups in group sampling



  string uniform;

  D Lf;

  D L2;

  D sumLi;


  D upper_bound;   // upper bound of F(x)-F^*

  D p; // the probability of changing w to x as in the paper
  L tau;  //number of threads on each node/computer

  D mu;

  D scaler;


  L current_nb_iters;

  L nb_iters;

  ofstream samp;

  D theta1;
  D theta2;
  D theta3;
  D eta;  //eta/L in the paper
  D oneovereta;

  D primal_value;

  std::vector<D> x;
  std::vector<D> w;
  std::vector<D> y;
  std::vector<D> z;
  std::vector<D> next_x;


public:


  gsl_rng * rng;
  virtual inline D value_of_phi_i(L) {return D(NULL);}
  virtual inline D value_of_g(){return D(NULL);}
  virtual inline void prox_of_g(D, vector<D> &, vector<D> &, vector<D> &){}   //prox_of_g(L,x, gr, y) computes y=argmin_u{L/2 ||u-(x-gr/L)||^2+g(u)}
  virtual inline void set_auxiliary_v(){}
  virtual inline void update_gradient(){}
  virtual inline void set_Li_Lf(){}
  virtual inline void compute_batch_delta_gradient(){}
  virtual inline void compute_full_gradient(vector<D> &, vector<D> &){}


  L_Katyusha()
  {

  }


  void update_x(){
    for(L i=0; i<nfeatures; i++){
      x[i]= theta1*z[i]+ theta2*w[i]+ theta3*y[i];
    }
  }

  void update_y(){
    for(L i=0; i<nfeatures; i++){
      y[i]=theta1*z[i]+ theta2*w[i]+ theta3*y[i];
    }
  }

    void update_gk(){
      batch_delta_gradient.clear();
      batch_delta_gradient.resize(nfeatures,0);
      compute_batch_delta_gradient();
      for(L i=0;i<nfeatures;i++){
        gk[i]=batch_delta_gradient[i]+gradient_f_w[i];
      }
    }

    void update_z(){
      prox_of_g(oneovereta,z,gk,z);
    }

    void update_w(){
      D yi=gsl_rng_uniform(rng);
      if(yi<=p){
        w=x;
        compute_full_gradient(w, gradient_f_w);
      }
    }



    void compute_upper_bound_of_optimality_gap(){
      compute_full_gradient(x, gradient_f_x);
      prox_of_g(Lf,x,gradient_f_x,next_x);
      D tmp=0;
      for(L i=0;i<nfeatures;i++){
        tmp+=(next_x[i]-x[i])*(next_x[i]-x[i]);
        x[i]=next_x[i];
      }
      upper_bound=tmp*Lf*Lf/mu;
   	}



    void set_rng()
    {
      gsl_rng_env_setup();
      const gsl_rng_type * T;
      T = gsl_rng_default;
      rng = gsl_rng_alloc(T);
      gsl_rng_set(rng,time(NULL));
      //gsl_rng_set(rng, 27432042);

    }



    void set_print_every_N(L i){print_every_N=i;}




    void set_optimal_probability()
    {
      tilde_proba_vector.clear();
      tilde_proba_vector.resize(nsamples,0);
      proba_vector.clear();
      proba_vector.resize(nsamples,0);
      D sum=0;
      for(L i=0; i<nsamples;i++)
      {
        tilde_proba_vector[i]=Li[i];
        sum+=Li[i];
      }
      max_p=0;
      for(L i=0; i<nsamples;i++)
      {
        tilde_proba_vector[i]=tilde_proba_vector[i]/sum;
        proba_vector[i]=1-pow(1-tilde_proba_vector[i],tau);
        if(max_p<tilde_proba_vector[i])
        {
          max_p=tilde_proba_vector[i];
        }
      }
      cout<<"sum="<<sum<<"; proba"<<tilde_proba_vector[0]<<endl;
    }


    void set_group_sampling_probability(){
      proba_vector.clear();
      proba_vector.resize(nsamples,0);
      std::vector<D> q(nsamples);
      D sumq=0;
      cout<<"start setting group sampling probablity"<<endl;
      for(L i=0; i<nsamples;i++)
      {
        q[i]=Li[i];
        sumq+=q[i];
      }
      D maxq=0;
      L nb=0;
      D tmp=0;
      for(L i=0; i<nsamples;i++)
      {
        q[i]=q[i]/sumq*tau;
        if(q[i]>maxq) maxq=q[i];
        if(q[i]>1) {
          nb++; //count the number of elements larger than 1
          tmp+=q[i]-1;
        }
      }
      cout<<"maxq="<<maxq<<endl;
      if(maxq<=1){
        for(L i=0; i<nsamples;i++)
        {
          proba_vector[i]=q[i];
        }
      }else{
        for(L i=0; i<nsamples;i++)
        {
          if(q[i]>1) q[i]=1;
          else {
            D deltaq=min(1-q[i],tmp);
            tmp=tmp-deltaq;
            q[i]+=deltaq;
          }
          proba_vector[i]=q[i];
        }
      }
    }



    void set_L2(L p_mod, L u){
      if(p_mod==0){ //batch sampling mode
        D tmp=0;
        D st;
        for(L i=0;i<nsamples;i++){
          st=Li[i]/tilde_proba_vector[i];
          if (st>tmp) tmp=st;
        }
        L2=tmp/nsamples/tau;
        cout<<"L2="<<L2<<"; "<<tilde_proba_vector[0]<<endl;
      }
      else{ //group sampling mode
        D tmp=0;
        D st;
        for(L i=0;i<nsamples;i++){
          st=Li[i]/proba_vector[i];
          if(isolated_I[i]==1) st=st-Li[i];
          if(st>tmp) tmp=st;
        }
        L2=tmp/nsamples;
      }
      cout<<"set_L2="<<L2<<endl;
    }

    void set_uniform_probability(L p_mod)
    {
      if(p_mod==0) {
        //batch sampling (stochastic process as defined in the paper)
        proba_vector.clear();
        tilde_proba_vector.clear();
        tilde_proba_vector.resize(nsamples,1./nsamples);
        proba_vector.resize(nsamples,1-pow(1-1.0/nsamples,tau));
        max_p=1.0/nsamples;
        cout<<"uniform="<<tilde_proba_vector[0]<<endl;
      }
      else{ //group sampling
        D pi=(tau*c+0.0)/nsamples;
        proba_vector.clear();
        proba_vector.resize(nsamples,pi);
      }
    }




    void sort_p(){
      vector<pair<D,L> >a;
      for (L i = 0 ;i < nsamples ; i++) {
        a.push_back(make_pair(proba_vector[i],i)); // k = value, i = original index
      }
      sort(a.begin(),a.end());
      group_C.clear();
      group_C.resize(nsamples);
      isolated_I.clear();
      isolated_I.resize(nsamples,0);
      D tmp=0;
      D maxpi=0;
      D previous_i=0;
      index_group_C.clear();
      index_group_C.push_back(0);
      maxp_group_C.clear();
      sump_group_C.clear();

      for (L i = 0 ;i < nsamples ; i++){
        group_C[i]=a[nsamples-1-i].second;
        D pi=a[nsamples-1-i].first;
        if(tmp+pi<=1){
          if(pi>maxpi) maxpi=pi;
          tmp+=pi;
        }
        else{
          index_group_C.push_back(i);
          maxp_group_C.push_back(maxpi);
          sump_group_C.push_back(tmp);
          tmp=pi;
          maxpi=pi;
          if(i-previous_i==1) isolated_I[group_C[previous_i]]=1;
          previous_i=i;
        }
      }
      index_group_C.push_back(nsamples);
      maxp_group_C.push_back(maxpi);
      sump_group_C.push_back(tmp);
      nb_groups=sump_group_C.size();
      if(previous_i==nsamples-1) isolated_I[group_C[previous_i]]=1;
      cout<<"size of group="<<nb_groups<<endl;
    }


    inline L sampling(L n)
    {
      L i=(floor)(gsl_rng_uniform(rng)*n);
      D y=gsl_rng_uniform(rng);
      while(y*max_p>tilde_proba_vector[i])
      {
        i=(floor)(gsl_rng_uniform(rng)*n);
        y=gsl_rng_uniform(rng);
      }
      return i;
    }


    void set_theta(){
      theta2=L2/2/max(Lf,L2);
      if(Lf<=L2/p){
        D tmp=mu/L2/p;
        cout<<"L2="<<L2<<endl;
        cout<<"tmp="<<tmp<<endl;
        if(tmp>=1){
          theta1=theta2;
        }
        else{
          theta1=sqrt(tmp)*theta2;
        }
      }else
      theta1=min(sqrt(mu/Lf),p/2);
      theta3=1-theta1-theta2;
      eta=1./(theta1*(Lf+2*max(L2,Lf)));
      oneovereta=theta1*(Lf+2*max(L2,Lf));
      cout<<"Lf="<<Lf<<"; L2="<<L2<<"; theta1="<<theta1<<"; theta2="<<theta2<<"; theta3="<<theta3<<endl;
      cout<<"eta="<<eta<<endl;
    }




    void set_p(D scal_p){
        p=scal_p*tau/(0.0+nsamples);
        cout<<"changing probablity: "<<p<<endl;
    }




    void initialize(vector<D> & x0, D val_mu, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
    {
      cout<<"start initializing loopless Katyusha"<<" u="<<u<<endl;
      c=nb_c;
      noverc=nsamples/c;
      if(nb_tau>noverc) perror("tau should be less than n over c");
      tau=nb_tau;
      nb_of_iters_per_loop=floor(nsamples/(c*(tau+0.0)));
      batch_i.clear();
      batch_i.resize(nsamples,0);


      gk.clear();
      gk.resize(nfeatures,0);
      upper_bound=std::numeric_limits<double>::max();

      /**setup parameters**/
      max_nb_loops=max_nb;
      mu=val_mu;
      cout<<"mu="<<mu<<endl;
      set_rng();

      set_Li_Lf();

      running_time=0;
      /**setup probability**/
      if(p_mod==0)   // batch sampling (stochastic process as defined in the paper)
      {
        if(u==0)   // uniform sampling
        {
          set_uniform_probability(p_mod);
          uniform="uniform";
        }
        else{
          set_optimal_probability();
          uniform="nonuniform";
        }
      }
      else{  //group sampling
        if(u==0){
          set_uniform_probability(p_mod);
          uniform="uniform";
        }else{
          set_group_sampling_probability();
          uniform="nonuniform";
        }
        sort_p();
        cout<<"sort p"<<endl;
      }


      set_L2(p_mod,u);
      scaler=scal_p;
      set_p(scal_p);
      set_thetaS(p_mod);
      cout<<"set thetaS"<<endl;


      x.clear();
      y.clear();
      z.clear();
      w.clear();
      x.resize(nfeatures);
      y.resize(nfeatures);
      z.resize(nfeatures);
      w.resize(nfeatures);
      next_x.clear();
      next_x.resize(nfeatures);

      gradient_f_w.clear();
      gradient_f_w.resize(nfeatures);
      gradient_f_x.clear();
      gradient_f_x.resize(nfeatures);

      for(L i=0;i<nfeatures;i++)
      {
        x[i]=x0[i];
        y[i]= x0[i];
        w[i]=x0[i];
        z[i]=x0[i];
      }

      set_theta();
      cout<<"set stepsizes"<<endl;

      current_nb_iters=0;
      compute_primal_value();
      cout<<"primal value="<<primal_value<<endl;
      cout<<"Initialization  loopless Katyusha is finished!"<<endl;
    }


void compute_primal_value()
{
  D res=0;
  for(L i=0;i<nsamples;i++)
  {
    res+=value_of_phi_i(i);
  }
  D res2=value_of_g();
  primal_value= res/nsamples+ res2;
}


void compute_and_record_result()
{
  if(nb_loops%print_every_N==0){
    compute_primal_value();
    cout<<setprecision(9)<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<";  "<<running_time<<" primal value="<<primal_value<<endl;
    samp<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<" "<<running_time<<" "<<primal_value<<endl;
  }
}

void compute_and_record_result2()
{
  if(nb_loops%print_every_N==0){
    compute_primal_value();
    compute_upper_bound_of_optimality_gap();
    cout<<setprecision(9)<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<";  "<<running_time<<" primal value="<<primal_value<<" upper bound of F(x)-F^*= "<< upper_bound<< endl;
    samp<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<" "<<upper_bound<<" "<<running_time<<" "<<primal_value<<endl;
  }

}


void batch_sampling()
{
  batch_size=0;
  L i;
  for(L it_t=0;it_t<tau;it_t++)
  {
    i=sampling(nsamples);
    batch_i[batch_size]=i;
    batch_size++;
  }
}


void set_thetaS(L p_mod){
  theta_S.clear();
  theta_S.resize(nsamples,0);
  if(p_mod==0){
    for(L i=0;i<nsamples;i++)
    theta_S[i]=1.0/(tau*tilde_proba_vector[i]);
  }
  else{
    for(L i=0;i<nsamples;i++)
    theta_S[i]=1.0/proba_vector[i];
  }
}

void group_sampling(){
  batch_size=0;
  std::vector<L> sampled_groups;
  for(L i=0;i<nb_groups;i++){
    D y=gsl_rng_uniform(rng);
    if(y<sump_group_C[i])
    {
      sampled_groups.push_back(i);
      batch_size++;
    }
  }
  for(L t=0;t<batch_size;t++){
    L group_i=sampled_groups[t];
    L s1=index_group_C[group_i];
    L s2=index_group_C[group_i+1];
    L i=s1+(floor)(gsl_rng_uniform(rng)*(s2-s1));
    D y=gsl_rng_uniform(rng);
    D maxpi=maxp_group_C[group_i];
    while(y*maxpi>proba_vector[group_C[i]])
    {
      i=s1+(floor)(gsl_rng_uniform(rng)*(s2-s1));
      y=gsl_rng_uniform(rng);
    }
    batch_i[t]=group_C[i];
  }
}




void loopless(vector<D> & x0, string filename, D val_mu, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
{
  initialize(x0, val_mu, max_nb,  nb_tau, nb_c,  u, p_mod, scal_p);
  string sampname="results/L_Katyusha"+filename+uniform;
  if(p_mod==0) sampname=sampname+"_batch";
  else sampname=sampname+"_group";
  string scal_str;
  stringstream  scal_convert;
  scal_convert<<scal_p;
  scal_str=scal_convert.str();
  sampname+=scal_str;
  cout<<"running Loopless Katyusha "<<" ; "<<sampname<<endl;
  samp.open(sampname.c_str());

  nb_loops=0;
  nb_iters=0;

  srand48(27432042);
  //srand(time(NULL));


  D start;
  while(nb_loops<max_nb_loops)
  {
    compute_and_record_result();

    nb_loops++;

    start = std::clock();

    for(L it=0;it<nb_of_iters_per_loop;it++)
    {
      if(p_mod==0)
      batch_sampling();
      else
      group_sampling();

      start = std::clock();
      update_x();
      update_gk();
      update_z();
      update_y();
      update_w();
      update_gradient();
      current_nb_iters++;
      nb_iters++;
      running_time+= ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    }
  }
}

void loopless2(vector<D> & x0, string filename, D val_mu, L max_nb, D epsilon, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
{
  cout<<"max_nb="<<max_nb<<endl;
  initialize(x0, val_mu, max_nb,  nb_tau, nb_c,  u, p_mod, scal_p);
  string sampname="results/L_"+filename+uniform;
  if(p_mod==0) sampname=sampname+"_batch";
  else sampname=sampname+"_group";
  string scal_str;
  stringstream  scal_convert;
  scal_convert<<scal_p;
  scal_str=scal_convert.str();
  sampname+=scal_str;
  cout<<"running Loopless Katyusha 2"<<" ; "<<sampname<<endl;
  samp.open(sampname.c_str());

  nb_loops=0;
  nb_iters=0;

  cout<<setprecision(9)<<"initial: "<< "primal: "<<primal_value << endl;
  samp<<((0.0+nb_iters)*c*tau/(nsamples))<<" "<<running_time<<" "<<primal_value<<endl;
  //srand48(27432042);
  srand(time(NULL));


  D start;
  cout<<"nb_loops="<<nb_loops<<"; max_nb_loops="<<max_nb_loops<<"; upper bound="<<upper_bound<<"; epsilon="<<epsilon<<endl;
  while(nb_loops<max_nb_loops && upper_bound> epsilon)
  {
    compute_and_record_result2();
    start = std::clock();

    nb_loops++;

    //cout<<"before the loop time elapsed="<<duration<<endl;
    start = std::clock();

    //cout<<"nb_of_iters_per_loop="<<nb_of_iters_per_loop<<endl;
    for(L it=0;it<nb_of_iters_per_loop;it++)
    {
      if(p_mod==0)
      batch_sampling();
      else
      group_sampling();

      start = std::clock();
      update_x();
      update_gk();
      update_z();
      update_y();
      update_w();
      update_gradient();
      current_nb_iters++;
      //cout<< "test 1"<< endl;
      //cout<< "test 2"<< endl;
      nb_iters++;
      running_time+= ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    }
    //cout<< "epoch "<< nb_loops<< " finished"<< endl;
    //cout<<"after the loop time elapsed="<<duration<<endl;
  }
}








};

#endif /* MIN_SMOOTH_CONVEX_H */
