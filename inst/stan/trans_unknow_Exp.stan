data{
  int<lower=1> p;
  int<lower=1> num_basis; // the degree of splines basis
  int<lower=1> K; // DP
  real<lower=0> alpha; //DP constant
  //int<lower=1> N; // number of time grinds
  int<lower=1> n; //numer of observations in this cluster
  int<lower=1> NN;
  matrix[num_basis, NN] Yspline;
  matrix[num_basis, NN] yspline;

  vector[n] Y;
  int location[n];
  vector[n] delta;
  matrix[n,p] X;
  real<lower = 0> eta;
  real<lower = 0> psi;
  real<lower = 0> v;
  vector[p] beta_init0;
  vector<lower=0, upper=1>[K] w_init0;
}

parameters {
  vector[p] beta;
  //real<lower=0> alpha;
  //real mu;
  vector<lower=0>[K] nu;
  vector<lower=0>[K] theta;
  vector<lower=0, upper=1>[K] w;


  row_vector<lower=0>[num_basis] a_raw;
}

transformed parameters{
  //row_vector [num_basis] tau_raw;


  simplex[K] DP_weights;

  DP_weights[1] = w[1];
  for (s in 2:(K-1))  { DP_weights[s] = w[s] * prod(1 - w[1:(s - 1)]);}
  DP_weights[K] = 1- sum(DP_weights[1:(K-1)]);
  //tau_raw = exp(a_raw);
}

model{
   beta ~ normal(beta_init0, 1000);
  //beta ~ cauchy(0, 2.5);
  nu ~ gamma(psi, v);
  theta ~ gamma(psi, v);
  //alpha ~ gamma(1, 1);
  //theta ~ gamma(.5, 1);
  w ~ beta(1, alpha);
  a_raw ~ exponential(eta);
  //a_raw ~ inv_gamma(.01, .01);
  //a_raw ~ exponential(1);

  for (i in 1:n) {
    real H_Y=0;
    real DH_Y=0;
    vector[K] lp_ik;
    H_Y =  a_raw  * to_vector(Yspline[:, location[i]]);
    DH_Y =  a_raw  * to_vector(yspline[:, location[i]]);
    //print(H_Y, " ",DH_Y);
    for(k in 1:K){
      if (delta[i]>0)
      {lp_ik[k] = log(DP_weights[k]) + weibull_lpdf(H_Y * exp((-X[i, 1:p] * beta))|nu[k], theta[k]) + log(DH_Y)
      -(X[i,1:p] * beta);}
      else
      {lp_ik[k] = log(DP_weights[k]) + weibull_lccdf(H_Y* exp((-X[i, 1:p]*beta))|nu[k], theta[k]);}

    }
    //print(lp_ik);
    target += log_sum_exp(lp_ik);
  }
}

generated quantities{
  vector[n] log_lik;
  vector[p] beta_trans;
  beta_trans = beta ./ sqrt(sum(beta .* beta));


for (qq in 1:n) {
    real H_Y=0;
    real DH_Y=0;
    vector[K] lp_kk;
    H_Y =  a_raw  * to_vector(Yspline[:, location[qq]]);
    DH_Y =  a_raw  * to_vector(yspline[:, location[qq]]);
    //print(H_Y, " ",DH_Y);
    for(k in 1:K){
      if (delta[qq]>0)
      {lp_kk[k] = log(DP_weights[k]) + weibull_lpdf(H_Y * exp((-X[qq, 1:p] * beta))|nu[k], theta[k]) + log(DH_Y)
      -(X[qq,1:p] * beta);}
      else
      {lp_kk[k] = log(DP_weights[k]) + weibull_lccdf(H_Y* exp((-X[qq, 1:p]*beta))|nu[k], theta[k]);}

    }

  log_lik[qq] = log_sum_exp(lp_kk);
}
}


