data{
  int<lower=1> p; // covariate basis
  int<lower=1> num_basis; // the number of splines basis
  int<lower=1> K; // DP truncation number
  real<lower=0> a; //DP constant
  //int<lower=1> N; // number of time grinds
  int<lower=1> n; // size of training set
  int<lower=1> NN; // number of grids of spline functions
  matrix[num_basis, NN] Yspline;
  matrix[num_basis, NN] yspline;
  vector[n] Y;
  int location[n];
  vector[n] delta;
  matrix[n,p] X;
  real<lower = 0> eta;
  real<lower = 0> zeta;
  real<lower = 0> v;
  vector[p] beta_init0;
  vector<lower=0, upper=1>[K] w_init0;
}

parameters {
  vector[p] beta;
  //real<lower=0> alpha;
  //real mu;
  vector<lower=0>[K] psi;
  vector<lower=0>[K] nu;
  vector<lower=0, upper=1>[K] w;

  row_vector<lower=0>[num_basis] alpha;
}


transformed parameters{

  simplex[K] DP_weights;

  DP_weights[1] = w[1];
  for (s in 2:(K-1))  { DP_weights[s] = w[s] * prod(1 - w[1:(s - 1)]);}
  DP_weights[K] = 1- sum(DP_weights[1:(K-1)]);
  //tau_raw = exp(alpha);
}

model{
  beta ~ normal(0, 1000);
  psi ~ exponential(zeta);
  nu ~ exponential(v);
  w ~ beta(1, a);
  alpha ~ exponential(eta);

  for (i in 1:n) {
    real H_Y=0;
    real DH_Y=0;
    vector[K] lp_ik;
    H_Y =  alpha  * to_vector(Yspline[:, location[i]]);
    DH_Y =  alpha  * to_vector(yspline[:, location[i]]);
    //print(H_Y, " ",DH_Y);
    for(k in 1:K){
      if (delta[i]>0)
      {lp_ik[k] = log(DP_weights[k]) + weibull_lpdf(H_Y * exp((-X[i, 1:p] * beta))|psi[k], nu[k]) + log(DH_Y)
      -(X[i,1:p] * beta);}
      else
      {lp_ik[k] = log(DP_weights[k]) + weibull_lccdf(H_Y* exp((-X[i, 1:p]*beta))|psi[k], nu[k]);}

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
    H_Y =  alpha  * to_vector(Yspline[:, location[qq]]);
    DH_Y =  alpha  * to_vector(yspline[:, location[qq]]);
    //print(H_Y, " ",DH_Y);
    for(k in 1:K){
      if (delta[qq]>0)
      {lp_kk[k] = log(DP_weights[k]) + weibull_lpdf(H_Y * exp((-X[qq, 1:p] * beta))|psi[k], nu[k]) + log(DH_Y)
      -(X[qq,1:p] * beta);}
      else
      {lp_kk[k] = log(DP_weights[k]) + weibull_lccdf(H_Y* exp((-X[qq, 1:p]*beta))|psi[k], nu[k]);}
    }
  log_lik[qq] = log_sum_exp(lp_kk);
}

}


