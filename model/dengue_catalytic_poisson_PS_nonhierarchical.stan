//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
////--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes

data {
  
  int nA; // N age groups
  int nT; // N time points
  array[nT,nA] int cases; // reported case data
  matrix[nT,nA] pop; // population data
  array[2,nA] int ageLims; // lower & upper bounds of age groups
  int hist_length; // length of historical lambda
  // row_vector[100] age;
  
}


parameters {
  
  real<lower=0,upper=0.15> lam_H_constant;  // constant historical FOI
  array[10] real<lower=0,upper=0.15> lam_H_recent; // recent historical FOI
  array[nT] real<lower=0,upper=0.15> lam_t; // time varying FOI
  real<lower=0,upper=1> rho;
  real<lower=0,upper=1> gamma;

}

transformed parameters {
  
  array[hist_length] real<lower=0,upper=0.15> lam_H;
  matrix<lower=0,upper=1>[nT+1,100] susc; // proportion susceptible
  matrix<lower=0,upper=1>[nT+1,100] mono; // proportion monotypic
  matrix<lower=0,upper=1>[nT,100] inc1; // incidence of primary infections
  matrix<lower=0,upper=1>[nT,100] inc2; // incidence of secondary infections
  array[nT] vector<lower=0>[nA] Ecases; // expected reported cases
  real sum_lamH;
  
  // Fill historical FOI array
  // Recent 10 years first
  for (t in 1:10) {
    lam_H[t] = lam_H_recent[t];
  }
  // Earlier constant period fills the rest
  for (t in 11:hist_length) {
    lam_H[t] = lam_H_constant;
  }
  
  //--- immune profiles at beginning of time series
  for (i in 1:100){
    
    if (i > hist_length) {
      sum_lamH = sum(lam_H[1:hist_length]);
    } else {
      // lam_H_vec[i] = lam_H;
      sum_lamH = sum(lam_H[1:i]);
    }
    
    susc[1,i] = exp(-4*sum_lamH);
    mono[1,i] = 4*exp(-3*sum_lamH)*(1-exp(-sum_lamH));
      
  }
  
  //--- subsequent time steps
  for (t in 2:(nT+1)){
    
    susc[t,1] = exp(-4*lam_t[t-1]);
    mono[t,1] = 4*exp(-3*lam_t[t-1])*(1-exp(-lam_t[t-1]));
    susc[t,2:100] = susc[t-1,1:99] - 4*lam_t[t-1]*susc[t-1,1:99];
    mono[t,2:100] = mono[t-1,1:99] + 4*lam_t[t-1]*susc[t-1,1:99] - 3*lam_t[t-1]*mono[t-1,1:99];
    inc1[t-1,1] = 4*lam_t[t-1]*1;
    inc2[t-1,1] = 3*lam_t[t-1]*0;
    inc1[t-1,2:100] = 4*lam_t[t-1]*susc[t-1,1:99];
    inc2[t-1,2:100] = 3*lam_t[t-1]*mono[t-1,1:99];
    
  }
  
  //--- reported cases
  // index of age + 1 as the age matrix starts from 0 and ends at 99
  for(t in 1:nT) for(a in 1:nA){
    Ecases[t,a] = rho*(mean(inc2[t,(ageLims[1,a]+1):(ageLims[2,a]+1)]) + gamma*mean(inc1[t,(ageLims[1,a]+1):(ageLims[2,a]+1)]))*pop[t,a];
  }

}

model {
  
  //--- priors
  lam_H_constant ~ normal(0,0.05);
  lam_H_recent ~ normal(0,0.05);
  lam_t ~ normal(0,0.05);
  rho ~ normal(0,0.15);
  gamma ~ normal(0,1);
  
  //--- likelihood 
  for(t in 1:nT) cases[t,] ~ poisson(Ecases[t,]); // poisson likelihood
  
}

generated quantities {
  array[nT, nA] real log_lik;  // log likelihoods
  array[nT, nA] int y_rep;     // posterior predictive samples
  // Add vector version of log_lik for LOO
  vector[nT * nA] log_lik_vector;
  int i = 1;
  
  for(t in 1:nT) {
    for(a in 1:nA) {
      log_lik[t,a] = poisson_lpmf(cases[t,a] | Ecases[t,a]);
      y_rep[t,a] = poisson_rng(Ecases[t,a]);
      log_lik_vector[i] = log_lik[t,a];
      i += 1;
    }
  }
}
