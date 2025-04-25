//--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes
// partial pooling at province-district level

data {
  
  int nA; // N age groups
  int nT; // N time points
  int nD; // N admin2 within a province
  array[nD, nT, nA] int cases; // reported case data
  array[nD, nT, nA] int pop; // population data
  array[2,nA] int ageLims; // lower & upper bounds of age groups
  int hist_length; // length of historical lambda
  
}

parameters {
  
  // Province level parameters
  real<lower=0,upper=0.15> lam_H_prov_constant;  // constant historical FOI
  array[10] real<lower=0,upper=0.15> lam_H_prov_recent;  // recent historical FOI
  array[nT] real<lower=0,upper=0.15> lam_t_prov;
  real<lower=0,upper=1> rho_prov;
  real<lower=0,upper=1> gamma_prov;
  
  // District level parameters
  array[nD] real<lower=0,upper=0.15> lam_H_constant;  // district constant historical FOI
  array[nD,10] real<lower=0,upper=0.15> lam_H_recent; // district recent historical FOI
  array[nD] vector<lower=0,upper=0.15>[nT] lam_t;
  array[nD] real<lower=0,upper=1> rho;
  array[nD] real<lower=0,upper=1> gamma;
  
  // Kappa parameters
  real<lower=0> kappa_lam_H_constant;
  real<lower=0> kappa_lam_H_recent;
  real<lower=0> kappa_lam_t;
  real<lower=0> kappa_rho;
  real<lower=0> kappa_gamma;

}

transformed parameters {
  
  array[nD,hist_length] real<lower=0,upper=0.15> lam_H;
  array[nD] matrix<lower=0,upper=1>[nT+1,100] susc; // proportion susceptible
  array[nD] matrix<lower=0,upper=1>[nT+1,100] mono; // proportion monotypic
  array[nD] matrix<lower=0,upper=1>[nT,100] inc1; // incidence of primary infections
  array[nD] matrix<lower=0,upper=1>[nT,100] inc2; // incidence of secondary infections
  array[nD] matrix<lower=0>[nT,nA] Ecases; // expected reported cases
  array[nD] real sum_lamH;
  
  // Construct full historical FOI for each district
  for(d in 1:nD) {
    // Recent 10 years first
    for (t in 1:10) {
      lam_H[d,t] = lam_H_recent[d,t];
    }
    // Earlier constant period fills the rest
    for (t in 11:hist_length) {
      lam_H[d,t] = lam_H_constant[d];
    }
  }
  
  //--- immune profiles at beginning of time series
  for(d in 1:nD) for (i in 1:100){
    
    if (i > hist_length) {
      sum_lamH[d] = sum(lam_H[d][1:hist_length]);
    } else {
      // lam_H_vec[i] = lam_H;
      sum_lamH[d] = sum(lam_H[d][1:i]);
    }
    
    susc[d][1,i] = exp(-4*sum_lamH[d]);
    mono[d][1,i] = 4*exp(-3*sum_lamH[d])*(1-exp(-sum_lamH[d]));
      
  }
  
  //--- subsequent time steps
  for(d in 1:nD) for (t in 2:(nT+1)){
    
    susc[d][t,1] = exp(-4*lam_t[d][t-1]);
    mono[d][t,1] = 4*exp(-3*lam_t[d][t-1])*(1-exp(-lam_t[d][t-1]));
    susc[d][t,2:100] = susc[d][t-1,1:99] - 4*lam_t[d][t-1]*susc[d][t-1,1:99];
    mono[d][t,2:100] = mono[d][t-1,1:99] + 4*lam_t[d][t-1]*susc[d][t-1,1:99] - 3*lam_t[d][t-1]*mono[d][t-1,1:99];
    inc1[d][t-1,1] = 4*lam_t[d][t-1]*1;
    inc2[d][t-1,1] = 3*lam_t[d][t-1]*0;
    inc1[d][t-1,2:100] = 4*lam_t[d][t-1]*susc[d][t-1,1:99];
    inc2[d][t-1,2:100] = 3*lam_t[d][t-1]*mono[d][t-1,1:99];
    
  }
  
  //--- reported cases
  // index of age + 1 as the age matrix starts from 0 and ends at 99
  for(d in 1:nD) for(t in 1:nT) for(a in 1:nA){
    Ecases[d][t,a] = rho[d]*(mean(inc2[d][t,(ageLims[1,a]+1):(ageLims[2,a]+1)]) + gamma[d]*mean(inc1[d][t,(ageLims[1,a]+1):(ageLims[2,a]+1)]))*pop[d][t,a];
  }

}

model {

  //--- priors
  // Province level
  lam_H_prov_constant ~ normal(0,0.05);
  lam_H_prov_recent ~ normal(0,0.05);
  lam_t_prov ~ normal(0,0.05);
  rho_prov ~ normal(0,0.15);
  gamma_prov ~ normal(0,1);
  
  // Kappa parameters
  kappa_lam_H_constant ~ normal(0,1);
  kappa_lam_H_recent ~ normal(0,1);
  kappa_lam_t ~ normal(0,1);
  kappa_rho ~ normal(0,1);
  kappa_gamma ~ normal(0,1);
  
  // District level with hierarchical structure
  for(d in 1:nD) {
    // Historical FOI
    lam_H_constant[d] ~ normal(lam_H_prov_constant, kappa_lam_H_constant);
    for(t in 1:10) {
      lam_H_recent[d,t] ~ normal(lam_H_prov_recent[t], kappa_lam_H_recent);
    }
    
    // Current FOI and reporting
    for(t in 1:nT) {
      lam_t[d][t] ~ normal(lam_t_prov[t], kappa_lam_t);
    }
    rho[d] ~ normal(rho_prov, kappa_rho);
    gamma[d] ~ normal(gamma_prov, kappa_gamma);
  }
  
  //--- likelihood 
  for(d in 1:nD) for(t in 1:nT) {
    cases[d][t,] ~ poisson(Ecases[d][t,]); // poisson likelihood
  }
  
}

generated quantities {
  array[nD, nT, nA] real log_lik;  // log likelihoods
  array[nD, nT, nA] int y_rep;     // posterior predictive samples
  // Add vector version of log_lik for LOO
  vector[nD * nT * nA] log_lik_vector;
  int i = 1;
  
  for(d in 1:nD) {
    for(t in 1:nT) {
      for(a in 1:nA) {
        log_lik[d,t,a] = poisson_lpmf(cases[d,t,a] | Ecases[d,t,a]);
        y_rep[d,t,a] = poisson_rng(Ecases[d,t,a]);
        log_lik_vector[i] = log_lik[d,t,a];
        i += 1;
      }
    }
  }
}
