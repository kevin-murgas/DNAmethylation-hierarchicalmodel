// Model 1

data {
  
  // This part defines variables used in model section
  
  int<lower=1> N;   // number of observations
  int<lower=1> P;   // number of patients
  
  
  real y[N];   // transformed beta values
  
  int<lower=1, upper=P> pID[N];   // Patient ID
}



parameters {
  
  // This part defines parameters we plan to use
  
  vector[P] b_pat;   /// random  effect of patients

  vector[N] b_T; // random effect of tumor tissue 
  
  real mu;    // grand mean
  
  real<lower=0> sigma_e;  // error sd
  real<lower=0> sigma_p;  // patient sd
  real<lower=0> sigma_t;  
  
}



model {
  //prior
  
  sigma_t ~ uniform(0,100);
  sigma_e ~ uniform(0,100);
  sigma_p ~ uniform(0,100);
 
  b_pat ~ normal(0,sigma_p);
  b_T ~ normal(0,sigma_t);
  
  mu ~ normal(0,100000);
  
  
  
  //posterior

  for (n in 1:N) {
  y[n] ~ normal(b_pat[pID[n]] + b_T[pID[n]] + mu, sigma_e);
  }
  
}
