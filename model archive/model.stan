data {
  
  // This part defines variables used in model section
  
  int<lower=1> N;   // number of observations
  int<lower=1> P;   // number of patients
  
  
  real y[N];   // transformed beta values
  int<lower=0, upper=1> tInd[N];   // tumor indicator
  
 
  int<lower=1, upper=P> pID[N];   // Patient ID
 }



parameters {
  
  // This part defines parameters we plan to use
  
  vector[P] b_pat;   /// random  effect of patients
  vector[P] bT_pat;  // random effect of tumor tissue (varying between patients)
  real betaT;   // fixed effect of tumor tissue
  real mu;    // grand mean
  
  real<lower=0> sigma_e;  // error sd
  real<lower=0> sigma_p;  // patient sd
  real<lower=0> sigma_t;  // tissue sd
  

}



model {
  //prior
  // uniform prior on mu, betaT, sigma_t, sigma_e and sigma_p when omitting assigning prior distribution on them (will change these later)
  
  b_pat ~ normal(0,sigma_p); 
  bT_pat ~ normal(0,sigma_t);
  
  
  for (n in 1:N){
    if (tInd[n]==0) {
      y[n] ~ normal(b_pat[pID[n]] + mu, sigma_e);
    } else {
      y[n] ~ normal(b_pat[pID[n]] + betaT + bT_pat[pID[n]] + mu, sigma_e);
    }
  }

}