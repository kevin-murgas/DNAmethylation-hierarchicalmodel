// This model introduces variation within tumor tissue

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
  real mu;    // grand mean
  real betaT;   // fixed effect of tumor tissue
  
  real<lower=0> sigma_p;  // patient sd
  real<lower=0> sigma_pt; // tumor sd at patient level
  real<lower=0> sigma_t;  // intratumor sd
  real<lower=0> sigma_e;  // error sd
  
  vector[P] b_pat;   /// random  effect of patients in normal
  vector[P] c_patT;  // random effect of tumor tissue (varying between patients)
  vector[N] d_T; // random effect of tumor tissue (intratumor)
}

model {
  //priors
  mu ~ normal(0,100000);
  betaT ~ normal(0,1.24242);
  
  sigma_p ~ gamma(4,8);
  sigma_pt ~ gamma(4,8);
  sigma_t ~ gamma(4,8);
  sigma_e ~ uniform(0,100);
  
  b_pat ~ normal(0,sigma_p);
  c_patT ~ normal(0,sigma_pt);
  d_T ~ normal(0,sigma_t);
  
  //posteriors
  for (n in 1:N){
    if (tInd[n]==0) {
      y[n] ~ normal(b_pat[pID[n]] + mu, sigma_e);
    } else {
      y[n] ~ normal(b_pat[pID[n]] + betaT + c_patT[pID[n]] + d_T[n] + mu, 
      sigma_e);
    }
  }
}
