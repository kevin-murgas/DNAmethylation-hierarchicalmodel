// DNA Methylation Hierarchical Variation Model
// with TCGA-data-generated priors
// Code By Kevin Murgas and Yanlin Ma
// Working with Dr. Marc Ryser, Dr. Darryl Shibata

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
  vector[N] summand;
  //priors 
  // mu is coded below using a bimodal mixture model
  target += log_mix(0.3557,normal_lpdf(mu|-3.0675, 0.75214) , normal_lpdf(mu|1.1290, 1.39197));
  betaT ~ cauchy(-0.04051, 0.19330);
  
  sigma_p ~ inv_gamma(6.3641, 1.6631);
  sigma_pt ~ inv_gamma(3.2624, 1.5397);
  sigma_t ~ inv_gamma(6.5392, 1.4894);
  sigma_e ~ inv_gamma(1,1); //non-informative inv-gamma for sigma_e
  
  b_pat ~ normal(0,sigma_p);
  c_patT ~ normal(0,sigma_pt); 
  d_T ~ normal(0,sigma_t);
  
  //posteriors
  for (n in 1:N){
    summand[n] = mu + b_pat[pID[n]] + tInd[n]*(betaT + c_patT[pID[n]] + d_T[n]);
  }
  y ~ normal(summand,sigma_e);
}
