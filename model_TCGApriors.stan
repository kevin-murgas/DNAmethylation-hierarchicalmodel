// DNA Methylation Hierarchical Variation Model
// with TCGA-data-generated priors
// Code By Kevin Murgas and Yanlin Ma
// Working with Dr. Marc Ryser, Dr. Darryl Shibata

data {
  // This part defines variables used in model section
  int<lower=1> N;   // number of observations
  int<lower=1> P;   // number of patients
  real y[N];   // logit-transformed beta values
  int<lower=0, upper=1> tInd[N];  // tumor indicator
  int<lower=1, upper=P> pID[N];   // Patient ID
}

parameters {
  real mu;      // grand mean
  real nu;      // fixed effect of tumor tissue
  
  real<lower=0> sigma_p;  // patient sd
  real<lower=0> sigma_pt; // tumor sd at patient level
  real<lower=0> sigma_t;  // intratumor sd
  real<lower=0> sigma_e;  // error sd
  
  vector[P] a_pat;  // random  effect of patients in normal
  vector[P] b_patT; // random effect of tumor tissue (varying between patients)
  vector[N] g_T;    // random effect of tumor tissue (intratumor)
}

model {
  vector[N] summand;
  // priors (based on TCGA empirical data fits)
  // mu is coded below using a bimodal mixture model
  target += log_mix(0.3557,normal_lpdf(mu|-3.0675, 1.3027) , normal_lpdf(mu|1.1290, 2.4110));
  nu ~ cauchy(-0.04051, 0.5799);
  
  sigma_p ~ gamma(2.7693, 7.0029);
  sigma_pt ~ gamma(2.1457, 2.4819);
  sigma_t ~ gamma(3.3643, 10.430);
  sigma_e ~ gamma(2,2); //non-informative gamma for sigma_e
  
  // random effect coefficient parameter distributions
  a_pat ~ normal(0,sigma_p);
  b_patT ~ normal(0,sigma_pt);
  g_T ~ normal(0,sigma_t);
  
  //posteriors
  for (n in 1:N){
    if (tInd[n]) {
      summand[n] = mu + a_pat[pID[n]] + nu + b_patT[pID[n]] + g_T[n];
    } else {
      summand[n] = mu + a_pat[pID[n]];
    }
  }
  y ~ normal(summand, sigma_e);
}
