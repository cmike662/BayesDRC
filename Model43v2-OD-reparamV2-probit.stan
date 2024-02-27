data{
    int N;
    int nPops;
    int nBlocks;
    array[N] int block;
    array[N] int tested;
    array[N] int surv;
    vector[N] scaledDose;
    array[N] int Population;
    real correction;
    int uselogit;
}
parameters{
    vector<lower=0>[nPops] a;
    vector[nPops] zB;
    vector<lower=0, upper=1>[nPops] nSurv;
    vector<lower=0, upper=2000>[nPops] phiB;
    real <lower=0>tauB;
    real <lower=0>tauA;
    real <lower=0> tauSurv;
    real <lower=0>muB;
    real <lower=0>muA;
    real <lower=0>muSurv;
}

transformed parameters{
    vector<lower=2.0>[nPops] theta;
    theta = 2. + phiB;
}

model{
    vector[N] p;
    tauB ~ normal(0,5);
    tauA ~ normal(0,5);
    tauSurv ~ normal(0,5);
    muB ~ normal(0,5);
    muA ~ normal(0, 5);
    muSurv ~ exponential(2); 
    zB ~ normal(0,1);
    a ~ normal( muA, tauA );
    //B = (muB+zB[Population[i]]*tauB);
    nSurv ~ normal(muSurv, tauSurv);

    for ( i in 1:N ) {
        if (scaledDose[i] == 0.0){
          p[i] = nSurv[Population[i]];
        } else { 
          if (uselogit){
            p[i] = inv_logit((log10(a[Population[i]]) - log10(scaledDose[i])) * (muB+zB[Population[i]]*tauB)) * nSurv[Population[i]];
          } else {
            p[i] = Phi((log10(a[Population[i]]) - log10(scaledDose[i])) * (muB+zB[Population[i]]*tauB)) * nSurv[Population[i]];
          }
        }
        surv[i] ~ beta_binomial( tested[i] , p[i]*theta[Population[i]] , (1-p[i])*theta[Population[i]] );
    }
}

generated quantities{
    vector[N] log_lik;
    vector[N] p;
    vector[nPops] B;
    vector[nPops] lc50;
    
    for ( i in 1:(nPops) ) {
      B[i] = (muB+zB[i]*tauB);
      lc50[i] = a[i] * correction;
    }
    
    
     for ( i in 1:(N) ) {
       if (scaledDose[i] == 0.0){
         p[i] = nSurv[Population[i]];
       } else {
         if (uselogit){
           p[i] = inv_logit((log10(a[Population[i]]) - log10(scaledDose[i])) * (muB+zB[Population[i]]*tauB)) * nSurv[Population[i]];
         } else {
           p[i] = Phi((log10(a[Population[i]]) - log10(scaledDose[i])) * (muB+zB[Population[i]]*tauB)) * nSurv[Population[i]];
         }
      }
    }
    
    for ( i in 1:(N) ) {
      log_lik[i] = beta_binomial_lpmf( surv[i] | tested[i] , p[i]*theta[Population[i]] , (1-p[i])*theta[Population[i]] );
    }
}
