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
    vector[nPops] B;
    vector<lower=0, upper=1>[nPops] nSurv;
    real <lower=0>tauB;
    real <lower=0>tauA;
    real <lower = 0> tauSurv;
    real <lower=0>muB;
    //real muB;
    real <lower=0>muA;
    real <lower=0>muSurv;
}
model{
    vector[N] p;
    tauB ~ normal( 0, 5 );
    tauA ~ normal(0,5);
    tauSurv ~ normal(0,5);

    muB ~ normal(0,5);
    muA ~ normal(0, 1);
    muSurv ~ exponential(2); //normal(0.5,0.5)
    B ~ normal( muB , tauB );
    a ~ normal( muA, tauA );
    nSurv ~ normal(muSurv, tauSurv);

    for ( i in 1:N ) {
        if (scaledDose[i] == 0.0){
          p[i] = nSurv[Population[i]];
        } else { 
          if (uselogit){
            p[i] = inv_logit((log10(a[Population[i]]) - log10(scaledDose[i])) * B[Population[i]]) * nSurv[Population[i]];
          } else {
            p[i] = Phi((log10(a[Population[i]]) - log10(scaledDose[i])) * B[Population[i]]) * nSurv[Population[i]];
          }
        }
    }
    surv ~ binomial( tested , p );
}

generated quantities{
    vector[N] log_lik;
    vector[N] pp;
    vector[nPops] lc50;
    
    for ( i in 1:(nPops) ) {
      lc50[i] = a[i] * correction;
    }
    
    for ( i in 1:(N) ) { 
       if (scaledDose[i] == 0.0){
         pp[i] = nSurv[Population[i]];
       } else {
         if (uselogit){
           pp[i] = inv_logit((log10(a[Population[i]]) - log10(scaledDose[i])) * B[Population[i]]) * nSurv[Population[i]];
         } else {
           pp[i] = Phi((log10(a[Population[i]]) - log10(scaledDose[i])) * B[Population[i]]) * nSurv[Population[i]];
         }
      }
    }
    
    for ( i in 1:(N) ) {
         log_lik[i] = binomial_lpmf( surv[i] | tested[i] , pp[i] );
    }
}
