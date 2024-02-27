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
    vector<lower=0, upper=2000>[nPops] phiB;
}

transformed parameters{
    vector<lower=2>[nPops] theta;
    theta = 2. + phiB;
}

model{
    vector[N] p;

    B ~ normal( 3 , 3 );
    a ~ normal(0.5, 1 );
    nSurv ~ normal(0.5, 1);

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
        surv[i] ~ beta_binomial( tested[i] , p[i]*theta[Population[i]] , (1-p[i])*theta[Population[i]] );
    }
}

generated quantities{
    vector[N] log_lik;
    vector[N] p;
    vector[nPops] lc50;
    
    for ( i in 1:(nPops) ) {
      lc50[i] = a[i] * correction;
    }
    
    for ( i in 1:(N) ) {
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
    

        for ( i in 1:(N) ) {
     log_lik[i] = beta_binomial_lpmf( surv[i] | tested[i] , p[i]*theta[Population[i]] , (1-p[i])*theta[Population[i]] );
    }
  
}
