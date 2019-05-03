data {
  int<lower=0> nrows;         // number of rows of data
  int<lower=0> nbugs;         // number of number of indicator species
  //real Coliforms[nrows]; 
  real Ecoli[nrows]; 
  int Recipe[nrows]; 
  int Time[nrows]; 
  int OL[nrows]; 
  real VS[nrows]; 
  real Ammonia[nrows]; 
  int Temp[nrows]; 
  real Methane[nrows]; 
  //vector[nrows] y;                         // pathogen load

  //real Enterococci[nrows]; 
}
parameters {
  real recipe_effect[3] ;
  //real load_effect[3] ;
  real temp_effect[3] ;
  real ol_effect[3] ;
  // ### error in coliforms, E. coli, Entrerococci
  //real<lower=0> sigma_co[nrows]; // standard error of effect estimates
  //real<lower=0> sigma_ec[nrows]; // standard error of effect estimates
  //real<lower=0> sigma_en[nrows]; // standard error of effect estimates
  real<lower=0> sigma; // measurement error
  //real<lower=0> eta; // measurement error
 }
transformed parameters {
  vector[nrows] bugs_hat; // the *real* pathogen load, without measurement error
  for (i in 1:nrows){
    bugs_hat[i] = ( temp_effect[Temp[i]] * Time[i]) + (recipe_effect[Recipe[i]] * ol_effect[OL[i]]) + (VS[i] * recipe_effect[Recipe[i]]);
  }
  // bugs_hat[i] = (recipe_effect[Recipe[i]] * Time[i]) + load_effect[Recipe[i]] + temp_effect[Temp[i]] + (recipe_effect[Recipe[i]] * OL[i]) + (VS[i] * recipe_effect[Recipe[i]]);
 // for (i in 1:nrows){
 //    bugs_hat[i] = (Time[i]) + load_effect[Recipe[i]] + temp_effect[Temp[i]] + ol_effect[Recipe[i]] * OL[i] + VS[i];
 //  }
}
model {
  // for (i in 1:nrows){
  //   ybugs[i] = cauchy(log_indicators[i, 1] + log_indicators[i, 2] + log_indicators[i, 3], )
  // }
  //eta ~ cauchy(1, 1);
  sigma ~ cauchy(0, 1);
  recipe_effect ~ cauchy(1, 1);
  //load_effect ~ cauchy(1, 1);
  temp_effect ~ cauchy(1, 1);
  //sigma_co ~ cauchy(sigma, eta);
  //sigma_ec ~ cauchy(sigma, eta);
  //sigma_en ~ cauchy(sigma, eta);
  //for(i in 1:nrows){
  Ecoli ~ normal(bugs_hat, sigma);
  //}
}

