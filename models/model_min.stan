data {
  int<lower=0> nrows;         // number of rows of data
  int<lower=0> nbugs;         // number of number of indicator species
  real Ecoli[nrows]; 
  real Enterococci[nrows]; 
  int Recipe[nrows]; 
  int Time[nrows]; 
  real sCOD[nrows]; 
  int OL[nrows]; 
  real VS[nrows]; 
  real Ammonia[nrows]; 
  int Temp[nrows]; 
  real Methane[nrows]; 
}
parameters {
  real recipe_effect_m[3] ;
  real temp_effect_m[3] ;
  real ol_effect_m[3] ;
  real recipe_effect_b[3] ;
  real temp_effect_b[3] ;
  real ol_effect_b[3] ;
  real<lower=0> sigma; // measurement error
  real<lower=0> sigma_b; // measurement error
 }
transformed parameters {
  vector[nrows] coli_hat; // the *real* pathogen load, without measurement error
  vector[nrows] cocci_hat; // the *real* pathogen load, without measurement error
  for (i in 1:nrows){
    coli_hat[i] = ((temp_effect_m[Temp[i]] + recipe_effect_m[Recipe[i]] +  ol_effect_m[OL[i]])  * Time[i]) + (recipe_effect_b[Recipe[i]] + temp_effect_b[Temp[i]] + ol_effect_b[OL[i]]);
    cocci_hat[i] = ((temp_effect_m[Temp[i]] + recipe_effect_m[Recipe[i]] + ol_effect_m[OL[i]])  * Time[i]) + (recipe_effect_b[Recipe[i]] + temp_effect_b[Temp[i]] + ol_effect_b[OL[i]]);
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
  sigma_b ~ cauchy(0, 1);
  recipe_effect_b ~ cauchy(0, 1);
  recipe_effect_m ~ cauchy(0, 1);
  ol_effect_m ~ cauchy(0, 1);
  ol_effect_b ~ cauchy(0, 1);
  temp_effect_m ~ cauchy(0, 1);
  temp_effect_b ~ cauchy(0, 1);
  //sigma_co ~ cauchy(sigma, eta);
  //sigma_ec ~ cauchy(sigma, eta);
  //sigma_en ~ cauchy(sigma, eta);
  //for(i in 1:nrows){
  Ecoli ~ normal(coli_hat, sigma);
  Enterococci ~ normal(cocci_hat, sigma_b);
  //}
}

