functions{
  real normalTrunc_lpdf(real y, real mu, real sigma, real lower, real upper){
    real result;
    
    if(y >= lower && y <= upper){
      result = normal_lpdf(y | mu, sigma) -
	log(normal_cdf(upper, mu, sigma) - normal_cdf(lower, mu, sigma));
    }else{
      result = negative_infinity();
    }
    return result;
  }
}

data{
  int<lower = 1> nObs;
  real<lower = 0> c24[nObs];
  real fxa24[nObs];
}

parameters{
  real<lower = 0, upper = 100> emax;
  real<lower = 0> ec50;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
 }

transformed parameters{
  real fxa24Hat[nObs];

  for(i in 1:nObs){
    fxa24Hat[i] = emax * c24[i]^gamma / (ec50^gamma + c24[i]^gamma);
  }
}

model{
  emax ~ normalTrunc(50, 25, 0, 100);
  ec50 ~ normal(0, 250);
  gamma ~ normal(0, 5);
  sigma ~ cauchy(0, 10);

  fxa24 ~ normal(fxa24Hat, sigma);
}

generated quantities{
  real fxa24Pred[nObs];
  
  for(i in 1:nObs){
    fxa24Pred[i] = normal_rng(fxa24Hat[i], sigma);
  }
}

