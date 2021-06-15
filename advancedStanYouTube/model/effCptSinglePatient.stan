data{
  int<lower = 1> nObs;
  real<lower = 0> dose;
  real<lower = 0> time[nObs];
  real<lower = 0> cObs[nObs];
  real respObs[nObs];
}

transformed data{
  real logCObs[nObs - 1] = log(cObs[2:nObs]);
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> ke0;
  real<lower = 0> EC50;
  real<lower = 0> sigma;
  real<lower = 0> sigmaResp;
}

transformed parameters{
  vector<lower = 0>[nObs] cHat;
  vector<lower = 0>[nObs] respHat;
  vector<lower = 0>[nObs] ceHat;
  matrix[nObs, 4] x;
  matrix[4, 4] K;

  K = rep_matrix(0, 4, 4);
  
  K[1, 1] = -ka;
  K[2, 1] = ka;
  K[2, 2] = -(CL + Q) / V1;
  K[2, 3] = Q / V2;
  K[3, 2] = Q / V1;
  K[3, 3] = -Q / V2;
  K[4, 2] = ke0;
  K[4, 4] = -ke0;

  for(i in 1:nObs)
    x[i, ] = (matrix_exp(time[i] * K) * to_vector({dose, 0, 0, 0}))';
                  
  cHat = 1000 * x[ ,2] ./ V1;
  ceHat = 1000 * x[ ,4] ./ V1;
  respHat = 100 * ceHat ./ (EC50 + ceHat);

}

model{
    CL ~ normal(0, 20);
    Q ~ normal(0, 40);
    V1 ~ normal(0, 150);
    V2 ~ normal(0, 150);
    ka ~ normal(0, 5);
    ke0 ~ normal(0, 2);
    EC50 ~ normal(0, 200);
    sigma ~ cauchy(0, 2);
    sigmaResp ~ cauchy(0, 5);

    logCObs ~ normal(log(cHat[2:nObs]), sigma); 
    respObs ~ normal(respHat, sigmaResp); 
}

generated quantities{
  real cObsPred[nObs];
  real respObsPred[nObs];
  
  for(i in 1:nObs) {
    if(time[i] == 0){
      cObsPred[i] = 0;
    }
    else{
      cObsPred[i] = exp(normal_rng(log(cHat[i]), sigma));
    }
    respObsPred[i] = normal_rng(respHat[i], sigmaResp);
  }

}

