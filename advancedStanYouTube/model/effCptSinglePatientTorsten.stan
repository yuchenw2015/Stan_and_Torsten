data{
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> time[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> amt[nt];
  real<lower = 0> rate[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> addl[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> ss[nt];
  real<lower = 0> cObs[nObs];
  real respObs[nt];
}

transformed data{
  real logCObs[nObs] = log(cObs);
  int<lower = 1> nCmt = 4;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  //  real<lower = 0> ka; // ka unconstrained
  real<lower = (CL / V1 + Q / V1 + Q / V2 +
		sqrt((CL / V1 + Q / V1 + Q / V2)^2 -
		     4 * CL / V1 * Q / V2)) / 2> ka; // ka > lambda_1
  real<lower = 0> ke0;
  real<lower = 0> EC50;
  real<lower = 0> sigma;
  real<lower = 0> sigmaResp;
}

transformed parameters{
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nt] respHat;
  vector<lower = 0>[nt] ceHat;
  matrix[nt, nCmt] x;
  matrix[nCmt, nCmt] K;

  K = rep_matrix(0, nCmt, nCmt);
  
  K[1, 1] = -ka;
  K[2, 1] = ka;
  K[2, 2] = -(CL + Q) / V1;
  K[2, 3] = Q / V2;
  K[3, 2] = Q / V1;
  K[3, 3] = -Q / V2;
  K[4, 2] = ke0;
  K[4, 4] = -ke0;

  x = linOdeModel(time, amt, rate, ii, evid, cmt, addl, ss,
		  K, F, tLag);
                  
  cHat = x[ ,2] ./ V1;
  ceHat = x[ ,4] ./ V1;
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

    logCObs ~ normal(log(cHat[iObs]), sigma); 
    respObs ~ normal(respHat, sigmaResp); 
}

generated quantities{
  real cObsPred[nt];
  real respObsPred[nt];
  
  for(i in 1:nt) {
    if(time[i] == 0){
      cObsPred[i] = 0;
    }
    else{
      cObsPred[i] = exp(normal_rng(log(cHat[i]), sigma));
    }
    respObsPred[i] = normal_rng(respHat[i], sigmaResp);
  }

}

