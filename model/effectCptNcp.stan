data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  real<lower = 0> rate[nt];
  int<lower = 0> ss[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  real<lower = 0> weight[nId];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] respObs;
  real<lower = 0> CLPrior;
  real<lower = 0> QPrior;
  real<lower = 0> V1Prior;
  real<lower = 0> V2Prior;
  real<lower = 0> dkaPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> QPriorCV;
  real<lower = 0> V1PriorCV;
  real<lower = 0> V2PriorCV;
  real<lower = 0> dkaPriorCV;
  real<lower = 0> ke0Prior;
  real<lower = 0> ke0PriorCV;
  real<lower = 0> E0Prior;
  real<lower = 0> E0PriorCV;
  real<lower = 0> EmaxPrior;
  real<lower = 0> EmaxPriorCV;
  real<lower = 0> EC50Prior;
  real<lower = 0> EC50PriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
}

transformed data{
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logRespObs = log(respObs);
  int<lower = 1> nRandom = 8;
  int<lower = 1> nCmt = 4;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> dkaHat;
  real<lower = 0> ke0Hat;
  real<lower = 0> E0Hat;
  real<lower = 0> EmaxHat;
  real<lower = 0> EC50;
  real<lower = 0> gamma;
  vector<lower = 0>[nRandom] omega;
  cholesky_factor_corr[nRandom] L;
  real<lower = 0> sigma;
  real<lower = 0> sigmaPD;
  matrix[nRandom, nId] etaStd;
}

transformed parameters{
  real<lower = 0> kaHat;
  real<lower = 0> k10;
  real<lower = 0> k12;
  real<lower = 0> k21;
  real<lower = 0> ksum;
  real<lower = 0> lambda1;
  vector<lower = 0>[nRandom] thetaHat;
  matrix<lower = 0>[nId, nRandom] theta;
  real<lower = 0> CL[nId];
  real<lower = 0> Q[nId];
  real<lower = 0> V1[nId];
  real<lower = 0> V2[nId];
  real<lower = 0> ka[nId];
  real<lower = 0> ke0[nId];
  real<lower = 0> Emax[nId];
  real<lower = 0> E0[nId];
  vector[nt] cHat;
  vector[nt] ceHat;
  vector[nt] respHat;
  vector[nObsPK] cHatObs;
  vector[nObsPD] respHatObs;
  matrix[nCmt, nCmt] K;
  matrix[nt, nCmt] x;

  k10 = CLHat / V1Hat;
  k12 = QHat / V1Hat;
  k21 = QHat / V2Hat;
  ksum = k10 + k12 + k21;
  lambda1 = (ksum + sqrt(ksum^2 - 4 * k10 * k21)) / 2;
  kaHat = dkaHat + lambda1;

  thetaHat = to_vector({CLHat, QHat, V1Hat, V2Hat, kaHat, ke0Hat, EmaxHat, E0Hat});

  theta = (rep_matrix(thetaHat, nId) .* 
	   exp(diag_pre_multiply(omega, L * etaStd)))';

  for(j in 1:nId){
    CL[j] = theta[j, 1] * (weight[j] / 70)^0.75;
    Q[j] = theta[j, 2] * (weight[j] / 70)^0.75;
    V1[j] = theta[j, 3] * weight[j] / 70;
    V2[j] = theta[j, 4] * weight[j] / 70;
    ka[j] = theta[j, 5];
    ke0[j] = theta[j, 6];
    Emax[j] = theta[j, 7];
    E0[j] = theta[j, 8];

    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -(CL[j] + Q[j]) / V1[j];
    K[2, 3] = Q[j] / V2[j];
    K[3, 2] = Q[j] / V1[j];
    K[3, 3] = -Q[j] / V2[j];
    K[4, 2] = ke0[j];
    K[4, 4] = -ke0[j];

    x[start[j]:end[j],] = linOdeModel(time[start[j]:end[j]], 
				      amt[start[j]:end[j]],
				      rate[start[j]:end[j]],
				      ii[start[j]:end[j]],
				      evid[start[j]:end[j]],
				      cmt[start[j]:end[j]],
				      addl[start[j]:end[j]],
				      ss[start[j]:end[j]],
				      K, F, tLag);

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] / V1[j];
    ceHat[start[j]:end[j]] = x[start[j]:end[j], 4] / V1[j];
    for(i in start[j]:end[j])
      respHat[i] = E0[j] + Emax[j] * ceHat[i]^gamma / (EC50^gamma + ceHat[i]^gamma);
  }

  cHatObs = cHat[iObsPK]; // predictions for observed data records
  respHatObs = respHat[iObsPD]; // predictions for observed data records
}

model{
  CLHat ~ lognormal(log(CLPrior), CLPriorCV);
  QHat ~ lognormal(log(QPrior), QPriorCV);
  V1Hat ~ lognormal(log(V1Prior), V1PriorCV);
  V2Hat ~ lognormal(log(V2Prior), V2PriorCV);
  dkaHat ~ lognormal(log(dkaPrior), dkaPriorCV);
  sigma ~ cauchy(0, 1);
  ke0Hat ~ lognormal(log(ke0Prior), ke0PriorCV);
  E0Hat ~ lognormal(log(E0Prior), E0PriorCV);
  EmaxHat ~ lognormal(log(EmaxPrior), EmaxPriorCV);
  EC50 ~ lognormal(log(EC50Prior), EC50PriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaPD ~ cauchy(0, 1);
  omega ~ cauchy(0, 1);
  L ~ lkj_corr_cholesky(1);

  // Inter-individual variability
  to_vector(etaStd) ~ normal(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); // observed data likelihood
  logRespObs ~ normal(log(respHatObs), sigmaPD);
}

generated quantities{
  real<lower = 0> CLPred[nId];
  real<lower = 0> QPred[nId];
  real<lower = 0> V1Pred[nId];
  real<lower = 0> V2Pred[nId];
  real<lower = 0> kaPred[nId];
  real<lower = 0> ke0Pred[nId];
  real<lower = 0> EmaxPred[nId];
  real<lower = 0> E0Pred[nId];
  vector[nt] cHatPred;
  vector[nt] ceHatPred;
  vector[nt] respHatPred;
  real cObsCond[nt];
  real respObsCond[nt];
  real cObsPred[nt];
  real respObsPred[nt];
  matrix<lower = 0>[nId, nRandom] thetaPred;
  matrix[nt, nCmt] xPred;
  matrix[nCmt, nCmt] KPred;
  corr_matrix[nRandom] rho;
  matrix[nRandom, nId] etaStdPred;

  rho = L * L';
  for(j in 1:nId) 
    for(i in 1:nRandom)
      etaStdPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nId) .* 
	       exp(diag_pre_multiply(omega, L * etaStdPred)))';

  for(j in 1:nId){
    CLPred[j] = thetaPred[j, 1] * (weight[j] / 70)^0.75;
    QPred[j] = thetaPred[j, 2] * (weight[j] / 70)^0.75;
    V1Pred[j] = thetaPred[j, 3] * weight[j] / 70;
    V2Pred[j] = thetaPred[j, 4] * weight[j] / 70;
    kaPred[j] = thetaPred[j, 5];
    ke0Pred[j] = thetaPred[j, 6];
    EmaxPred[j] = thetaPred[j, 7];
    E0Pred[j] = thetaPred[j, 8];

    KPred = rep_matrix(0, nCmt, nCmt);
    
    KPred[1, 1] = -kaPred[j];
    KPred[2, 1] = kaPred[j];
    KPred[2, 2] = -(CLPred[j] + QPred[j]) / V1Pred[j];
    KPred[2, 3] = QPred[j] / V2Pred[j];
    KPred[3, 2] = QPred[j] / V1Pred[j];
    KPred[3, 3] = -QPred[j] / V2Pred[j];
    KPred[4, 2] = ke0Pred[j];
    KPred[4, 4] = -ke0Pred[j];

    xPred[start[j]:end[j],] = linOdeModel(time[start[j]:end[j]], 
					  amt[start[j]:end[j]],
					  rate[start[j]:end[j]],
					  ii[start[j]:end[j]],
					  evid[start[j]:end[j]],
					  cmt[start[j]:end[j]],
					  addl[start[j]:end[j]],
					  ss[start[j]:end[j]],
					  KPred, F, tLag);

    cHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 2] / V1Pred[j];
    ceHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 4] / V1Pred[j];
    for(i in start[j]:end[j])
      respHatPred[i] = E0Pred[j] + EmaxPred[j] * 
        ceHatPred[i]^gamma / (EC50^gamma + ceHatPred[i]^gamma);
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
    }else{
      cObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])),
        sigma));
      cObsPred[i] = exp(normal_rng(log(fmax(machine_precision(),
        cHatPred[i])), sigma));
    }
    respObsCond[i] = exp(normal_rng(log(fmax(machine_precision(),
      respHat[i])), sigmaPD));
    respObsPred[i] = exp(normal_rng(log(fmax(machine_precision(),
      respHatPred[i])), sigmaPD));
  }

}
