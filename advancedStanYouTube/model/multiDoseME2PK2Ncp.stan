data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 0> nBql;
  int<lower = 1> iObs[nObs];
  int<lower = 1> iBql[nBql]; ## indices of BQL data records
  real<lower = 0> loq; ## limit of quantitation
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> weight[nSubjects];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int<lower = 1> nRandom = 5;
  int<lower = 1> nCmt = 3;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
//  corr_matrix[nRandom] rho;
  cholesky_factor_corr[nRandom] L;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
//  vector[nRandom] logtheta[nSubjects];
  matrix[nRandom, nSubjects] eta;
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
//  cov_matrix[nRandom] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> Q[nSubjects];
  real<lower = 0> V1[nSubjects];
  real<lower = 0> V2[nSubjects];
  real<lower = 0> ka[nSubjects];
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  vector<lower = 0>[nBql] cHatBql;
  matrix[nt, nCmt] x;
  matrix<lower = 0>[nSubjects, nRandom] theta;
  real<lower = 0> parms[5];

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

//  Omega = quad_form_diag(rho, omega); ## diag_matrix(omega) * rho * diag_matrix(omega)
  theta = (rep_matrix(thetaHat, nSubjects) .* 
          exp(diag_pre_multiply(omega, L * eta)))';

  for(j in 1:nSubjects){
    CL[j] = theta[j, 1] * (weight[j] / 70)^0.75;
    Q[j] = theta[j, 2] * (weight[j] / 70)^0.75;
    V1[j] = theta[j, 3] * weight[j] / 70;
    V2[j] = theta[j, 4] * weight[j] / 70;
    ka[j] = theta[j, 5];

    parms[1] = CL[j];
    parms[2] = Q[j];
    parms[3] = V1[j];
    parms[4] = V2[j];
    parms[5] = ka[j];

    x[start[j]:end[j],] = PKModelTwoCpt(time[start[j]:end[j]], 
					amt[start[j]:end[j]],
					rate[start[j]:end[j]],
					ii[start[j]:end[j]],
					evid[start[j]:end[j]],
					cmt[start[j]:end[j]],
					addl[start[j]:end[j]],
					ss[start[j]:end[j]],
					parms, F, tLag);

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] ./ V1[j];
  }

  cHatObs = cHat[iObs]; ## predictions for observed data records
  cHatBql = cHat[iBql]; ## predictions for BQL data records

}

model{
  CLHat ~ normal(0, 25);
  QHat ~ normal(0, 25);
  V1Hat ~ normal(0, 100);
  V2Hat ~ normal(0, 300);
  kaHat ~ normal(0, 5);
  omega ~ cauchy(0, 1);
  //  rho ~ lkj_corr(1); 
  L ~ lkj_corr_cholesky(1);
  sigma ~ cauchy(0, 1);

  ## Inter-individual variability
//  logtheta ~ multi_normal(log(thetaHat), Omega);
  to_vector(eta) ~ normal(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); ## observed data likelihood
  target += normal_lcdf(log(loq) | log(cHatBql), sigma); ## BQL data likelihood
}

generated quantities{
//  vector[nRandom] logthetaPred[nSubjects];
  vector<lower = 0>[nt] cHatPred;
  vector[nt] cObsCond;
  vector[nt] cObsPred;
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> QPred[nSubjects];
  real<lower = 0> V1Pred[nSubjects];
  real<lower = 0> V2Pred[nSubjects];
  real<lower = 0> kaPred[nSubjects];
  matrix[nt, nCmt] xPred;
  corr_matrix[nRandom] rho;
  matrix[nRandom, nSubjects] etaPred;
  matrix<lower = 0>[nSubjects, nRandom] thetaPred;
  real<lower = 0> parmsPred[5];

  rho = L * L';
  for(j in 1:nSubjects) 
    for(i in 1:nRandom)
      etaPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nSubjects) .* 
              exp(diag_pre_multiply(omega, L * etaPred)))';

  for(j in 1:nSubjects){
//    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    CLPred[j] = thetaPred[j, 1] * (weight[j] / 70)^0.75;
    QPred[j] = thetaPred[j, 2] * (weight[j] / 70)^0.75;
    V1Pred[j] = thetaPred[j, 3] * weight[j] / 70;
    V2Pred[j] = thetaPred[j, 4] * weight[j] / 70;
    kaPred[j] = thetaPred[j, 5];

    parmsPred[1] = CLPred[j];
    parmsPred[2] = QPred[j];
    parmsPred[3] = V1Pred[j];
    parmsPred[4] = V2Pred[j];
    parmsPred[5] = kaPred[j];

    xPred[start[j]:end[j],] = PKModelTwoCpt(time[start[j]:end[j]], 
					    amt[start[j]:end[j]],
					    rate[start[j]:end[j]],
					    ii[start[j]:end[j]],
					    evid[start[j]:end[j]],
					    cmt[start[j]:end[j]],
					    addl[start[j]:end[j]],
					    ss[start[j]:end[j]],
					    parmsPred, F, tLag);

    cHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 2] ./ V1Pred[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
    }else{
      cObsCond[i] = exp(normal_rng(log(cHat[i]), sigma));
      cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
    }
  }
}
