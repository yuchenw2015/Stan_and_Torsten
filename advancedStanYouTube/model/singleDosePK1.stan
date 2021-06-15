functions{

  real twoCptModel1(real time, real dose, real CL, real Q, real V1, real V2, real ka){

    real k10;
    real k12;
    real k21;
    real ksum;
    vector[3] a;
    vector[3] alpha;

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;
    ksum = k10 + k12 + k21;
    alpha[1] = (ksum + sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[2] = (ksum - sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[3] = ka;

    a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1]) * (alpha[2] - alpha[1]));
    a[2] = ka * (k21 - alpha[2]) / ((ka - alpha[2]) * (alpha[1] - alpha[2]));
    a[3] = -(a[1] + a[2]);

    return (dose / V1) * sum(a .* exp(-alpha * time));
  }

}

data{
  int nSubjects;
  int nObs;
  real dose[nSubjects];
  real weight[nSubjects];
  int subject[nObs];
  vector[nObs] time;
  vector[nObs] cObs;
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nRandom = 5;
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  vector[5] logtheta[nSubjects];
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[5] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> Q[nSubjects];
  real<lower = 0> V1[nSubjects];
  real<lower = 0> V2[nSubjects];
  real<lower = 0> ka[nSubjects];
  vector<lower = 0>[nObs] cHat;

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)

  for(j in 1:nSubjects){
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
  }

  for(i in 1:nObs)
    cHat[i] = twoCptModel1(time[i], dose[subject[i]], CL[subject[i]],
			   Q[subject[i]], V1[subject[i]], V2[subject[i]],
			   ka[subject[i]]);  
}

model{
    CLHat ~ normal(0, 25);
    QHat ~ normal(0, 50);
    V1Hat ~ normal(0, 100);
    V2Hat ~ normal(0, 200);
    kaHat ~ normal(0, 5);
    omega ~ cauchy(0, 1);
    rho ~ lkj_corr(1); 
    sigma ~ cauchy(0, 1);

    // Inter-individual variability
    logtheta ~ multi_normal(log(thetaHat), Omega);
            
    logCObs ~ normal(log(cHat), sigma);
}

generated quantities{
  vector[nRandom] logthetaPred[nSubjects];
  real<lower = 0> cHatPred[nObs];
  real<lower = 0> cObsCond[nObs];
  real<lower = 0> cObsPred[nObs];
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> QPred[nSubjects];
  real<lower = 0> V1Pred[nSubjects];
  real<lower = 0> V2Pred[nSubjects];
  real<lower = 0> kaPred[nSubjects];

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    CLPred[j] = exp(logthetaPred[j, 1]) * (weight[j] / 70)^0.75;
    QPred[j] = exp(logthetaPred[j, 2]) * (weight[j] / 70)^0.75;
    V1Pred[j] = exp(logthetaPred[j, 3]) * weight[j] / 70;
    V2Pred[j] = exp(logthetaPred[j, 4]) * weight[j] / 70;
    kaPred[j] = exp(logthetaPred[j, 5]);
  }

  for(i in 1:nObs){
    cObsCond[i] = exp(normal_rng(log(cHat[i]), sigma));

    cHatPred[i] = twoCptModel1(time[i], dose[subject[i]], CLPred[subject[i]],
			    QPred[subject[i]], V1Pred[subject[i]], V2Pred[subject[i]],
			    kaPred[subject[i]]);
    cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
  }

}

