functions{

  vector oneCptModel1(real dt, vector init, real amt, int cmt, int evid,
		    real CL, real V, real ka){
    vector[2] x;
    vector[2] a;
    vector[2] alpha;

    alpha[1] = CL / V;
    alpha[2] = ka;

    x = rep_vector(0.0, 2);

    if(init[1] != 0.0){
      x[1] = init[1] * exp(-alpha[2] * dt);
      a[1] = alpha[2] / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[2] = init[1] * sum(a .* exp(-alpha * dt));
    }
    
    if(init[2] != 0)
      x[2] = x[2] + init[2] * exp(-alpha[1] * dt);

    if(evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }

  matrix oneCptModel(real[] time, real[] amt, int[] cmt, int[] evid, 
		     real CL, real V, real ka){
    vector[2] init;
    real dt;
    real t0;
    matrix[size(time), 2] result;
    int nt;

    nt = size(time);

    init = rep_vector(0, 2);
    t0 = time[1];
    for(i in 1:nt){
      dt = time[i] - t0;
      init = oneCptModel1(dt, init, amt[i], cmt[i], evid[i],
			   CL, V, ka);
      for(j in 1:2) result[i, j] = init[j];
      t0 = time[i];
    }
    
    return result;
  }
  
}

data{
  int nSubjects;
  int nt;
  int nObs;
  int iObs[nObs];
  real amt[nt];
  int cmt[nt];
  int evid[nt];
  int start[nSubjects];
  int end[nSubjects];
  real time[nt];
  vector[nObs] cObs;
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nti[nSubjects];
  int nRandom = 3;
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> VHat;
  real<lower = 0> kaHat;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  vector[nRandom] logtheta[nSubjects];
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> V[nSubjects];
  real<lower = 0> ka[nSubjects];
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, 2] x;

  thetaHat[1] = CLHat;
  thetaHat[2] = VHat;
  thetaHat[3] = kaHat;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)

  for(j in 1:nSubjects){
    CL[j] = exp(logtheta[j, 1]);
    V[j] = exp(logtheta[j, 2]);
    ka[j] = exp(logtheta[j, 3]);

    x[start[j]:end[j],] = oneCptModel(time[start[j]:end[j]],
				       amt[start[j]:end[j]],
				       cmt[start[j]:end[j]],
				       evid[start[j]:end[j]],
				       CL[j], V[j], ka[j]);
    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] ./ V[j];
  }

  cHatObs = cHat[iObs];

}

model{
    CLHat ~ normal(0, 20);
    VHat ~ normal(0, 100);
    kaHat ~ lognormal(log(0.45), 0.2);
    omega[1] ~ cauchy(0, 2);
    omega[2] ~ cauchy(0, 2);
    omega[3] ~ lognormal(log(0.25), 0.3);
    rho ~ lkj_corr(1); 
    sigma ~ cauchy(0, 2);

    // Inter-individual variability
    logtheta ~ multi_normal(log(thetaHat), Omega);

    logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
  vector[nRandom] logthetaPred[nSubjects];
  vector<lower = 0>[nt] cHatPred;
  real cObsCond[nt];
  real cObsPred[nt];
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> VPred[nSubjects];
  real<lower = 0> kaPred[nSubjects];
  matrix[nt, 2] xPred;

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    CLPred[j] = exp(logthetaPred[j, 1]);
    VPred[j] = exp(logthetaPred[j, 2]);
    kaPred[j] = exp(logthetaPred[j, 3]);

    xPred[start[j]:end[j],] = oneCptModel(time[start[j]:end[j]],
				       amt[start[j]:end[j]],
				       cmt[start[j]:end[j]],
				       evid[start[j]:end[j]],
				       CLPred[j], VPred[j], kaPred[j]);
    cHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 2] ./ V[j];
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

