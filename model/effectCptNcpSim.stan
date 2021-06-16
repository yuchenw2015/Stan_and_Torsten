////////////////////////////////////////////////////////////////////////
//// Adapted by Yuchen Wang                     
//// Scripts adapted due to the updates of Torsten built-in functions                      
//// Function names and the related matrix/vector dimensions apdated
//// R scripts adapted to call these Torsten functions via cmdstan, see .R files
//// Date: June/15/2021
//// email: yuchenw2015@gmail.com
//// Based on the PKPD Stan course by Bill Gillespie
//// Link of the original materials: 
//// https://www.metrumrg.com/course/advanced-use-stan-rstan-torsten-
//// pharmacometric-applications/
///////////////////////////////////////////////////////////////////////

data{
  int<lower = 1> nId;
  int<lower = 1> nt;
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
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  real<lower = 0> ke0Hat;
  real<lower = 0> EmaxHat;
  real<lower = 0> EC50;
  real<lower = 0> E0Hat;
  real<lower = 0> gamma;
  int<lower = 0> nRandom;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  real<lower = 0> sigmaPD;
}

transformed data{
  int<lower = 1> nCmt = 4;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
  vector<lower = 0>[nRandom] thetaHat =
    to_vector({CLHat, QHat, V1Hat, V2Hat, kaHat, ke0Hat, EmaxHat, E0Hat});
  cov_matrix[nRandom] Omega;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
}

parameters{
}

transformed parameters{
}

model{
}

generated quantities{
  vector[nRandom] logtheta[nId];
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
  real cObs[nt];
  real respObs[nt];
  matrix[nCmt, nCmt] K;
  matrix[nt, nCmt] x;

  for(j in 1:nId){
    logtheta[j] = multi_normal_rng(log(thetaHat), Omega);
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
    ke0[j] = exp(logtheta[j, 6]);
    Emax[j] = exp(logtheta[j, 7]);
    E0[j] = exp(logtheta[j, 8]);

    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -(CL[j] + Q[j]) / V1[j];
    K[2, 3] = Q[j] / V2[j];
    K[3, 2] = Q[j] / V1[j];
    K[3, 3] = -Q[j] / V2[j];
    K[4, 2] = ke0[j];
    K[4, 4] = -ke0[j];

    // x[start[j]:end[j],] = linOdeModel(time[start[j]:end[j]],  
    x[start[j]:end[j],] = (pmx_solve_linode(time[start[j]:end[j]], // adapt function names  
				      amt[start[j]:end[j]],
				      rate[start[j]:end[j]],
				      ii[start[j]:end[j]],
				      evid[start[j]:end[j]],
				      cmt[start[j]:end[j]],
				      addl[start[j]:end[j]],
				      ss[start[j]:end[j]],
				      //K, F, tLag);
				      K, F, tLag))'; //adapt dimension

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] / V1[j];
    ceHat[start[j]:end[j]] = x[start[j]:end[j], 4] / V1[j];
    for(i in start[j]:end[j])
      respHat[i] = E0[j] + Emax[j] * ceHat[i]^gamma / (EC50^gamma + ceHat[i]^gamma);
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObs[i] = 0;
    }else{
      cObs[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])), sigma));
    }
    respObs[i] = exp(normal_rng(log(fmax(machine_precision(), respHat[i])), sigmaPD));
  }

}
