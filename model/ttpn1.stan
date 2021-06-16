////////////////////////////////////////////////////////////////////////
//// Adapted by Yuchen Wang                     
//// Scripts adapted due to the updates of Torsten built-in functions                      
//// Function names and the related matrix/vector dimensions apdated
//// R scripts adapted to call these Torsten functions, see .R files
//// Date: June/15/2021
//// email: yuchenw2015@gmail.com
//// Based on the PKPD Stan course by Bill Gillespie
//// Link of the original materials: 
//// https://www.metrumrg.com/course/advanced-use-stan-rstan-torsten-pharmacometric-applications/
///////////////////////////////////////////////////////////////////////
functions{

    real[] oneCptPNODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy){
    real dxdt[3];
    real CL = parms[1];
    real V = parms[2];
    real ke0 = parms[3];
    real alpha = parms[4];
    real beta = parms[5];
    real Edrug = fmax(machine_precision(), alpha * x[2]);
    real hazard;

    dxdt[1] = -(CL / V) * x[1];
    dxdt[2] = ke0 * (x[1] / V - x[2]);
    if(t <= 0 || beta <= 0){
      hazard = 0;
    }else{
      hazard = beta * Edrug^beta * t^(beta - 1);
    }
    dxdt[3] = hazard;
    
    return dxdt;
  }

}

data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nPKObs;
  int<lower = 1> nPNObs;
  int<lower = 1> nPNCens;
  int<lower = 1> iPKObs[nPKObs];
  int<lower = 1> iPNObs[nPNObs];
  int<lower = 1> iPNCens[nPNCens];
  real<lower = 0> amt[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  real<lower = 0> time[nt];
  vector<lower = 0>[nPKObs] cObs;
}

transformed data{
  int<lower = 0> ss[nt] = rep_array(0, nt);
  vector[nPKObs] logCObs = log(cObs);
  int<lower = 1> nRandom = 2;
  int<lower = 1> nCmt = 3;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);

  int ntPred = 253;
  real dt = 12;
  real tPred[ntPred];
  real ratePred[ntPred] = rep_array(0.0, ntPred);
  real iiPred[ntPred];
  int addlPred[ntPred];
  int cmtPred[ntPred] = rep_array(1, ntPred);
  int evidPred[ntPred];   
  int ssPred[ntPred] = rep_array(0, ntPred);

  iiPred[1] = 21 * 24;
  iiPred[2:ntPred] = rep_array(0.0, ntPred - 1);
  addlPred[1] = 5;
  addlPred[2:ntPred] = rep_array(0, ntPred - 1);
  evidPred[1] = 1;
  evidPred[2:ntPred] = rep_array(0, ntPred - 1);
  for(i in 1:ntPred) tPred[i] = dt * (i - 1);
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> VHat;
  real<lower = 0> ke0;
  real<lower = 0> alpha;
  real<lower = 0> beta;
  vector<lower = 0>[nRandom] omega;
  cholesky_factor_corr[nRandom] L;
  real<lower = 0> sigma;
  matrix[nRandom, nId] eta;

}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat = to_vector({CLHat, VHat});
  real<lower = 0> CL[nId];
  real<lower = 0> V[nId];
  matrix<lower = 0>[nId, nRandom] theta;

  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nPKObs] cHatObs;
  vector<lower = 0>[nPNObs] survObs;
  vector<lower = 0>[nPNObs] EdrugObs;
  vector<lower = 0>[nPNObs] hazardObs;
  vector<lower = 0>[nPNCens] survCens;
  matrix<lower = -10 * machine_precision()>[nt, 3] x;
  real<lower = 0> parms[5];

  theta = (rep_matrix(thetaHat, nId) .* 
	   exp(diag_pre_multiply(omega, L * eta)))';

  for(j in 1:nId){
    CL[j] = theta[j, 1];
    V[j] = theta[j, 2];

    parms = {CL[j], V[j], ke0, alpha, beta};
//    print(parms)

    // x[start[j]:end[j],] = generalOdeModel_rk45(oneCptPNODE, 3,
    x[start[j]:end[j],] = (pmx_solve_rk45(oneCptPNODE, 3,
			     time[start[j]:end[j]], 
			     amt[start[j]:end[j]],
			     rate[start[j]:end[j]],
			     ii[start[j]:end[j]],
			     evid[start[j]:end[j]],
			     cmt[start[j]:end[j]],
			     addl[start[j]:end[j]],
			     ss[start[j]:end[j]],
			     parms, F, tLag,
			     1e-6, 1e-6, 1e8))'; // adapt function name and dimension
			     
		for(i in start[j]:end[j]){
		  for(k in 1:nCmt){
		    x[i, k] = x[i, k] < 0 ? machine_precision() : x[i, k];
		  }
		}

    cHat[start[j]:end[j]] = x[start[j]:end[j], 1] / V[j];
  }

  cHatObs = cHat[iPKObs]; // predictions for observed data records
  survObs = exp(-x[iPNObs, 3]);
  EdrugObs = alpha * x[iPNObs, 2];
  for(i in 1:nPNObs)
    hazardObs[i] = beta * EdrugObs[i]^beta * 
      time[iPNObs[i]]^(beta - 1);
  survCens = exp(-x[iPNCens, 3]);
  //  print("cHatObs:", cHatObs);
  //  print("hazardObs .* survObs:", hazardObs .* survObs);
  //  print("survCens:", survCens);
}

model{
  CLHat ~ normal(0, 0.005);
  VHat ~ normal(0, 0.25);
  ke0 ~ normal(0, 0.0005);
  alpha ~ normal(0, 0.000003);
  beta ~ normal(0, 1.5);
  sigma ~ cauchy(0, 1);
  omega ~ cauchy(0, 1);
  L ~ lkj_corr_cholesky(1);

  to_vector(eta) ~ normal(0, 1);
  
  logCObs ~ normal(log(cHatObs), sigma); // observed concentration log likelihood
  target += log(hazardObs .* survObs); // observed PN event log likelihood
  target += log(survCens); // censored PN event log likelihood
}

generated quantities{
  real<lower = 0> CLPred[nId];
  real<lower = 0> VPred[nId];
  matrix<lower = 0>[nId, nRandom] thetaPred;
  vector<lower = 0>[nt] cHatPred;
  vector[nt] cObsCond;
  vector[nt] cObsPred;
  matrix<lower = -10 * machine_precision()>[nt, 3] xPred;
  real<lower = 0> parmsPred[5];
  corr_matrix[nRandom] rho;
  matrix[nRandom, nId] etaPred;
  vector<lower = 0, upper = 1>[ntPred] cdfPred[nId];
  real<lower = 0> amtPred[ntPred] = rep_array(0.0, ntPred);

  rho = L * L';
  for(j in 1:nId) 
    for(i in 1:nRandom)
      etaPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nId) .* 
	   exp(diag_pre_multiply(omega, L * etaPred)))';

  for(j in 1:nId){
    CLPred[j] = thetaPred[j, 1];
    VPred[j] = thetaPred[j, 2];

    parmsPred = {CLPred[j], VPred[j], ke0, alpha, beta};

    // xPred[start[j]:end[j],] = generalOdeModel_rk45(oneCptPNODE, 3,
    xPred[start[j]:end[j],] = (pmx_solve_rk45(oneCptPNODE, 3,
			     time[start[j]:end[j]], 
			     amt[start[j]:end[j]],
			     rate[start[j]:end[j]],
			     ii[start[j]:end[j]],
			     evid[start[j]:end[j]],
			     cmt[start[j]:end[j]],
			     addl[start[j]:end[j]],
			     ss[start[j]:end[j]],
			     parmsPred, F, tLag,
			     1e-6, 1e-6, 1e8))'; // adapt function name and dimension
			     
		for(i in start[j]:end[j]){
		  for(k in 1:nCmt){
		    xPred[i, k] = xPred[i, k] < 0 ? machine_precision() : xPred[i, k];
		  }
		}

    cHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 1] / VPred[j];
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

  // Simulate CDF for time to PN

  for(j in 1:nId){
    parmsPred = {CLPred[j], VPred[j], ke0, alpha, beta};
    amtPred[1] = amt[start[j]];

    //cdfPred[j] =  1 - exp(-generalOdeModel_rk45(oneCptPNODE, 3,
    cdfPred[j] =  1 - exp(-(pmx_solve_rk45(oneCptPNODE, 3,
						tPred, 
						amtPred,
						ratePred,
						iiPred,
						evidPred,
						cmtPred,
						addlPred,
						ssPred,
						parmsPred, F, tLag,
						1e-6, 1e-6, 1e8))'[,3]); // adapt function name and dimension

  }

}
