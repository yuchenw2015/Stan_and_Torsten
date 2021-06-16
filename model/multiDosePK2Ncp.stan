////////////////////////////////////////////////////////////////////////
//// Adapted by Yuchen Wang                     
//// Scripts adapted due to the updates of Torsten built-in functions                      
//// Function names and the related matrix/vector dimensions apdated
//// R scripts adapted to call these Torsten functions, see .R files
//// Date: June/15/2021
//// email: yuchenw2015@gmail.com
//// Based on the PKPD Stan course by Bill Gillespie
//// Link of the original materials: 
//// https://www.metrumrg.com/course/advanced-use-stan-rstan-torsten-
//// pharmacometric-applications/
///////////////////////////////////////////////////////////////////////
data{
  // General data items
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  int<lower = 0> nBql;
  int<lower = 1> iBql[nBql]; // indices of BQL data records
  real<lower = 0> loq; // limit of quantitation
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;

  // Data items for posterior predictions
  int<lower = 1> ntPred;
  real<lower = 0> amtPred[ntPred];
  int<lower = 1> cmtPred[ntPred];
  int<lower = 0> evidPred[ntPred];
  real<lower = 0> ratePred[ntPred];
  real<lower = 0> iiPred[ntPred];
  int<lower = 0> addlPred[ntPred];
  int<lower = 0> ssPred[ntPred];
  int<lower = 1> startPred[nSubjects];
  int<lower = 1> endPred[nSubjects];
  real<lower = 0> tPred[ntPred];
}

transformed data{
  // Transformed dependent variable, e.g.,

  vector[nObs] logCObs = log(cObs);

  // Integers required to specify dimensions
  int<lower = 1> nRandom = 3; // Number of random effects
  int<lower = 1> nCmt = 2; // Number of model compartments
  int<lower = 1> nParms = 3; // Number of parameters passed to Torsten function

  // Fixed value parameters, e.g.,
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  // Population-level model parameters
  // These are the parameters for which you specify prior distributions
  // and initial estimates, e.g., 
  real<lower = 0> CLHat;
  real<lower = 0> VHat;
  real<lower = 0> kaHat;
  
  //  corr_matrix[nRandom] rho;
  cholesky_factor_corr[nRandom] L;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  
  // Individual-level model parameters directly sampled from the IIV
  // distribution
  //  vector[nRandom] logtheta[nSubjects];
  matrix[nRandom, nSubjects] eta;
}

transformed parameters{
  // Vector of PK parameter typical values -- only those with IIV
  vector<lower = 0>[nRandom] thetaHat = to_vector({CLHat, VHat, kaHat});

  // Matrix of individual-level model parameters
  matrix<lower = 0>[nSubjects, nRandom] theta;

  // Individual-level model parameters with recognizable names, e.g.,
  real<lower = 0> CL[nSubjects];
  real<lower = 0> V[nSubjects];
  real<lower = 0> ka[nSubjects];
  
  // Covariance matrix
  //  cov_matrix[nRandom] Omega;

  // Predicted concentrations (without residual variation)
  vector<lower = 0>[nt] cHat; // All events
  vector<lower = 0>[nObs] cHatObs; // Observation events
  vector<lower = 0>[nBql] cHatBql; // BQL events

  // Amounts in each compartment at each event
  matrix[nt, nCmt] x;

  // Array used to pass parameters to the Torsten function
  real<lower = 0> parms[nParms];

  //  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  theta = (rep_matrix(thetaHat, nSubjects) .* 
          exp(diag_pre_multiply(omega, L * eta)))';

  for(j in 1:nSubjects){
    
    // Calculation of individual parameter values given logtheta and covariates, e.g.
    CL[j] = theta[j, 1];
    V[j] = theta[j, 2];
    ka[j] = theta[j, 3];

    // Pack individual PK parameters into parms array, e.g.
    
    parms = {CL[j], V[j], ka[j]};

    //x[start[j]:end[j],] = PKModelOneCpt(time[start[j]:end[j]], 
    x[start[j]:end[j],] = (pmx_solve_onecpt(time[start[j]:end[j]], 
					amt[start[j]:end[j]],
					rate[start[j]:end[j]],
					ii[start[j]:end[j]],
					evid[start[j]:end[j]],
					cmt[start[j]:end[j]],
					addl[start[j]:end[j]],
					ss[start[j]:end[j]],
					parms, F, tLag))'; //adapt function name and dimension

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] ./ V[j];
  }

  cHatObs = cHat[iObs]; // predictions for observed data records
  cHatBql = cHat[iBql]; // predictions for BQL data records
}

model{
  // Priors for PK model-specific parameters
  CLHat ~ normal(0, 20);
  VHat ~ normal(0, 100);
  kaHat ~ lognormal(log(0.45), 0.2);
  omega[1] ~ cauchy(0, 2);
  omega[2] ~ cauchy(0, 2);
  omega[3] ~ lognormal(log(0.25), 0.3);
  
  //  rho ~ lkj_corr(1); 
  L ~ lkj_corr_cholesky(1);
  sigma ~ cauchy(0, 2);

  // Inter-individual variability
  //  logtheta ~ multi_normal(log(thetaHat), Omega);
  to_vector(eta) ~ normal(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); // observed data likelihood
  target += normal_lcdf(log(loq) | log(cHatBql), sigma); // BQL data likelihood

}

generated quantities{
  //  vector[nRandom] logthetaPred[nSubjects];
  matrix[nRandom, nSubjects] etaPred;
  matrix<lower = 0>[nSubjects, nRandom] thetaPred;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[ntPred] cHatCond;
  vector<lower = 0>[ntPred] cHatPred;
  vector[ntPred] cObsCond;
  vector[ntPred] cObsPred;

  // Individual-level model parameters with recognizable names, e.g.,
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> VPred[nSubjects];
  real<lower = 0> kaPred[nSubjects];

  matrix[ntPred, nCmt] xCond;
  matrix[ntPred, nCmt] xPred;
  real<lower = 0> parmsPred[nParms];

  rho = L * L';
  for(j in 1:nSubjects) 
    for(i in 1:nRandom)
      etaPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nSubjects) .* 
              exp(diag_pre_multiply(omega, L * etaPred)))';

  for(j in 1:nSubjects){

    // Individual predictions
    
    // Pack individual PK parameters into parms array, e.g.
    parmsPred = {CL[j], V[j], ka[j]};

   // xCond[startPred[j]:endPred[j],] = PKModelOneCpt(tPred[startPred[j]:endPred[j]], 
    xCond[startPred[j]:endPred[j],] = (pmx_solve_onecpt(tPred[startPred[j]:endPred[j]], 
					    amtPred[startPred[j]:endPred[j]],
					    ratePred[startPred[j]:endPred[j]],
					    iiPred[startPred[j]:endPred[j]],
					    evidPred[startPred[j]:endPred[j]],
					    cmtPred[startPred[j]:endPred[j]],
					    addlPred[startPred[j]:endPred[j]],
					    ssPred[startPred[j]:endPred[j]],
					    parmsPred, F, tLag))'; //adapt function name and dimension

    // Calculate target concentration for specified compartment.

    cHatCond[startPred[j]:endPred[j]] = xCond[startPred[j]:endPred[j], 2] ./ V[j];

    // Population predictions

    //    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);

    // Calculation of individual parameter values given logtheta and covariates, e.g.
    CLPred[j] = thetaPred[j, 1];
    VPred[j] = thetaPred[j, 2];
    kaPred[j] = thetaPred[j, 3];

    // Pack individual PK parameters into parms array, e.g.
    parmsPred = {CLPred[j], VPred[j], kaPred[j]};

    //xPred[startPred[j]:endPred[j],] = PKModelOneCpt(tPred[startPred[j]:endPred[j]], 
    xPred[startPred[j]:endPred[j],] = (pmx_solve_onecpt(tPred[startPred[j]:endPred[j]], 
					    amtPred[startPred[j]:endPred[j]],
					    ratePred[startPred[j]:endPred[j]],
					    iiPred[startPred[j]:endPred[j]],
					    evidPred[startPred[j]:endPred[j]],
					    cmtPred[startPred[j]:endPred[j]],
					    addlPred[startPred[j]:endPred[j]],
					    ssPred[startPred[j]:endPred[j]],
					    parmsPred, F, tLag))'; //adapt function name and dimension

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHatPred[startPred[j]:endPred[j]] = xPred[startPred[j]:endPred[j], 2] ./ VPred[j];
  }

  for(i in 1:ntPred){
    if(tPred[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
    }else{
      cObsCond[i] = exp(normal_rng(log(cHatCond[i]), sigma));
      cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
    }
  }
}
