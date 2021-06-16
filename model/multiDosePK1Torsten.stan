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
data{
  // General data items
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
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
  
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  
  // Individual-level model parameters directly sampled from the IIV
  // distribution
  vector[nRandom] logtheta[nSubjects];
}

transformed parameters{
  // Vector of PK parameter typical values -- only those with IIV
  vector<lower = 0>[nRandom] thetaHat = to_vector({CLHat, VHat, kaHat});

  // Individual-level model parameters with recognizable names, e.g.,
  real<lower = 0> CL[nSubjects];
  real<lower = 0> V[nSubjects];
  real<lower = 0> ka[nSubjects];
  
  // Covariance matrix
  cov_matrix[nRandom] Omega;

  // Predicted concentrations (without residual variation)
  vector<lower = 0>[nt] cHat; // All events
  vector<lower = 0>[nObs] cHatObs; // Observation events

  // Amounts in each compartment at each event
  matrix[nt, nCmt] x;

  // Array used to pass parameters to the Torsten function
  real<lower = 0> parms[nParms];

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)

  for(j in 1:nSubjects){
    
    // Calculation of individual parameter values given logtheta and covariates, e.g.
    CL[j] = exp(logtheta[j, 1]);
    V[j] = exp(logtheta[j, 2]);
    ka[j] = exp(logtheta[j, 3]);

    // Pack individual PK parameters into parms array, e.g.
    
    parms = {CL[j], V[j], ka[j]};

    // x[start[j]:end[j],] = PKModelOneCpt(time[start[j]:end[j]], 
    x[start[j]:end[j],] = (pmx_solve_onecpt(time[start[j]:end[j]], 
					amt[start[j]:end[j]],
					rate[start[j]:end[j]],
					ii[start[j]:end[j]],
					evid[start[j]:end[j]],
					cmt[start[j]:end[j]],
					addl[start[j]:end[j]],
					ss[start[j]:end[j]],
					parms, F, tLag))'; // adapt function name and dimension

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] ./ V[j];
  }

  cHatObs = cHat[iObs]; // predictions for observed data records
}

model{
  // Priors for PK model-specific parameters
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

  logCObs ~ normal(log(cHatObs), sigma); // observed data likelihood`
}

generated quantities{
  vector[nRandom] logthetaPred[nSubjects];
  vector<lower = 0>[nt] cHatPred;
  vector[nt] cObsCond;
  vector[nt] cObsPred;

  // Individual-level model parameters with recognizable names, e.g.,
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> VPred[nSubjects];
  real<lower = 0> kaPred[nSubjects];

  matrix[nt, nCmt] xPred;
  real<lower = 0> parmsPred[nParms];

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);

    // Calculation of individual parameter values given logtheta and covariates, e.g.
    CLPred[j] = exp(logthetaPred[j, 1]);
    VPred[j] = exp(logthetaPred[j, 2]);
    kaPred[j] = exp(logthetaPred[j, 3]);

    // Pack individual PK parameters into parms array, e.g.
    
    parmsPred = {CLPred[j], VPred[j], kaPred[j]};

    // xPred[start[j]:end[j],] = PKModelOneCpt(time[start[j]:end[j]], 
    xPred[start[j]:end[j],] = (pmx_solve_onecpt(time[start[j]:end[j]], 
					    amt[start[j]:end[j]],
					    rate[start[j]:end[j]],
					    ii[start[j]:end[j]],
					    evid[start[j]:end[j]],
					    cmt[start[j]:end[j]],
					    addl[start[j]:end[j]],
					    ss[start[j]:end[j]],
					    parmsPred, F, tLag))'; // adapt function name and dimension

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 2] ./ VPred[j];
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
