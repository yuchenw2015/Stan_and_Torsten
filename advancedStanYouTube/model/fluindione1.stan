functions{

  real[] fluindioneODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy){
    real k10;
    real k12;
    real k21;
    real CL;
    real Q;
    real V1;
    real V2;
    real ka;
    real Emax;
    real EC50;
    real inr0;
    real kin;
    real kout;
    real gamma;
    real dxdt[4];
    real conc;
    real EDrug;

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
    Emax = parms[6];
    EC50 = parms[7];
    inr0 = parms[8];
    kout = parms[9];
    gamma = parms[10];

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;
    kin = kout / inr0;

    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = fmax(0, x[2]/V1);
    EDrug = Emax * conc^gamma / (EC50^gamma + conc^gamma);
    // x[4] = 1/INR - 1/INR(0)
    dxdt[4] = kin * (1 - EDrug) - kout * (x[4] + kin / kout);

    return dxdt;
    }
  
}

data{
  int<lower = 1> nt;
  real<lower = 0> dose;
  real<lower = 0> time[nt];
  vector<lower = 0>[nt] cObs;
  vector<lower = 0>[nt] inrObs;
}

transformed data{
  int<lower = 1> nCmt = 4;
  int idummy[0];
  real rdummy[0];
  real<lower = 0> init[nCmt] = {dose, 0, 0, 0};
  vector[nt] logcObs = log(cObs);
  vector[nt] loginrObs = log(inrObs);
  
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0, upper = 1> Emax;
  real<lower = 0> EC50;
  real<lower = 0> inr0;
  real<lower = 0> kout;
  real<lower = 0, upper = 10> gamma;
  real<lower = 0> sigmaPK;
  real<lower = 0> sigmaPD;
}

transformed parameters{
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nt] inrHat;
  real x[nt, nCmt];
  real<lower = 0> parms[10];

  parms = {CL, Q, V1, V2, ka, Emax, EC50, inr0, kout, gamma};

  x[1,] = init;
  x[2:nt,] =  integrate_ode_bdf(fluindioneODE, init, 
				 time[1], time[2:nt],
				 parms, rdummy, idummy,
				 1E-6, 1E-6, 1E8);

  cHat = to_vector(x[, 2]) / V1;
  inrHat = 1 ./ (to_vector(x[, 4]) + 1 / inr0);

}

model{
  CL ~ normal(0, 0.5); 
  Q ~ normal(0, 0.5);
  V1 ~ normal(0, 10);
  V2 ~ normal(0, 10);
  ka ~ normal(0, 5);
  Emax ~ beta(3, 1);
  EC50 ~ normal(0, 10);
  inr0 ~ normal(0, 2);
  kout ~ normal(0, 1);
  gamma ~ normal(0, 5);

  sigmaPK ~ cauchy(0, 1);
  sigmaPD ~ cauchy(0, 1);

  logcObs[2:nt] ~ normal(log(cHat[2:nt]), sigmaPK); // observed data likelihood
  loginrObs ~ normal(log(inrHat), sigmaPD);
}

generated quantities{
  real cObsPred[nt];
  real inrObsPred[nt];

  for(i in 1:nt){
    if(time[i] == 0){
      cObsPred[i] = 0;
    }else{
      cObsPred[i] = exp(normal_rng(log(cHat[i]), sigmaPK));
    }
    inrObsPred[i] = exp(normal_rng(log(inrHat[i]), sigmaPD));
  }

}
