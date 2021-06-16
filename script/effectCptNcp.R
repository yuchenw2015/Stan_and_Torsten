#######################################################################
#### Adapted by Yuchen Wang                     
#### Scripts adapted to build and fit Stan model from cmdstans                      
#### Fit object converted from cmdsran fit to stan fit
#### Stan model adapted due to the function update, see .stan files
#### Date: June/15/2021
#### email: yuchenw2015@gmail.com
#### Based on the PKPD Stan course by Bill Gillespie
#### Link of the original materials: 
#### https://www.metrumrg.com/course/advanced-use-stan-rstan-torsten-pharmacometric-applications/
#######################################################################

rm(list = ls())
gc()

modelName <- "effectCptNcp"
simModelName <- "effectCptNcpSim"

## Relative paths assuming the working directory is the Stan_and_Torsten directory

scriptDir <- getwd()
projectDir <- scriptDir
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")

library(rstan)
library(bayesplot)
library(ggplot2)
## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())
library(tidyverse)
library(parallel)
source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "functions.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(11191951) ## not required but assures repeatable results

#Load cmdstanr package and the cmdstan path
#Change the cmdstan path accordingly
library(cmdstanr)
set_cmdstan_path(path = "/pkpd/Torsten/cmdstan")
################################################################################################

### Simulate ME-2 plasma concentrations and ANC values

## Parameter values

ka <- 2.0
CL <- 10 # L/h
V1 <- 35 # L
V2 <- 105 # L
Q <- 15
sigma <- 0.15

ke0 <- 0.5
E0 <- 80
Emax <- 40
EC50 <- 250
gamma <- 2

sigmaPD <- 0.05

omega <- c(0.25, 0.4, 0.25, 0.4, 0.25, 0.25, 0.25, 0.25)
rho <- diag(8)

## Observation and dosing times
doses = c(10, 20, 40, 80)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
xpd <- xpk[xpk %% 1 == 0]
time <- sort(unique(c(xpk, xpd)))

nPerDose <- 5
nId <- nPerDose * length(doses) ## Number of individuals
weight = rnorm.trunc(nId, 70, 15, 50, 100)

## Assemble data set for simulation using Stan
obsData <- data.frame(id = 1:nId,
                      amt = rep(doses * 1000, nPerDose)) %>%
  merge(data.frame(time = time)) %>%
  mutate(rate = 0,
         cmt = 2,
         evid = 0,
         ii = 0,
         addl = 0,
         ss = 0)

doseData <- data.frame(id = 1:nId,
                       amt = rep(doses * 1000, nPerDose)) %>%
    mutate(time = 0,
           rate = 0,
           cmt = 1,
           evid = 1,
           ii = 12,
           addl = 14,
           ss = 0)

allData <- doseData %>%
  bind_rows(obsData) %>%
    arrange(id, time, desc(evid))

nt <- nrow(allData)
start <- (1:nt)[!duplicated(allData$id)]
end <- c(start[-1] - 1, nt)

dataSim <- with(allData,
                list(nId = nId,
                     nt = nt,
                     amt = amt,
                     rate = rate,
                     cmt = cmt,
                     evid = evid,
                     ii = ii,
                     addl = addl,
                     rate = rate,
                     ss = ss,
                     time = time,
                     start = start,
                     end = end,
                     weight = weight,
                     CLHat = CL,
                     QHat = Q,
                     V1Hat = V1,
                     V2Hat = V2,
                     kaHat = ka,
                     ke0Hat = ke0,
                     EmaxHat = Emax,
                     EC50 = EC50,
                     E0Hat = E0,
                     gamma = gamma,
                     nRandom = length(omega),
                     omega = omega,
                     rho = rho,
                     sigma = sigma,
                     sigmaPD = sigmaPD))

### Simulate using Stan
#######################Simulate via rstan#################################
#sim <- stan(file = file.path(modelDir, paste(simModelName, ".stan", sep = "")),
#            data = dataSim,
#            algorithm = "Fixed_param",
#            iter = 1,
#            chains = 1)
#######################Simulate vis cmdstanr #############################
mod.sim <- cmdstan_model(file.path(modelDir, paste(simModelName, ".stan", sep = "")),quiet=F)
sim.cmdstan <- mod.sim$sample(data = dataSim,
                  chains=1,
                  #iter_warmup = 0,
                  iter_sampling = 1,
                  thin = 1,
                  fixed_param = T)
#convert output to a stan fit object
sim <- read_stan_csv(sim.cmdstan$output_files())
################################################################################################
### Assemble data set for fitting via Stan

xdata <- allData %>%
    bind_cols(as.data.frame(sim, pars = "cObs") %>%
                  gather(factor_key = TRUE) %>%
                      select(cObs = value)) %>%
                          bind_cols(as.data.frame(sim, pars = "respObs") %>%
                                        gather(factor_key = TRUE) %>%
                                            select(respObs = value))

xdata <- xdata %>%
    mutate(cObs = ifelse(time %in% xpk & time != 0 & evid == 0, cObs, NA),
           respObs = ifelse(time %in% xpd & evid == 0, respObs, NA))

head(xdata)

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

for(thisDose in doses){
  xplot <- xdata %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot %>% filter(!is.na(cObs)), aes(x = time, y = cObs))
  p1 <- p1 + geom_point() + geom_line() +
    labs(title = paste("dose =", thisDose, "mg"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1)
  
  p1 <- ggplot(xplot %>% filter(!is.na(respObs)), 
               aes(x = time, y = respObs))
  p1 <- p1 + geom_point() + geom_line() +
    labs(title = paste("dose =", thisDose, "mg"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1)
}

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(cObs) & evid == 0])
nObsPK <- length(iObsPK)
## Indices of records containing observed response
iObsPD <- with(xdata, (1:nrow(xdata))[!is.na(respObs) & evid == 0])
nObsPD <- length(iObsPD)

## Parameters for informative priors

CLPrior <- 10
QPrior <- 15
V1Prior <- 35
V2Prior <- 105
kaPrior <- 2
CLPriorCV <- 1 ## 0.10
QPriorCV <- 1 ## 0.18
V1PriorCV <- 1 ## 0.14
V2PriorCV <- 1 ## 0.17
kaPriorCV <- 1 ## 0.16

k10 <- CLPrior / V1Prior
k12 <- QPrior / V1Prior
k21 <- QPrior / V2Prior
ksum <- k10 + k12 + k21
lambda1 <- (ksum + sqrt(ksum^2 - 4 * k10 * k21)) / 2
dkaPrior <- kaPrior - lambda1
dkaPriorCV <- kaPriorCV

ke0Prior <- 0.5
ke0PriorCV <- 1
E0Prior <- 80
E0PriorCV <- 1
EmaxPrior <- 40
EmaxPriorCV <- 1
EC50Prior <- 250
EC50PriorCV <- 1
gammaPrior <- 2
gammaPriorCV <- 1

## create data set
data <- with(xdata,
             list(nId = nId,
                  nt = nt,
                  nObsPK = nObsPK,
                  iObsPK = iObsPK,
                  nObsPD = nObsPD,
                  iObsPD = iObsPD,
                  amt = amt,
                  cmt = cmt,
                  evid = evid,
                  ii = ii,
                  addl = addl,
                  rate = rate,
                  ss = ss,
                  time = time,
                  start = start,
                  end = end,
                  weight = weight,
                  cObs = cObs[iObsPK],
                  respObs = respObs[iObsPD],
                  CLPrior = CLPrior,
                  QPrior = QPrior,
                  V1Prior = V1Prior,
                  V2Prior = V2Prior,
                  dkaPrior = dkaPrior,
                  CLPriorCV = CLPriorCV,
                  QPriorCV = QPriorCV,
                  V1PriorCV = V1PriorCV,
                  V2PriorCV = V2PriorCV,
                  dkaPriorCV = dkaPriorCV,
                  ke0Prior = ke0Prior,
                  ke0PriorCV = ke0PriorCV,
                  E0Prior = E0Prior,
                  E0PriorCV = E0PriorCV,
                  EmaxPrior = EmaxPrior,
                  EmaxPriorCV = EmaxPriorCV,
                  EC50Prior = EC50Prior,
                  EC50PriorCV = EC50PriorCV,
                  gammaPrior = gammaPrior,
                  gammaPriorCV = gammaPriorCV
             ))

### create initial estimates
init <- function(){
    list(CLHat = exp(rnorm(1, log(CLPrior), CLPriorCV)),
         QHat = exp(rnorm(1, log(QPrior), QPriorCV)),
         V1Hat = exp(rnorm(1, log(V1Prior), V1PriorCV)),
         V2Hat = exp(rnorm(1, log(V2Prior), V2PriorCV)),
         dkaHat = exp(rnorm(1, log(dkaPrior), dkaPriorCV)),
         sigma = 0.2,
         ke0Hat = exp(rnorm(1, log(ke0Prior), ke0PriorCV)),
         E0Hat = exp(rnorm(1, log(E0Prior), E0PriorCV)),
         EmaxHat = exp(rnorm(1, log(EmaxPrior), EmaxPriorCV)),
         EC50 = exp(rnorm(1, log(EC50Prior), EC50PriorCV)),
         gamma = exp(rnorm(1, log(gammaPrior), gammaPriorCV)),
         sigmaPD = 0.2,
         omega = exp(rnorm(8, rep(log(0.25), 8), 0.5)),
         L = diag(8),
         etaStd = matrix(rep(0, 8 * nId), nrow = 8)
    )
}

### Specify the variables for which you want history and density plots

parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat",
                      "sigma", "ke0Hat", "E0Hat", "EmaxHat", "EC50", 
                      "gamma", "sigmaPD", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "respObsCond", "cObsPred", "respObsPred",
              "theta")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

####################### build and fit Stan model via rstan #########################
#fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
#            data = data,
#            pars = parameters,
#            iter = nIter,
#            warmup = nBurnin,
#            thin = nThin, 
#            init = init,
#            chains = nChains,
#            control = list(adapt_delta = 0.9))

####################### build and fit Stan model via cmdstanr #########################
nSample <- nIter - nBurnin
#build the model via cmdstan_model
#quiet = F will show all details while building the Stan model
mod <- cmdstan_model(file.path(modelDir, paste(modelName, ".stan", sep = "")), quiet = T)
##fit the model and do MCMC sampling
fit.cmdstan <- mod$sample(data = data,
                          seed = 11191951,
                          refresh = 25, ##the process script will refresh every 25 samplings
                          init = init,
                          chains = nChains,
                          parallel_chains = min(nChains, detectCores()),
                          iter_warmup = nBurnin,
                          iter_sampling = nSample,
			  adapt_delta = 0.9,
                          thin = nThin)
##convert the fit output to a stan fit object
fit <- read_stan_csv(fit.cmdstan$output_files())

################################################################################################

##save the output
dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

## Remove diagonal & redundant elements of rho
dimRho <- nrow(init()$L)
parametersToPlot <- c(parametersToPlot,
                     paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho), 1, paste, collapse = ","),
                                          ncol = dimRho)[upper.tri(diag(dimRho), diag = FALSE)], "]", sep = ""))
parametersToPlot <- setdiff(parametersToPlot, "rho")

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "brightblue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 5, myTheme = myTheme)
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, byChain = TRUE, 
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, 
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))

pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
### posterior predictive distributions

#### prediction of future observations in the same studies, i.e., posterior predictions conditioned on observed data from the same study

pred <- as.data.frame(fit, pars = "cObsCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                      median = quantile(value, probs = 0.5, na.rm = TRUE),
                      ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
                          bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = cObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("individual predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

pred <- as.data.frame(fit, pars = "respObsCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                      median = quantile(value, probs = 0.5, na.rm = TRUE),
                      ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
                          bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = respObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("individual predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

################################################################################################
### posterior predictive distributions

#### prediction of future observations in a new study, i.e., posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = cObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("population predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

pred <- as.data.frame(fit, pars = "respObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = respObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("population predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

dev.off()
