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

modelName <- "ttpn1"
simModelName <- "ttpn1Sim"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
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
library(sampling)
##library(survival)
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

CLHat = 0.0052
omegaCL = 0.25
VHat = 0.25
omegaV = 0.25
sigma = 0.15
omega = c(omegaCL, omegaV)
rho = diag(2)

ke0 = 3.60E-4
alpha = 2.26E-6
beta = 1.37

## Observation and dosing times
dt = 12
tObs <- seq(0, 6 * 21 * 24, by = dt)

## Assemble data set for simulation using Stan
obsData <- data.frame(time = tObs) %>%
    mutate(cmt = 1,
           evid = 0,
           addl = 0,
           ii = 0)

doseData <- data.frame(time = 0) %>%
    mutate(cmt = 1,
           evid = 1,
           addl = 5,
           ii = 21 * 24)

nId = 20 * 3
allData <- doseData %>%
    bind_rows(obsData) %>%
    merge(data.frame(id = 1:nId, amt = rep(c(1.2, 1.8, 2.4), ea = nId / 3))) %>%
    arrange(id, time, desc(evid)) %>%
    mutate(amt = ifelse(evid == 0, 0, amt),
           rate = amt)

nt <- nrow(allData)
start <- (1:nt)[!duplicated(allData$id)]
end <- c(start[-1] - 1, nt)

## Generate sparse PK sampling times by randomly selecting 6 sample times per patient
iPKObs = sampling::strata(allData, stratanames = "id", method = "srswor",
                size = rep(6, length(unique(allData$id))))$ID_unit
nPKObs = length(iPKObs)

dataSim <- with(allData,
                list(nId = nId,
                     nt = nt,
                     nPKObs = nPKObs,
                     iPKObs = iPKObs,
                     amt = 1000 * amt,
                     rate = 0 * rate,
                     cmt = cmt,
                     evid = evid,
                     ii = ii,
                     addl = addl,
                     time = time,
                     start = start,
                     end = end,
                     CLHat = CLHat,
                     VHat = VHat,
                     ke0 = ke0,
                     alpha = alpha,
                     beta = beta,
                     nRandom = length(omega),
                     omega = omega,
                     rho = rho,
                     sigma = sigma))

parameters = c("cObs", "cHat", "ce", "cdf", "CL", "V", "x")

### Simulate using Stan
#######################Simulate via rstan#################################
#sim <- stan(file = file.path(modelDir, paste(simModelName, ".stan", sep = "")),
#            data = dataSim,
#            pars = parameters,
#            algorithm = "Fixed_param",
#            iter = 1,
#            chains = 1)
#######################Simulate vis cmdstanr #############################
mod.sim <- cmdstan_model(file.path(modelDir, paste(simModelName, ".stan", sep = "")), quiet=T)
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

pkData <- allData[iPKObs,] %>%
  bind_cols(as.data.frame(sim, pars = "cObs") %>%
            gather(factor_key = TRUE) %>%
            select(cObs = value)) %>%
  mutate(type = "pk",
         censor = FALSE)

cdf <- allData[c("id", "time")] %>%
  bind_cols(as.data.frame(sim, pars = "cdf") %>%
              gather(factor_key = TRUE) %>%
              select(cdf = value))

## Simulate time to PN using simulated cdf.
tpn <- cdf %>%
  group_by(id) %>%
  summarize(tpn = approx(cdf, time, xout = runif(1))$y,
            censor = is.na(tpn),
            tpn = ifelse(censor, max(time), tpn))

tpn = pkData %>%
  filter(!duplicated(id)) %>%
  select(-time, -cObs, -censor) %>%
  left_join(tpn) %>%
  rename(time = tpn) %>%
  mutate(cObs = NA,
         type = "tpn")
  
doseData <- allData %>%
  filter(evid == 1) %>%
  mutate(cObs = NA,
         censor = FALSE,
         type = "dose")

xdata <- pkData %>%
  bind_rows(tpn) %>%
  bind_rows(doseData) %>%
  arrange(id, time, desc(evid))
  
## Indices of records containing observed concentrations
iPKObs <- with(xdata, (1:nrow(xdata))[type == "pk"])
nPKObs <- length(iPKObs)
## Indices of records containing observed time to PN
iPNObs <- with(xdata, (1:nrow(xdata))[type == "tpn" & !censor])
nPNObs <- length(iPNObs)
## Indices of records containing censored time to PN
iPNCens <- with(xdata, (1:nrow(xdata))[type == "tpn" & censor])
nPNCens <- length(iPNCens)

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$id)]
end <- c(start[-1] - 1, nt)

data <- with(xdata,
                list(nId = nId,
                     nt = nt,
                     nPKObs = nPKObs,
                     iPKObs = iPKObs,
                     nPNObs = nPNObs,
                     iPNObs = iPNObs,
                     nPNCens = nPNCens,
                     iPNCens = iPNCens,
                     amt = 1000 * amt,
                     rate = 0 * rate,
                     cmt = cmt,
                     evid = evid,
                     ii = ii,
                     addl = addl,
                     time = time,
                     start = start,
                     end = end,
                     cObs = cObs[iPKObs]))

### create initial estimates
init <- function(){
  list(CLHat = exp(rnorm(1, log(0.005), 0.25)),
       VHat = exp(rnorm(1, log(0.25), 0.25)),
       sigma = 1,
       alpha = exp(rnorm(1, log(2.0E-6), 0.25)),
       beta = 1 + exp(rnorm(1, log(0.5), 0.25)),
       ke0 = exp(rnorm(1, log(4.0E-4), 0.25)),
       omega = exp(rnorm(2, rep(log(0.25), 2), 0.25)),
       L = diag(2),
       eta = matrix(rep(0, 2 * nId), nrow = 2)
  )
}

### Specify the variables for which you want history and density plots

parametersToPlot <- c("CLHat", "VHat", "sigma", "ke0", "alpha", 
                      "beta", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("CL", "V", "cObsCond", "cObsPred", "cdfPred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 250 ## Number of post-burn-in samples per chain after thinning
nBurn <- 250 ## Number of burn-in samples per chain after thinning
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
#            cores = nChains,
#            refresh = 10,
#            control = list(adapt_delta = 0.95, stepsize = 0.01))
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
                          thin = nThin,
			              step_size = 0.01,
			              adapt_delta = 0.95)
##convert the fit output to a stan fit object
fit <- read_stan_csv(fit.cmdstan$output_files())
################################################################################################

## save the output
dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "brightblue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

dimRho <- nrow(init()$L)
parametersToPlot <- c(parametersToPlot,
                      paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho), 1, paste, collapse = ","),
                                           ncol = dimRho)[upper.tri(diag(dimRho), diag = FALSE)], "]", sep = ""))
parametersToPlot <- setdiff(parametersToPlot, "rho")

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

dev.off()

