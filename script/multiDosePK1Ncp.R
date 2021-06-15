#######################################################################
#### Adapted by Yuchen Wang                     
#### Scripts adapted to build and fit Stan model from cmdstans                      
#### Fit object converted from cmdsran fit to stan fit
#### Stan model adapted due to the function update, see .stan files
#### Date: June/15/2021
#### email: yuchenw2015@gmail.com
#### Based on the PKPD Stan course by Bill Gillespie
#### Link of the original materials: 
#### https://www.metrumrg.com/course/advanced-use-stan-rstan-torsten-
#### pharmacometric-applications/
#######################################################################

rm(list = ls())
gc()

modelName <- "multiDosePK1Ncp"

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
## Go back to default ggplot2 theme that was overridden by bayesplot
library(ggplot2)
theme_set(theme_gray())
library(tidyverse)
library(parallel)
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(10271998) ## not required but assures repeatable results

#Load cmdstanr package and the cmdstan path
#Assuming the cmdstan is under the local Documents/Trosten folder
library(cmdstanr)
set_cmdstan_path(path = "~/User/Documents/Torsten/cmdstan")
################################################################################################

## get data file
xdata <- read.csv(file.path(dataDir, "prob_2fixed2.csv"))
xdata <- xdata %>%
  mutate(EVID = ifelse(AMT > 0, 1, 0),
         CMT = ifelse(AMT > 0, 1, 2),
         DV = ifelse(EVID, NA, DV),
         RATE = 0,
         SS = 0)

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$ID)]
end <- c(start[-1] - 1, nt)
nSubjects = length(unique(xdata$ID))

## Indices of records containing observed concentrations
iObs <- with(xdata, (1:nrow(xdata))[!is.na(DV) & EVID == 0])
nObs <- length(iObs)

## Create data set for simulations
ntPred <- 129
tPredMin <- 0
tPredMax <- 64
dtPred <- (tPredMax - tPredMin) / (ntPred - 1)
tPred <- tPredMin + ((1:ntPred) - 1) * dtPred
data1 <- data.frame(TIME = tPred) %>%
    mutate(AMT = 0,
           RATE = 0,
           CMT = 2,
           EVID = 0,
           II = 0,
           ADDL = 0, SS = 0)

xdataPred <- xdata %>%
    bind_rows(crossing(ID = unique(xdata$ID), data1)) %>%
    arrange(ID, TIME, desc(EVID))

ntPred <- nrow(xdataPred)
startPred <- (1:ntPred)[!duplicated(xdataPred$ID)]
endPred <- c(startPred[-1] - 1, ntPred)

## create data set
data <- c(with(xdata,
             list(
                 nSubjects = nSubjects,
                 nt = nt,
                 nObs = nObs,
                 iObs = iObs,
                 amt = AMT,
                 cmt = CMT,
                 evid = EVID,
                 rate = RATE,
                 ii = II,
                 addl = ADDL,
                 ss = SS,
                 start = start,
                 end = end,
                 time = TIME,
                 cObs = DV[iObs])),
          with(xdataPred,
             list(
                 ntPred = ntPred,
                 amtPred = AMT,
                 cmtPred = CMT,
                 evidPred = EVID,
                 ratePred = RATE,
                 iiPred = II,
                 addlPred = ADDL,
                 ssPred = SS,
                 startPred = startPred,
                 endPred = endPred,
                 tPred = TIME)))

## create initial estimates
init <- function(){
    list(CLHat = exp(rnorm(1, log(8), 0.2)),
         VHat = exp(rnorm(1, log(50), 0.2)),
         kaHat = exp(rnorm(1, log(0.45), 0.2)),
         omega = exp(rnorm(3, log(0.2), 0.5)),
         L = diag(3),
         sigma = runif(1, 0.5, 2),
         eta = matrix(rep(0, 3 * nSubjects), nrow = 3))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "VHat", "kaHat",
                      "sigma", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", "CL", "V", "ka")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
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
                          seed = 10271998,
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
 
dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

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
## posterior predictive distributions

# prediction of future observations in the same studies, i.e., posterior predictions
# conditioned on observed data from the same study

pred <- as.data.frame(fit, pars = "cObsCond") %>%
  gather(factor_key = TRUE) %>%
  mutate(value = ifelse(value == -99, NA, value)) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdataPred)

p1 <- ggplot(pred, aes(x = TIME, y = DV))
p1 <- p1 + geom_point() +
    labs(title = "individual predictions",
         x = "time (h)",
         y = "plasma drug concentration (mcg/mL)") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8)) +
                       facet_wrap(~ ID)
p1 + geom_line(aes(x = TIME, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

# prediction of future observations in a new study, i.e., posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  mutate(value = ifelse(value == -99, NA, value)) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdataPred)

p1 <- ggplot(pred, aes(x = TIME, y = DV))
p1 <- p1 + geom_point() +
    labs(title = "population predictions",
         x = "time (h)",
         y = "plasma drug concentration (mcg/mL)") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8)) +
                       facet_wrap(~ ID)
p1 + geom_line(aes(x = TIME, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
