rm(list = ls())
gc()

modelName <- "fluindione1"

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

.libPaths("lib")

library(rstan)
library(bayesplot)
## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())
library(tidyverse)
library(parallel)
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(11191951) ## not required but assures repeatable results

################################################################################################

## get data file
xdata <- read.csv(file.path(dataDir, "fluindione.csv"), as.is = TRUE)

## create data set
data <- with(xdata,
             list(nt = nrow(xdata),
                  dose = 20,
                  time = time,
                  cObs = cObs,
                  inrObs = inrObs
                  ))

## create initial estimates
init <- function(){
    list(CL = exp(rnorm(1, log(0.15), 1)),
         Q = exp(rnorm(1, log(0.2), 1)),
         V1 = exp(rnorm(1, log(5), 1)),
         V2 = exp(rnorm(1, log(5), 1)),
         ka = exp(rnorm(1, log(2), 1)),
         Emax = rbeta(1, 3, 1),
         EC50 = exp(rnorm(1, log(2), 1)),
         inr0 = exp(rnorm(1, log(1), 0.25)),
         kout = exp(rnorm(1, log(0.03), 1)),
         gamma = exp(rnorm(1, log(2), 0.5)),
         sigmaPK = runif(1, 0.5, 2),
         sigmaPD = runif(1, 0.5, 2))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka",
                      "Emax", "EC50", "inr0", "kout", "gamma",
                      "sigmaPK", "sigmaPD")

## Additional variables to monitor
otherRVs <- c("cObsPred", "inrObsPred") 

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            refresh = 10,
            chains = nChains,
            control = list(adapt_delta = 0.95, stepsize = 0.01))

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

pred <- as.data.frame(fit, pars = "cObsPred") %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
              median = quantile(value, probs = 0.5, na.rm = TRUE),
              ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
    bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
    labs(x = "time (h)",
         y = "fluindione plasma concentration (mg/L)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
      geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)


pred <- as.data.frame(fit, pars = "inrObsPred") %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
              median = quantile(value, probs = 0.5, na.rm = TRUE),
              ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
    bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = inrObs))
p1 <- p1 + geom_point() +
    labs(x = "time (h)",
         y = "INR") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
      geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
