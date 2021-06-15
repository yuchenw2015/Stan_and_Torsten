## Nonlinear regression example

rm(list = ls())
gc()

modelName <- "fxaInhibitAvg"

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

set.seed(10271998) ## not required but assures repeatable results

################################################################################################

## get data file
xdata <- read.csv(file.path(dataDir, "fxa.data.avg.csv"))

## create data set
data <- with(xdata,
             list(
                 nObs = nrow(xdata),
                 c24 = cavg,
                 fxa24 = fxa.inh.avg
                 ))

## create initial estimates
init <- function()
    list(emax = runif(1, 0, 100),
         ec50 = exp(rnorm(1, log(100), 0.5)),
         gamma = runif(1, 0.1, 4),
         sigma = exp(rnorm(1, log(10), 0.5)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("emax", "ec50", "gamma", "sigma")

## Additional variables to monitor
otherRVs <- c("fxa24Pred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

dir.create(outDir)

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            control = list(adapt_delta = 0.9),
            chains = nChains)

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
            myTheme = theme(text = element_text(size = 12), 
                            axis.text = element_text(size = 10)))
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, 
            myTheme = theme(text = element_text(size = 12), 
                            axis.text = element_text(size = 10)))

pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive distributions

pred <- as.data.frame(fit, pars = "fxa24Pred") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05),
                      median = quantile(value, probs = 0.5),
                      ub = quantile(value, probs = 0.95)) %>%
                          bind_cols(xdata %>% dplyr:::select(cavg, fxa.inh.avg))

p1 <- ggplot(pred, aes(x = cavg, y = fxa.inh.avg))
p1 <- p1 + geom_point() +
    labs(x = "average ME-2 plasma concentration (ng/mL)",
         y = "average factor Xa inhibition (%)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = cavg, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
