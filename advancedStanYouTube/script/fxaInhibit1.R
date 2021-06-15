rm(list = ls())
gc()

modelName <- "fxaInhibit1"

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
xdata <- read.csv(file.path(dataDir, "fxa.data.csv"))

## create data set
data <- with(xdata,
             list(
                 nSubjects = max(subject),
                 nObs = nrow(xdata),
                 subject = subject,
                 cObs = cobs,
                 fxa = fxa.inh.obs,
                 cmin = 0, cmax = 1600, nsim = 201))

## create initial estimates
init <- function() list(
    emax = runif(1, 40, 100),
    ec50Hat = exp(rnorm(1, log(100), 1)),
    gamma = runif(1, 0.25, 3),
    sigma = exp(rnorm(1, log(10), 1)),
    omegaEc50 = exp(rnorm(1, log(0.2), 1)),
    logEc50 = rep(log(100), data$nSubjects))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("emax", "ec50Hat", "gamma", "sigma", "omegaEc50")

## Additional variables to monitor
otherRVs <- c("ec50", "fxaCond", "fxaPred", "fxaPred2")

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
            chains = nChains)#,
            #control = list(adapt_delta = 0.95))

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

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

pairs(fit, pars = parametersToPlot)

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive distributions

# prediction of future observations in the same studies, i.e., posterior predictions
# conditioned on observed data from the same study

pred <- as.data.frame(fit, pars = "fxaCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05),
                      median = quantile(value, probs = 0.5),
                      ub = quantile(value, probs = 0.95)) %>%
                          bind_cols(xdata)

doses <- sort(unique(xdata$dose))
for(thisDose in doses){
    plotdata <- subset(pred, dose == thisDose)
    p1 <- ggplot(plotdata, aes(x = cobs, y = fxa.inh.obs))
    p1 <- p1 + geom_point() +
        labs(title = paste(thisDose, "mg\n individual predictions"),
             x = "ME-2 plasma concentration (ng/mL)",
             y = "factor Xa inhibition (%)") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                       legend.position = "none", strip.text = element_text(size = 8)) +
                           facet_wrap(~ subject)
    print(p1 + geom_line(aes(x = cobs, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

for(thisDose in doses){
    plotdata <- subset(pred, dose == thisDose)
    p1 <- ggplot(plotdata, aes(x = time, y = fxa.inh.obs))
    p1 <- p1 + geom_point() +
        labs(title = paste(thisDose, "mg\n individual predictions"),
             x = "time (h)",
             y = "factor Xa inhibition (%)") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                       legend.position = "none", strip.text = element_text(size = 8)) +
                           facet_wrap(~ subject)
    print(p1 + geom_line(aes(x = time, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

# prediction of future observations in a new study, i.e., posterior predictive distributions

pred <- as.data.frame(fit, pars = "fxaPred") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05),
                      median = quantile(value, probs = 0.5),
                      ub = quantile(value, probs = 0.95)) %>%
                          bind_cols(xdata)

doses <- sort(unique(xdata$dose))
for(thisDose in doses){
    plotdata <- subset(pred, dose == thisDose)
    p1 <- ggplot(plotdata, aes(x = cobs, y = fxa.inh.obs))
    p1 <- p1 + geom_point() +
        labs(title = paste(thisDose, "mg\n population predictions"),
             x = "ME-2 plasma concentration (ng/mL)",
             y = "factor Xa inhibition (%)") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                       legend.position = "none", strip.text = element_text(size = 8)) +
                           facet_wrap(~ subject)
    print(p1 + geom_line(aes(x = cobs, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

for(thisDose in doses){
    plotdata <- subset(pred, dose == thisDose)
    p1 <- ggplot(plotdata, aes(x = time, y = fxa.inh.obs))
    p1 <- p1 + geom_point() +
        labs(title = paste(thisDose, "mg\n population predictions"),
             x = "time (h)",
             y = "factor Xa inhibition (%)") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                       legend.position = "none", strip.text = element_text(size = 8)) +
                           facet_wrap(~ subject)
    print(p1 + geom_line(aes(x = time, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

## Posterior predictions plotted as a function of concentration

cPred <- with(data, cmin + (0:(nsim-1)) * (cmax - cmin) / (nsim - 1))
pred <- as.data.frame(fit, pars = "fxaPred2") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05),
                      median = quantile(value, probs = 0.5),
                      ub = quantile(value, probs = 0.95)) %>%
                          mutate(cPred = cPred)
                          
p1 <- ggplot(xdata, aes(x = cobs, y = fxa.inh.obs))
p1 <- p1 + geom_point() +
    labs(title = "population predictions",
         x = "ME-2 plasma concentration (ng/mL)",
         y = "factor Xa inhibition (%)") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(data = pred, aes(x = cPred, y = median)) +
          geom_ribbon(data = pred, aes(x = cPred, y = median, ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
