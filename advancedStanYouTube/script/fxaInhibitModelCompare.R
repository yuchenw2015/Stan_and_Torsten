rm(list = ls())
gc()

model1Name <- "fxaInhibit1Ncp"
model2Name <- "fxaInhibit2Ncp"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- file.path(projectDir, "model")
out1Dir <- file.path(modelDir, model1Name)
out2Dir <- file.path(modelDir, model2Name)
toolsDir <- file.path(scriptDir, "tools")

.libPaths("lib")

library(rstan)
library(loo)
source(file.path(toolsDir, "stanTools.R"))

set.seed(10271998) ## not required but assures repeatable results

load(file.path(out1Dir, paste(model1Name, "Fit.Rsave", sep = "")))
logLik1 <- extract_log_lik(fit)

load(file.path(out2Dir, paste(model2Name, "Fit.Rsave", sep = "")))
logLik2 <- extract_log_lik(fit)

waic1 <- waic(logLik1)
waic2 <- waic(logLik2)
waic1
waic2
compare(waic1, waic2)

loo1 <- loo(logLik1)
loo2 <- loo(logLik2)
loo1
loo2
compare(loo1, loo2)
