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

modelName <- "multiDoseME2PK1Ncp"

## Set to TRUE if you want to analyze a small subset of the data for
## demonstration purposes
demo <- TRUE

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
library(tidyverse)
library(parallel)
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(11191951) ## not required but assures repeatable results

#Load cmdstanr package and the cmdstan path
#Assuming the cmdstan is under the local Documents/Trosten folder
library(cmdstanr)
set_cmdstan_path(path = "~/User/Documents/Torsten/cmdstan")
################################################################################################

## get data file
xdata <- read.csv(file.path(dataDir, "fxaNONMEMData.csv"), as.is = TRUE)
xdata <- xdata %>%
  mutate(DV = as.numeric(DV))

## Augment data set with new dosing records based on the ii and addl entries
addl <- xdata %>%
    filter(EVID == 1 & ADDL > 0) %>%
    rowwise() %>%
    do(as.data.frame(.)[rep(1, .$ADDL),]) %>%
    group_by(ID) %>%
    mutate(TIME = TIME + (1:ADDL[1]) * II,
           ADDL = 0,
           II = 0)

xdata <- xdata %>%
    bind_rows(addl) %>%
    arrange(ID, TIME, desc(EVID))

if(demo){
    ## Reduce data set for a demo that runs more quickly
    
    ## IDs for 1st 2 individuals in each dose group of studies 1 & 2
    study12 <- xdata %>%
        filter(STUDY %in% 1:2) %>%
        distinct(ID, .keep_all = TRUE) %>%
        select(ID, AMT, STUDY) %>%
        group_by(AMT, STUDY) %>%
        slice(1:2)
    
    ## IDs for 1st 20 individuals in study 3
    study3 <- xdata %>%
        filter(STUDY %in% 3) %>%
        distinct(ID, .keep_all = TRUE) %>%
        select(ID, AMT, STUDY) %>%
        group_by(AMT, STUDY) %>%
        slice(1:20)
    
    idKeep <- bind_rows(study12, study3)
    
    xdata <- xdata %>%
        filter(ID %in% idKeep$ID)
}

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$ID)]
end <- c(start[-1] - 1, nt)
nSubjects = length(unique(xdata$ID))

## Indices of records containing observed concentrations (excludes dosing records and BQL data)
iObs <- with(xdata, (1:nrow(xdata))[!is.na(DV) & EVID == 0])
nObs <- length(iObs)

## create data set
data <- with(xdata,
             list(
                 nSubjects = nSubjects,
                 nt = nt,
                 nObs = nObs,
                 iObs = iObs,
                 amt = 1000 * AMT,
                 cmt = CMT,
                 evid = EVID,
                 rate = RATE,
                 ii = II,
                 addl = ADDL,
                 ss = SS,
                 start = start,
                 end = end,
                 time = TIME,
                 cObs = DV[iObs],
                 weight = WEIGHT[!duplicated(ID)]
))

## create initial estimates
init <- function(){
    list(CLHat = exp(rnorm(1, log(10), 0.2)),
         QHat = exp(rnorm(1, log(15), 0.2)),
         V1Hat = exp(rnorm(1, log(35), 0.2)),
         V2Hat = exp(rnorm(1, log(105), 0.2)),
         kaHat = exp(rnorm(1, log(2), 0.2)),
         omega = exp(rnorm(5, log(0.25), 0.5)),
         L = diag(5),
         sigma = runif(1, 0.5, 2),
         eta = matrix(rep(0, 5 * nSubjects), nrow = 5))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat",
                      "sigma", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", "CL", "Q", "V1", "V2", "ka")

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
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
              median = quantile(value, probs = 0.5, na.rm = TRUE),
              ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
    bind_cols(xdata)

xplot <- subset(pred, STUDY == 1)
doses <- sort(unique(xplot$DOSE))
for(thisDose in doses){
    p1 <- ggplot(subset(xplot, DOSE == thisDose), aes(x = TIME, y = DV))
    p1 <- p1 + geom_point() +
        labs(title = paste("study 1", thisDose, "mg\n individual predictions"),
             x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) +
        facet_wrap(~ ID)
    print(p1 + geom_line(aes(x = TIME, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, STUDY == 2)
doses <- sort(unique(xplot$DOSE))
for(thisDose in doses){
    p1 <- ggplot(subset(xplot, DOSE == thisDose), aes(x = TIME, y = DV))
    p1 <- p1 + geom_point() +
        labs(title = paste("study 2", thisDose, "mg\n individual predictions"),
             x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) +
        facet_wrap(~ ID)
    print(p1 + geom_line(aes(x = TIME, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, STUDY == 3)

nPerPage = 25
subjects <- sort(unique(xplot$ID))
nSubjects <- length(subjects)
nPages <- ceiling(nSubjects / nPerPage)
subjects <- data.frame(ID = subjects,
                       page = sort(rep(1:nPages, length = nSubjects)),
                       stringsAsFactors = FALSE)
xplot <- xplot %>% left_join(subjects)

for(thisPage in 1:nPages){
    p1 <- ggplot(subset(xplot, page == thisPage), aes(x = TIME, y = DV))
    p1 <- p1 + geom_point() +
        labs(title = "study 3 20 mg\n individual predictions",
             x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) +
        facet_wrap(~ ID)
    print(p1 + geom_line(aes(x = TIME, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

# prediction of future observations in a new study, i.e., posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
              median = quantile(value, probs = 0.5, na.rm = TRUE),
              ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
    bind_cols(xdata)

xplot <- subset(pred, STUDY == 1)
doses <- sort(unique(xplot$DOSE))
for(thisDose in doses){
    p1 <- ggplot(subset(xplot, DOSE == thisDose), aes(x = TIME, y = DV))
    p1 <- p1 + geom_point() +
        labs(title = paste("study 1", thisDose, "mg\n population predictions"),
             x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) +
        facet_wrap(~ ID)
    print(p1 + geom_line(aes(x = TIME, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, STUDY == 2)
doses <- sort(unique(xplot$DOSE))
for(thisDose in doses){
    p1 <- ggplot(subset(xplot, DOSE == thisDose), aes(x = TIME, y = DV))
    p1 <- p1 + geom_point() +
        labs(title = paste("study 2", thisDose, "mg\n population predictions"),
             x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) +
        facet_wrap(~ ID)
    print(p1 + geom_line(aes(x = TIME, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, STUDY == 3)

nPerPage = 25
subjects <- sort(unique(xplot$ID))
nSubjects <- length(subjects)
nPages <- ceiling(nSubjects / nPerPage)
subjects <- data.frame(ID = subjects,
                       page = sort(rep(1:nPages, length = nSubjects)),
                       stringsAsFactors = FALSE)
xplot <- xplot %>% left_join(subjects)

for(thisPage in 1:nPages){
    p1 <- ggplot(subset(xplot, page == thisPage), aes(x = TIME, y = DV))
    p1 <- p1 + geom_point() +
        labs(title = "study 3 20 mg\n population predictions",
             x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
        theme(text = element_text(size = 12), axis.text = element_text(size = 12),
              legend.position = "none", strip.text = element_text(size = 8)) +
        facet_wrap(~ ID)
    print(p1 + geom_line(aes(x = TIME, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

dev.off()
