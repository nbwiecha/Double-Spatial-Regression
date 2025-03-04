# Two-Stage Estimators for Spatial Confounding with Point-Referenced Data
# Code by Nate Wiecha, North Carolina State University

# Real data analysis: Fayetteville participants of GenX Cohort Study
# Data is not included to protect privacy of participants
# Available on reasonable request to corresponding author

# Although the data is not available, providing this code as example of how 
# DSR can be used in applied studies.

rm(list=ls())
library(tidyverse)
library(GpGp)
library(car)
library(viridis)
setwd("~/GitHub/Double-Spatial-Regression")
# load("data//genx//final_samples.Rdata")
# geocodes <- read.csv("data//genx//Nate_FayPitt_Geocodes.csv")
source("code//dml_function_vector.R")
source("code//dml_function_vector_nosplit.R")

source("code//predict_functions.R")

# Select data from Fayetteville Private Well Community, limit to women, join 
# with lat/lon geocodes, discretize race variable

fay.samples <- data.all %>%
  filter(community=="Fayetteville") %>%
  left_join(geocodes) %>%
  filter(gender=="F") %>%
  mutate(white=race=="White") %>%
  drop_na()
# colnames(fay.samples)

# Get PFAS detection rates in sample
pfas.detection.rates <- fay.samples %>%
  left_join(data2)
colnames(pfas.detection.rates)

pfas.detection.rates <- pfas.detection.rates[c(2:6, 41:79)]
pfas.mins <- apply(pfas.detection.rates, FUN=min, MARGIN=2)
for(j in 1:ncol(pfas.detection.rates)){
  pfas.detection.rates[,j] <- pfas.detection.rates[,j] > pfas.mins[j]
}
colMeans(pfas.detection.rates)

# Summary stats
summary(fay.samples[,c("PFHpS", "PFHxS", "PFNA", "PFOS", "PFOA", "TSH", "age", "race", "smoke", "drinker_status")])
apply(fay.samples[,c("PFHpS", "PFHxS", "PFNA", "PFOS", "PFOA", "TSH", "age")], FUN=sd, MARGIN=2)

summary(log(fay.samples$TSH))
sd(log(fay.samples$TSH))

# Geographic plots
ggplot(data=fay.samples, aes(x=long, y=lat)) +
  geom_point(aes(color=TSH)) +
  scale_color_viridis() +
  ggtitle("Free T4")

ggplot(data=fay.samples, aes(x=long, y=lat)) +
  geom_point(aes(color=PFNA)) +
  scale_color_viridis() +
  ggtitle("PFNA")

ggplot(data=fay.samples, aes(x=long, y=lat)) +
  geom_point(aes(color=PFOA)) +
  scale_color_viridis() +
  ggtitle("PFOA")

ggplot(data=fay.samples, aes(x=long, y=lat)) +
  geom_point(aes(color=PFHpS)) +
  scale_color_viridis() +
  ggtitle("PFHpS")

ggplot(data=fay.samples, aes(x=long, y=lat)) +
  geom_point(aes(color=PFHxS)) +
  scale_color_viridis() +
  ggtitle("PFHxS")

ggplot(data=fay.samples, aes(x=long, y=lat)) +
  geom_point(aes(color=PFOS)) +
  scale_color_viridis() +
  ggtitle("PFOS")

# Check linear model assumptions using non-spatial model (OLS)

lm.fit <- lm(log(TSH) ~ PFNA + PFOA + PFOS + PFHpS + PFHxS +
               age +
               smoke +drinker_status +race ,
               # white,
               # smoke,
             data=fay.samples)
round(summary(lm.fit)$coefficients,2)
residualPlots(lm.fit, quadratic=FALSE)
plot(lm.fit)
plot(hatvalues(lm.fit))
which.max(hatvalues(lm.fit))
plot(cooks.distance(lm.fit))
which.max(cooks.distance(lm.fit))

# Checking OLS results excluding a high-ish leverage point
lm.fit.sens <- lm(log(TSH) ~ PFNA + PFOA + PFOS + PFHpS + PFHxS +
                    age + smoke + drinker_status + race, data=fay.samples[-31,])

summary(lm.fit.sens)
residualPlots(lm.fit.sens)
plot(lm.fit.sens)
plot(hatvalues(lm.fit.sens))
which.max(hatvalues(lm.fit.sens))


lm.fit2 <- lm(TSH ~ PFNA + PFOA + PFOS + PFHpS + PFHxS +
                age + smoke + drinker_status + race, data=fay.samples[-31,])

summary(lm.fit2)
residualPlots(lm.fit2)
plot(lm.fit2)
plot(hatvalues(lm.fit2))
plot(cooks.distance(lm.fit2))
hist(residuals(lm.fit2))
hist(residuals(lm.fit))

# Get results using spatial linear mixed model (LMM)
# fit using package GpGp

Y <- log(fay.samples$TSH)

AX <- model.matrix(lm.fit)
A <- AX[,2:6,drop=FALSE]
X <- AX[,-(1:6),drop=FALSE]
S <- as.matrix(fay.samples[, c("lat", "long")])
n <- length(Y)

round(cor(A), 2)

lmm.fit <- fit_model(y=Y, locs=S, X=AX)
summary(lmm.fit)
#(variance, range, smoothness, nugget)
lmm.fit$covparms

# Vulnerability to spatial confounding is highest when both treatment and 
# response exhibit spatial variation. To check this, determine how much
# of the residual variance in a linear model can be attributed to the spatial
# random effect. If a decent proportion of the residual variance is spatial,
# there is potential vulnerability to spatial confounding. 

# Proportion of residual variance that is spatial, from a model fit:
get_prop_spatial <- function(lmm.fit){
  nug <- lmm.fit$covparms[4]
  spat <- lmm.fit$covparms[1]
  
  1 - nug*spat / (spat + nug*spat)
  
}

get_prop_spatial(lmm.fit) # about 0.41

# Check this proportion for each PFAS
p <- ncol(A)
fitD <- list(length=p)
for(j in 1:p){
  fitD[[j]] <- fit_model(y=A[,j], locs=S, X=cbind(rep(1, n), X), silent=TRUE)
}

fitD[[1]]$covparms 
get_prop_spatial(fitD[[1]]) # PFNA: 0.33
fitD[[2]]$covparms 
get_prop_spatial(fitD[[2]]) # PFOA: 0.46
fitD[[3]]$covparms
get_prop_spatial(fitD[[3]]) # PFOS: 0.17
fitD[[4]]$covparms 
get_prop_spatial(fitD[[4]]) # PFHpS: 0.70
fitD[[5]]$covparms 
get_prop_spatial(fitD[[5]]) # PFHxS: 0.34

# DSR estimation
# Lots of folds and repeated splits to make sure don't use a weird random split

# Let's do 11 times and take the median
M <- 11
beta_hats <- matrix(nrow=M, ncol=5)
sd_beta_hats <- matrix(nrow=M, ncol=5)
for(m in 1:M){
  print(m)
  # Use the DSR function with K=45 folds, refitting the entire model for each fold
  # With only 98 observations, not too time consuming
  dml.fit.re <- dml_fit(Y, D=A, S=S, X=X, K=45, refit=TRUE)
  beta_hats[m,] <- as.numeric(dml.fit.re[[1]])
  sd_beta_hats[m,] <- sqrt(diag(dml.fit.re[[2]]))
}
save.image("outputs//analysis_06202024.Rdata")
load("outputs//analysis_06202024.Rdata")
beta_hat_median <- apply(beta_hats, 2, median)
sd1 <- (sd_beta_hats[beta_hats[,1]==beta_hat_median[1],1])
sd2 <- (sd_beta_hats[beta_hats[,2]==beta_hat_median[2],2])
sd3 <- (sd_beta_hats[beta_hats[,3]==beta_hat_median[3],3])
sd4 <- (sd_beta_hats[beta_hats[,4]==beta_hat_median[4],4])
sd5 <- (sd_beta_hats[beta_hats[,5]==beta_hat_median[5],5])

sd_beta_hat_median <- c(sd1, sd2, sd3, sd4, sd5)
ts_median <- beta_hat_median / sd_beta_hat_median
ts_median

# Forest plot

# DSR results
forplot.lb.dml <- beta_hat_median - 2*sd_beta_hat_median
forplot.ub.dml <- beta_hat_median + 2*sd_beta_hat_median
forplot.dat.dml <- data.frame(beta=beta_hat_median, lb=forplot.lb.dml, ub=forplot.ub.dml, 
                          PFAS=as.factor(c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS")),
                          method="DSR")

# LMM results
beta_hat_lmm <- lmm.fit$betahat[colnames(AX) %in% c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS") ]
forplot.lb.lmm <- beta_hat_lmm - 2*lmm.fit$sebeta[colnames(AX) %in%c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS") ]
forplot.ub.lmm <- beta_hat_lmm + 2*lmm.fit$sebeta[colnames(AX) %in%c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS") ]
forplot.dat.lmm <- data.frame(beta=beta_hat_lmm, lb=forplot.lb.lmm, ub=forplot.ub.lmm, 
                              PFAS=as.factor(c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS")),
                              method="LMM")

# OLS results
beta_hat_ols <- lm.fit$coefficients[names(lm.fit$coefficients) %in% c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS")]
forplot.lb.ols <- confint(lm.fit)[rownames(confint(lm.fit)) %in% c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS"), 1]
forplot.ub.ols <- confint(lm.fit)[rownames(confint(lm.fit)) %in% c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS"), 2]
forplot.dat.ols <- data.frame(beta=beta_hat_ols, lb=forplot.lb.ols, ub=forplot.ub.ols, 
                              PFAS=as.factor(c("PFNA", "PFOA", "PFOS", "PFHpS", "PFHxS")),
                              method="OLS")

# Forest plot
forplot.dat <- bind_rows(forplot.dat.dml, forplot.dat.lmm, forplot.dat.ols) %>%
  mutate(Method=as.factor(method))

forest.plot <- ggplot(forplot.dat) +
  geom_point(aes(x=beta, y=PFAS, color=Method), position=position_dodge(width=.5), size=4) +
  geom_linerange(aes(xmin=lb, xmax=ub, y=PFAS, color=Method, linetype=Method), position=position_dodge(width=.5),
                 linewidth=1.25) +
  geom_vline(xintercept=0) +
  xlab(bquote(beta)) +
  theme_classic() +
  theme(text=element_text(size=25),
        legend.key.size = unit(1.5, "cm"))

png("outputs//forest_plot.png", width = 900, height = 600)
forest.plot
dev.off()
