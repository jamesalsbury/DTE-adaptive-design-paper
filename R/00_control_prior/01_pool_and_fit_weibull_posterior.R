#!/usr/bin/env Rscript
#
# Control-arm prior, step 1 of 2 (Section 5.1).
#
# Pools the reconstructed docetaxel control-arm IPD from ZODIAC, REVEL and
# INTEREST and fits a Weibull model to the pooled data by MCMC (JAGS).
# The posterior samples are saved for 02_fit_beta_approximation.R.
#
# Run from this folder:  Rscript 01_pool_and_fit_weibull_posterior.R
#
# Produces:
#   - weibull_posterior_samples.rds (lambda2sample, gamma2sample)
#   - Rplots.pdf (KM curves per trial and pooled; MCMC trace plots)

library(survival)
library(rjags)

# Calculating the control distribution (through MCMC) for the example ----------------------------------------------------------------
#png("ExampleCombinedControl.png", units="in", width=8, height=5, res=700)

Herbst <- read.csv(file = "../../data/zodiac_control_arm.csv")
kmfit1 <- survfit(Surv(Survival.time, Status)~1, data = Herbst)
plot(kmfit1, col = "blue", conf.int = F, xlab = "Time (months)", ylab = "Overall Survival")

Garon <- read.csv(file = "../../data/revel_control_arm.csv")
kmfit2 <- survfit(Surv(Survival.time, Status)~1, data = Garon)
lines(kmfit2, col = "red", conf.int = F)

Kim <- read.csv(file = "../../data/interest_control_arm.csv")
kmfit3 <- survfit(Surv(Survival.time, Status)~1, data = Kim)
lines(kmfit3, col = "yellow", conf.int = F)

legend("topright", legend = c("ZODIAC", "REVEL", "INTEREST"), lty = 1, col = c("blue", "red", "yellow"))
#dev.off()
combinedDoce <- rbind(Herbst, Garon, Kim)

kmfit4 <- survfit(Surv(Survival.time, Status)~1, data = combinedDoce)
lines(kmfit4, col = "black", conf.int = F)


#Performing MCMC on this data set

modelstring="

data {
  for (j in 1:n){
    zeros[j] <- 0
  }
}

model {
  C <- 10000
  for (i in 1:n){
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, log(gamma2)+gamma2*log(lambda2*datTimes[i])-(lambda2*datTimes[i])^gamma2-log(datTimes[i]), -(lambda2*datTimes[i])^gamma2)
  }

    lambda2 ~ dnorm(1,1/10000)T(0,)
    gamma2 ~ dnorm(1,1/10000)T(0,)

    }
"

model = jags.model(textConnection(modelstring), data = list(datTimes = combinedDoce$Survival.time, datEvent = combinedDoce$Status, n= nrow(combinedDoce)), quiet = T)

update(model, n.iter=1000)
output=coda.samples(model=model, variable.names=c("lambda2", "gamma2"), n.iter = 10000)

plot(output)

lambda2sample <- as.numeric(unlist(output[,2]))
gamma2sample <- as.numeric(unlist(output[,1]))

saveRDS(list(lambda2sample = lambda2sample, gamma2sample = gamma2sample),
        "weibull_posterior_samples.rds")
cat("Saved posterior samples to weibull_posterior_samples.rds\n")

#MCMCSample <- data.frame(shape = gamma2sample, scale = lambda2sample)

#write.csv(MCMCSample, "MCMCSample.csv", row.names=FALSE)



weibfit <- survreg(Surv(Survival.time, Status)~1, data = combinedDoce, dist = "weibull")
fixedgammac <- as.numeric(exp(-weibfit$icoef[2]))
fixedlambdac <- as.numeric(1/(exp(weibfit$icoef[1])))
cat(sprintf("Pooled Weibull MLE: lambda_c = %.8f, gamma_c = %.6f\n", fixedlambdac, fixedgammac))
