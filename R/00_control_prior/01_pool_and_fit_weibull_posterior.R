# Calculating the control distribution (through MCMC) for the example ----------------------------------------------------------------
#png("ExampleCombinedControl.png", units="in", width=8, height=5, res=700)

Herbst <- read.csv(file = "Papers/DTEPaper/data/Herbst/Doce2.csv")
kmfit1 <- survfit(Surv(Survival.time, Status)~1, data = Herbst)
plot(kmfit1, col = "blue", conf.int = F, xlab = "Time (months)", ylab = "Overall Survival")

Garon <- read.csv(file = "Papers/DTEPaper/data/Garon/Doce1.csv")
kmfit2 <- survfit(Surv(Survival.time, Status)~1, data = Garon)
lines(kmfit2, col = "red", conf.int = F)

Kim <- read.csv(file = "Papers/DTEPaper/data/Kim/Doce3.csv")
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

#MCMCSample <- data.frame(shape = gamma2sample, scale = lambda2sample)

#write.csv(MCMCSample, "MCMCSample.csv", row.names=FALSE)



weibfit <- survreg(Surv(Survival.time, Status)~1, data = combinedDoce, dist = "weibull")
fixedgammac <- as.numeric(exp(-weibfit$icoef[2]))
fixedlambdac <- as.numeric(1/(exp(weibfit$icoef[1])))