
#try(setwd('/home/ditte/Documents/PhD/IV-Structural/StructuralCox/Paper/Biostatistics/SupMat/'), silent=TRUE)

###################################################################################################
# PACKAGES NEEDED

library(rootSolve)
library(timereg)

###################################################################################################
# OPTIMAL TIME

findTimes <- function(m, dd) {
# function to make vector of timepoints in which the optimal time point is to be found
tt <- m$cum[-1,1] # jump times in the data from the model object (m)
tt.max <- min(c(  max(dd$time[dd$X==0&dd$G==0]),
                  max(dd$time[dd$X==0&dd$G==1]),
                  max(dd$time[dd$X==1&dd$G==0]),
                  max(dd$time[dd$X==1&dd$G==1])))
tt.min <- max(c(  min(dd$time[dd$X==0&dd$G==0 & dd$status==1]),
                  min(dd$time[dd$X==0&dd$G==1 & dd$status==1]),
                  min(dd$time[dd$X==1&dd$G==0 & dd$status==1]),
                min(dd$time[dd$X==1&dd$G==1 & dd$status==1])))
tt <- tt[tt<=tt.max & tt>=tt.min]
tt <- tt[unique(round(as.numeric(quantile(1:length(tt), 1:100/101), 0)))] # simplifies the vector to have at most 100 entries
tt
}

optimalTime <- function(dd=d, m) {
# find the optimal time by estimating psi and its asymptotic variance at a range of timepoints, and finding the time with the smallest variance
tt <- findTimes(m, dd=dd)
var <- psivec <- rep(0, length(tt))
i <- 1
for(time in tt) {
    tmp <- try(EstAal(m=m, t.st=time, d=dd), silent=TRUE)
if(length(tmp)!=1) {
    psivec[i] <- var[i] <- NA}
if(length(tmp)==1) {
    psivec[i] <- tmp
    ep <- epU(dd=dd, time=time, psi=psivec[i], max.timepoint.sim=NULL, m=m)
    V <- mean(ep$epU^2)/ep$EDpsi^2/dim(dd)[1] # the variance is calculated
    var[i] <- V
}
i <- i+1
}
index <- !is.na(var)
return(c(tt[index][which.min(var[index])], psivec[index][which.min(var[index])])) # returns the optimal time, and the corresponding psi
}


###################################################################################################
# ESTIMATION

EstimationFct <- function(d=data, S.hat) {
# function to estimate psi from predicted survival probabilities
est.function=function(est){
        val = sum((d$G-mean(d$G))*S.hat^(exp(-est*d$X)))
	return(val)
}
est.function.Vec <- Vectorize(est.function)
psihat.cox <- uniroot.all(est.function.Vec,interval=c(-10,10))
vec <- seq(-15,15, by=0.5)
val <- est.function.Vec(vec)
vec <- seq(vec[which.min(abs(val))]-.1, vec[which.min(abs(val))]+.1, by=0.01)
val <- est.function.Vec(vec)
vec <- seq(vec[which.min(abs(val))]-.1, vec[which.min(abs(val))]+.1, by=0.001)
val <- est.function.Vec(vec)
psihat.cox <- vec[which.min(abs(val))]
return(psihat.cox)
}

EstAal <- function(t.st =1, d=data, m) {
# function to estimate psi from an aalen model (m) and a timepoint (t)
timeIndex <- findInterval(t.st, m$cum[, 'time'])
MM <- model.matrix(formula(m), d)
Gam0 <- as.numeric(m$cum[timeIndex,-1])
Gam0 <- MM %*% matrix(Gam0, ncol=1)
pred <- as.numeric(exp(-Gam0))
est <- EstimationFct(S.hat = pred, d=d)
return(est) }

###################################################################################################
# ASYMPTOTIC VARIANCE

epU <- function(dd=d, time = 5, psi, max.timepoint.sim=101, m) {
# function to get the IID of U
muG <- mean(dd$G)
timeIndex <- findInterval(time, m$cum[, 'time'])
MM <- model.matrix(formula(m), dd)
Gam0 <- as.numeric(m$cum[timeIndex,-1])
Gam0 <- MM %*% matrix(Gam0, ncol=1)
S <- exp(-Gam0)
A <- (dd$G - muG)*S^{exp(-dd$X*psi)}*exp(-psi*dd$X)
Dpsi <- A * Gam0 * dd$X
Index <- findInterval(time, m$cum[, 'time'])
Biid <- do.call(rbind, lapply(m$B.iid, function(x) x[Index, ]))
ep0 <- (dd$G - muG)*(S^(exp(-dd$X*psi)) - mean(S^(exp(-dd$X*psi))))
epB <- mapply(FUN=function(i) mean(A*MM[, i]) * Biid[,i], 1:dim(MM)[2])
epB <- -epB * dim(dd)[1]
epU <- (ep0+rowSums(epB))
EDpsi <- mean(Dpsi)
l <- list()
l$epU <- epU
l$EDpsi <- EDpsi
l$time <- time
return(l)}

CI <- function(dd=d, t.st = 5, psi=0, m) {
# function to calculate variance (and confidence intervals) of the estimated psi
ep <- epU(dd=dd, time=t.st, psi=psi, m=m)
V <- mean(ep$epU^2)
sd0 <- sqrt(V /ep$EDpsi^2 /dim(dd)[1])
ConfInt <- c(psi, psi-qnorm(0.975)*sd0, psi+qnorm(0.975)*sd0)
retur <- matrix(c(ConfInt, sd0), nrow=1)
colnames(retur) <- c('est', 'lower ci', 'upper ci', 'as sd')
return(retur)}

###################################################################################################
# DATA ANALYSIS

load('data.Rdata')

n <- dim(d)[1]
n # sample size

# Association model :
model <- aalen(Surv(time,status==1)~X*G, n.sim=51, data=d, resample.iid=1, max.clust=n)

###################################################################################################
# Estimation at a fixed time (t.fixed)
t.fixed <- 1
psi.fixed <- EstAal(m=model, t=t.fixed, d=d)
psi.fixed # psi estimated at the time given in t.fixed
CI.fixed <- CI(psi=psi.fixed, dd=d, t.st=t.fixed, m=model)
CI.fixed # estimate at the time given in t.fixed, confidence interval and asymptotic sd

###################################################################################################
# Estimating at the optimal time (t.opt)
tt <- optimalTime(dd=d, model)
t.opt <- tt[1]
t.opt # optimal time
psi.opt <- EstAal(m=model, t=t.opt, d=d)
psi.opt # psi estimated at the optimal time
CI.opt <- CI(psi=psi.opt, dd=d, t.st=t.opt, m=model)
CI.opt # estimate at the optimal time, confidence interval and asymptotic sd




