###########################################################################################################################
# Updated R code associated to the paper 
# "On models for the estimation of the excess mortality hazard in case of insufficiently stratified life tables"
# See also: http://rpubs.com/FJRubio/FGH
###########################################################################################################################

rm(list=ls())

# Required packages
library(numDeriv)
library(knitr)
library(compiler)
library(HazReg)


# Reading the simulated data (available at: https://sites.google.com/site/fjavierrubio67/dataFGH.txt)
dat <- read.table("dataFGH.txt", header = FALSE)

# Difference of cumulative hazards
Hp.diff.new <- as.vector(dat[,5])

# Variables used in the modelling
hp.as = as.vector(dat[,4])  # expected rates
ind.cens = as.vector(dat[,6])  # censoring status: 0-alive, 1-dead
surv.time = as.vector(dat[,7]) # survival times
x2 = as.matrix(dat[,1:3]) # design matrix
colnames(x2) <- c("agec","sex","TTT") # names of the covariates
dim.x22 <- dim(x2)[2]

##  Model M1: Classical Model

OPTPGW <- GEHMLE(init = c(0,0,1,rep(0,ncol(x2)*2)), times = surv.time, status = ind.cens, 
                 hp =hp.as, hstr = "GH", dist = "PGW",des = x2, des_t = x2, 
                 method = "nlminb", maxit = 1e4  )

MLEPGW <- c(exp(OPTPGW$OPT$par[1:3]),OPTPGW$OPT$par[-c(1:3)])
AICPGW <- 2*OPTPGW$OPT$objective + 2*length(MLEPGW)

CI <- Conf_Int(OPTPGW$log_lik,OPTPGW$OPT$par,level=0.95)
kable(CI)


## Model M2: Single parameter correction

loglik_pgwC = function(par){
  if(any(is.na(par))) return(Inf)
  alpha <- exp(par[1]); ae0 <- exp(par[2]); be0 <- exp(par[3]); ce0 <- exp(par[4]); beta0 <- par[5:(5+dim(x2)[2]-1)]; beta1 <- tail(par,n = dim(x2)[2]);
  exp.x.beta0 <- exp(x2%*%beta0)
  exp.x.beta1 <- exp(x2%*%beta1)
  exp.x.dif <- exp(x2%*%(beta1-beta0))
  lhaz0= log(alpha*hp.as + hpgw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.beta1)
  lSurv0 = -alpha*(Hp.diff.new) -chpgw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.dif
  val = sum(- ind.cens*lhaz0 -lSurv0)
  ifelse(is.na(val),return(Inf),return(val))
}

# Optimisation
initC <- c(log(3),log(MLEEWG[c(1,2,3)]),MLEEWG[-c(1,2,3)])
OPTC1 = nlminb(initC,loglik_pgwC,control=list(iter.max=10000))
OPTC2 = optim(initC,loglik_pgwC,control=list(maxit=10000))

if(OPTC1$objective <=  OPTC2$value){
  MLEC <- c(exp(OPTC1$par[1:4]),OPTC1$par[-c(1:4)])
  AICC <- 2*OPTC1$objective + 2*length(MLEC)
}

if(OPTC1$objective >  OPTC2$value){
  MLEC <- c(exp(OPTC2$par[1:4]),OPTC2$par[-c(1:4)])
  AICC <- 2*OPTC2$value + 2*length(MLEC)
}

# MLE
kable(MLEC)

# Confidence interval
CIC <- Conf_Int(loglik_pgwC,c(log(MLEC[1:4]),MLEC[-c(1:4)]),level=0.95)
kable(CIC)


## Model M3: Frailty Model

loglik_pgwF = function(par){
  if(any(is.na(par))) return(Inf)
  theta = exp(par[1]); kappa = exp(par[2])/exp(par[1]); ae0 = exp(par[3]); be0 = exp(par[4]); ce0 <- exp(par[5]);
  beta0 <- par[6:(6+dim(x2)[2]-1)]; beta1 <- tail(par,n = dim(x2)[2]);
  exp.x.beta0 <- exp(x2%*%beta0)
  exp.x.beta1 <- exp(x2%*%beta1)
  exp.x.dif <- exp(x2%*%(beta1-beta0))
  haz0= kappa*theta*hp.as/(1+theta*Hp.diff.new) + hpgw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.beta1
  log.Surv0 = -chpgw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.dif - kappa*log(1+theta*Hp.diff.new)
  val = sum(- ind.cens*log(haz0) - log.Surv0)
  ifelse(is.na(val),return(Inf),return(val))
}

# Optimisation step
initE <- c(log(10),log(6.5),log(MLEEWG[c(1,2,3)]),MLEEWG[-c(1,2,3)])
OPTE1 = nlminb(initE,loglik_pgwF,control=list(iter.max=10000))
OPTE2 = optim(initE,loglik_pgwF,control=list(maxit=10000))

if(OPTE1$objective <=  OPTE2$value){
  MLEE <- c(exp(OPTE1$par[1:5]),OPTE1$par[-c(1:5)])
  AICE <- 2*OPTE1$objective + 2*length(MLEE)
}

if(OPTE1$objective >  OPTE2$value){
  MLEE <- c(exp(OPTE2$par[1:5]),OPTE2$par[-c(1:5)])
  AICE <- 2*OPTE2$value + 2*length(MLEE)
}

# MLE
kable(MLEE)

# Confidence interval
CIE <- Conf_Int(loglik_pgwF,c(log(MLEE[1:5]),MLEE[-c(1:5)]),level=0.95)
kable(CIE)



## Comparison of the models

# True values of the parameters
beta1 <- c(0.1,0.1,0.1)
beta2 <- c(0.05,0.2,0.25)
ae <- 1.75; be <- 0.6; ce <- 2.5;

# Comparison with the MLEs
MLES <- cbind(c(NA,NA,MLEEWG), c(NA,MLEC), MLEE, c(10,6.5,ae,be,ce,beta1,beta2))
colnames(MLES) <- c("M1", "M2", "M3", "True")
kable(MLES, digits = 2)

# Comparison using AIC
AICS <- cbind( c("M1", "M2", "M3"), round(c(AICEWG, AICC, AICE),2))
kable(AICS, col.names = c("Model", "AIC"))

# Comparison of the fitted excess baseline hazards
true.haz <- Vectorize(function(t) hew(t,ae,be,ce))
fit1 <- Vectorize(function(t) hpgw(t, MLEEWG[1], MLEEWG[2], MLEEWG[3]))
fit2 <- Vectorize(function(t) hpgw(t, MLEC[2], MLEC[3], MLEC[4]))
fit3 <- Vectorize(function(t) hpgw(t, MLEE[3], MLEE[4], MLEE[5]))

curve(true.haz, 0, 5, lwd = 2, xlab = "Time", ylab ="Baseline Excess Hazard", main = "Fitted Excess Hazards",
      cex.axis = 1.5, cex.lab = 1.5, ylim= c(0.05,0.475), n = 1000)
curve(fit1, 0, 5, lwd = 2, lty = 2, col = "red",add = T, n = 1000)
curve(fit2, 0, 5, lwd = 2, lty = 2, col = "orange",add = T, n = 1000)
curve(fit3, 0, 5, lwd = 2, lty = 2, col = "blue",add = T, n = 1000)
legend(4,0.475, c("True","M1","M2","M3"),
       text.col = c("black","red","orange","blue"), col = c("black", "red", "orange", "blue"), lty = c(1, 2, 2, 2), lwd=2,
       merge = TRUE, bg = "white")
