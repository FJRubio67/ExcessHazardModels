###########################################################################################################################
# R code associated to the paper 
# "On models for the estimation of the excess mortality hazard in case of insufficiently stratified life tables"
# See also: http://rpubs.com/FJRubio/FGH
###########################################################################################################################

rm(list=ls())

# Required packages
library(numDeriv)
library(knitr)
library(compiler)

# Exponentiated Weibull Distribution
# Hazard
hexpw <- function(t,lambda,kappa,alpha){
  pdf0 <-  alpha*dweibull(t,scale=lambda,shape=kappa)*pweibull(t,scale=lambda,shape=kappa)^(alpha-1) 
  cdf0 <- pweibull(t,scale=lambda,shape=kappa)^alpha
  cdf0 <- ifelse(cdf0==1,0.9999999,cdf0)
  return(pdf0/(1-cdf0))
}                                                                                      

# Cumulative hazard
Hexpw <- function(t,lambda,kappa,alpha,log.p=FALSE){
  cdf <- pweibull(t,scale=lambda,shape=kappa)^alpha  
  return(-log(1-cdf))
} 

###########################################################################################
# Function to calculate the normal confidence intervals
# The parameters indicated with "index" are transformed to the real line using log()
###########################################################################################
# FUN   : minus log-likelihood function to be used to calculate the confidence intervals
# MLE   : maximum likelihood estimator of the parameters of interest
# level : confidence level
# index : position of the positive parameters under the original parameterisation

Conf.Int <- function(FUN,MLE,level=0.95){
  sd.int <- abs(qnorm(0.5*(1-level)))
  HESS <- hessian(FUN,x=MLE)
  Fisher.Info <- solve(HESS, tol = 1e-18)
  Sigma <- sqrt(diag(Fisher.Info))
  U<- MLE + sd.int*Sigma
  L<- MLE - sd.int*Sigma
  C.I <- cbind(L,U,MLE, Sigma)
  rownames(C.I)<- paste0("par", seq_along(1:length(MLE)))
  colnames(C.I)<- c("Lower","Upper","Transf MLE", "Std. Error")
  return(C.I)
}

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

log.likewG <- cmpfun(function(par){
  if(anyNA(par)) return(Inf)
  ae0 <- exp(par[1]); be0 <- exp(par[2]); ce0=exp(par[3]); beta0 <- par[4:(4+dim(x2)[2]-1)]; beta1 <- tail(par,n = dim(x2)[2]);
  exp.x.beta0 <- exp(x2%*%beta0)
  exp.x.beta1 <- exp(x2%*%beta1)
  exp.x.dif <- exp(x2%*%(beta1-beta0))
  haz0 <- hp.as + hexpw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.beta1
  CH0 <- Hexpw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.dif + Hp.diff.new
  if(anyNA(CH0)) return(Inf)
  if(anyNA(haz0)) return(Inf)
  else  val <- - ind.cens*log(haz0) + CH0
  return(sum(val))
})

#----------------------------------------------------------------------------------------------------
# nlminb within Coordinate Descent algorithm
# Used to calculate a better initial point
#----------------------------------------------------------------------------------------------------

NCD <- 1 # Number of coordinate descent iterations

init.cd <- c(0,0,0.5,rep(0,3),rep(0,3)) # initial point for coordinate descent

for(j in 1:NCD){
  #print(j)
  for(i in 1:length(init.cd)){
    tempf <- Vectorize(function(par){
      val <- replace(init.cd, i, par)
      return(log.likewG(val))
    })
    
    new.val <- nlminb(init.cd[i],tempf)$par
    init.cd <- replace(init.cd, i,  new.val )
  }}  


#---------------------------------------------------------------------------------------------------- 
# Optimisation

if(any(is.na(init.cd))) {
  initG <- c(0,0,0.5,rep(0,3),rep(0,3))}  else initG <- init.cd

OPTEWG1 <- nlminb(initG, log.likewG,control=list(iter.max=10000))
OPTEWG2 <- optim(initG, log.likewG,control=list(maxit=10000))

if(OPTEWG1$objective <=  OPTEWG2$value){
  MLEEWG <- c(exp(OPTEWG1$par[1:3]),OPTEWG1$par[-c(1:3)])
  AICEWG <- 2*OPTEWG1$objective + 2*length(MLEEWG)
}

if(OPTEWG1$objective >  OPTEWG2$value){
  MLEEWG <- c(exp(OPTEWG2$par[1:3]),OPTEWG2$par[-c(1:3)])
  AICEWG <- 2*OPTEWG2$value + 2*length(MLEEWG)
}

# MLE
kable(MLEEWG)

# Confidence interval (the parameters of the Exponentiated Weibull baseline hazard are in the log-scale)
CI <- Conf.Int(log.likewG,c(log(MLEEWG[1:3]),MLEEWG[-c(1,2,3)]),level=0.95)
kable(CI)


## Model M2: Single parameter correction

log.likEWC = function(par){
  if(any(is.na(par))) return(Inf)
  alpha <- exp(par[1]); ae0 <- exp(par[2]); be0 <- exp(par[3]); ce0 <- exp(par[4]); beta0 <- par[5:(5+dim(x2)[2]-1)]; beta1 <- tail(par,n = dim(x2)[2]);
  exp.x.beta0 <- exp(x2%*%beta0)
  exp.x.beta1 <- exp(x2%*%beta1)
  exp.x.dif <- exp(x2%*%(beta1-beta0))
  lhaz0= log(alpha*hp.as + hexpw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.beta1)
  lSurv0 = -alpha*(Hp.diff.new) -Hexpw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.dif
  val = sum(- ind.cens*lhaz0 -lSurv0)
  ifelse(is.na(val),return(Inf),return(val))
}

# Optimisation
initC <- c(log(3),log(MLEEWG[c(1,2,3)]),MLEEWG[-c(1,2,3)])
OPTC1 = nlminb(initC,log.likEWC,control=list(iter.max=10000))
OPTC2 = optim(initC,log.likEWC,control=list(maxit=10000))

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
CIC <- Conf.Int(log.likEWC,c(log(MLEC[1:4]),MLEC[-c(1:4)]),level=0.95)
kable(CIC)


## Model M3: Frailty Model

log.lik.E = function(par){
  if(any(is.na(par))) return(Inf)
  theta = exp(par[1]); kappa = exp(par[2])/exp(par[1]); ae0 = exp(par[3]); be0 = exp(par[4]); ce0 <- exp(par[5]);
  beta0 <- par[6:(6+dim(x2)[2]-1)]; beta1 <- tail(par,n = dim(x2)[2]);
  exp.x.beta0 <- exp(x2%*%beta0)
  exp.x.beta1 <- exp(x2%*%beta1)
  exp.x.dif <- exp(x2%*%(beta1-beta0))
  haz0= kappa*theta*hp.as/(1+theta*Hp.diff.new) + hexpw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.beta1
  log.Surv0 = -Hexpw(surv.time*exp.x.beta0,ae0,be0,ce0)*exp.x.dif - kappa*log(1+theta*Hp.diff.new)
  val = sum(- ind.cens*log(haz0) - log.Surv0)
  ifelse(is.na(val),return(Inf),return(val))
}

# Optimisation step
initE <- c(log(10),log(6.5),log(MLEEWG[c(1,2,3)]),MLEEWG[-c(1,2,3)])
OPTE1 = nlminb(initE,log.lik.E,control=list(iter.max=10000))
OPTE2 = optim(initE,log.lik.E,control=list(maxit=10000))

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
CIE <- Conf.Int(log.lik.E,c(log(MLEE[1:5]),MLEE[-c(1:5)]),level=0.95)
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
true.haz <- Vectorize(function(t) hexpw(t,ae,be,ce))
fit1 <- Vectorize(function(t) hexpw(t, MLEEWG[1], MLEEWG[2], MLEEWG[3]))
fit2 <- Vectorize(function(t) hexpw(t, MLEC[2], MLEC[3], MLEC[4]))
fit3 <- Vectorize(function(t) hexpw(t, MLEE[3], MLEE[4], MLEE[5]))

curve(true.haz, 0, 5, lwd = 2, xlab = "Time", ylab ="Baseline Excess Hazard", main = "Fitted Excess Hazards",
      cex.axis = 1.5, cex.lab = 1.5, ylim= c(0.05,0.475), n = 1000)
curve(fit1, 0, 5, lwd = 2, lty = 2, col = "red",add = T, n = 1000)
curve(fit2, 0, 5, lwd = 2, lty = 2, col = "orange",add = T, n = 1000)
curve(fit3, 0, 5, lwd = 2, lty = 2, col = "blue",add = T, n = 1000)
legend(4,0.475, c("True","M1","M2","M3"),
       text.col = c("black","red","orange","blue"), col = c("black", "red", "orange", "blue"), lty = c(1, 2, 2, 2), lwd=2,
       merge = TRUE, bg = "white")
