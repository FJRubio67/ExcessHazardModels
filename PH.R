###########################################################################################################################
# R code associated to the paper 
# "On a general structure for hazard-based regression models: an application to population-based cancer research"
# See also: http://rpubs.com/FJRubio/GHPH
#           http://rpubs.com/FJRubio/GHGH
#           http://rpubs.com/FJRubio/EWD
###########################################################################################################################


rm(list=ls())
# library to compile the functions
library(compiler)
enableJIT(3)

library(knitr)

# True values of the parameters
ae <- 1.75; be <- 0.5; ce <- 2.5; beta <- c(0.035,0.2,0.3); rate0 <- 0; betaH = c(0,0,0)
true.par <- c(ae,be,ce,betaH,beta)
true.parc <- c(ae,be,ce,beta)

# read the data (available at: https://sites.google.com/site/fjavierrubio67/dataPH.txt)
mydata <- read.table("dataPH.txt")
colnames(mydata)

# create the variables for the models
x <- as.matrix(mydata[,1:3]) # covariates
sim <- mydata[,4] # survival time
ind.cens <- mydata[,5] # vital status
hp.as <- mydata[,6] # Population (expected) mortality hazard

# Summary of the data
summary(sim) # simulated survival times
apply(x,2,summary)  # covariates
mean(ind.cens) # censoring rate

##############################################################################################################
# Hazard functions
##############################################################################################################

# Weibull Hazard
hw <- function(t,lambda,kappa){
  pdf0 <-  dweibull(t,scale=lambda,shape=kappa)
  cdf0 <- pweibull(t,scale=lambda,shape=kappa)
  return(pdf0/(1-cdf0))
}                                                                                      

# Weibull Cumulative hazard
Hw <- function(t,lambda,kappa){
  cdf <- pweibull(t,scale=lambda,shape=kappa)  
  return(-log(1-cdf))
}

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


###################################################################################################################
# Estimation
###################################################################################################################
# Implemented using nlminb(), alternatively, use optim()
#--------------------------------------------------------------------------------------------------------
# Proportional Hazards
#--------------------------------------------------------------------------------------------------------

#----------------------------------------------
# log-likelihood: Weibull model
#----------------------------------------------

log.likw1 <- cmpfun(function(par){
  if(anyNA(par)) return(Inf)
  ae0 <- par[1]; be0 <- par[2]; beta0 <- par[3:(3+dim(x)[2]-1)];
  if(par[1]>0 & par[2]>0 ){
    exp.x.beta <- exp(x%*%beta0)
    haz0 <- hp.as + hw(sim,ae0,be0)*exp.x.beta
    val <- - ind.cens*log(haz0) + Hw(sim,ae0,be0)*exp.x.beta
    return(sum(val))
  }
  else return(Inf)
})
# Alternatively use optim(c(ae,be,beta),log.likw1,control=list(maxit=10000))
OPTW1 <- nlminb(c(ae,be,beta),log.likw1,control=list(iter.max=10000))
MLEW1 <- OPTW1$par
AICW1 <- 2*OPTW1$objective + 2*length(MLEW1)

#----------------------------------------------
# log-likelihood: EW Model
#----------------------------------------------
log.likew1 <- cmpfun(function(par){
  if(anyNA(par)) return(Inf)
  ae0 <- par[1]; be0 <- par[2]; ce0 <- par[3]; beta0 <- par[4:(4+dim(x)[2]-1)];
  if(par[1]>0 & par[2]>0 & par[3]>0){
    exp.x.beta <- exp(x%*%beta0)
    haz0 <- hp.as + hexpw(sim,ae0,be0,ce0)*exp.x.beta
    CH0 <- Hexpw(sim,ae0,be0,ce0)*exp.x.beta
    if(anyNA(CH0)) return(Inf)
    if(anyNA(haz0)) return(Inf)
    else  val <- - ind.cens*log(haz0) + CH0
    return(sum(val))
  }
  else return(Inf)
})

OPTEW1 <- nlminb(c(ae,be,ce,beta),log.likew1,control=list(iter.max=10000))
MLEEW1 <- OPTEW1$par
AICEW1 <- 2*OPTEW1$objective + 2*length(MLEEW1)

#--------------------------------------------------------------------------------------------------------
# Accelerated Hazards
#--------------------------------------------------------------------------------------------------------

#----------------------------------------------
# log-likelihood: EW Model
#----------------------------------------------
log.likew2 <- cmpfun(function(par){
  if(anyNA(par)) return(Inf)
  ae0 <- par[1]; be0 <- par[2]; ce0 <- par[3]; beta0 <- par[4:(4+dim(x)[2]-1)];
  if(par[1]>0 & par[2]>0 & par[3]>0){
    exp.x.beta <- exp(x%*%beta0)
    haz0 <- hp.as + hexpw(sim*exp.x.beta,ae0,be0,ce0)
    CH0 <- Hexpw(sim*exp.x.beta,ae0,be0,ce0)/exp.x.beta
    if(anyNA(CH0)) return(Inf)
    if(anyNA(haz0)) return(Inf)
    else  val <- - ind.cens*log(haz0) + CH0
    return(sum(val))
  }
  else return(Inf)
})

OPTEW2 <- nlminb(c(ae,be,ce,beta),log.likew2,control=list(iter.max=10000))
MLEEW2 <- OPTEW2$par
AICEW2 <- 2*OPTEW2$objective + 2*length(MLEEW2)

#--------------------------------------------------------------------------------------------------------
# Accelerated Failure Time
#--------------------------------------------------------------------------------------------------------

#----------------------------------------------
# log-likelihood: EW Model
#----------------------------------------------
log.likew3 <- cmpfun(function(par){
  if(anyNA(par)) return(Inf)
  ae0 <- par[1]; be0 <- par[2]; ce0=par[3]; beta0 <- par[4:(4+dim(x)[2]-1)];
  if(par[1]>0 & par[2]>0 & par[3]>0){
    exp.x.beta <- exp(x%*%beta0)
    haz0 <- hp.as + hexpw(sim*exp.x.beta,ae0,be0,ce0)*exp.x.beta
    CH0 <- Hexpw(sim*exp.x.beta,ae0,be0,ce0)
    if(anyNA(CH0)) return(Inf)
    if(anyNA(haz0)) return(Inf)
    else  val <- - ind.cens*log(haz0) + CH0
    return(sum(val))
  }
  else return(Inf)
})

OPTEW3 <- nlminb(c(ae,be,ce,beta),log.likew3,control=list(iter.max=10000))
MLEEW3 <- OPTEW3$par
AICEW3 <- 2*OPTEW3$objective + 2*length(MLEEW3)

#--------------------------------------------------------------------------------------------------------
# General Model
#--------------------------------------------------------------------------------------------------------

#----------------------------------------------
# log-likelihood: EW Model
#----------------------------------------------
log.likewG <- cmpfun(function(par){
  if(anyNA(par)) return(Inf)
  ae0 <- par[1]; be0 <- par[2]; ce0=par[3]; beta0 <- par[4:(4+dim(x)[2]-1)]; beta1 <- tail(par,n = dim(x)[2]);
  if(par[1]>0 & par[2]>0 & par[3]>0){
    exp.x.beta0 <- exp(x%*%beta0)
    exp.x.beta1 <- exp(x%*%beta1)
    exp.x.dif <- exp(x%*%(beta1-beta0))
    haz0 <- hp.as + hexpw(sim*exp.x.beta0,ae0,be0,ce0)*exp.x.beta1
    CH0 <- Hexpw(sim*exp.x.beta0,ae0,be0,ce0)*exp.x.dif
    if(anyNA(CH0)) return(Inf)
    if(anyNA(haz0)) return(Inf)
    else  val <- - ind.cens*log(haz0) + CH0
    return(sum(val))
  }
  else return(Inf)
})

OPTEWG <- nlminb(c(ae,be,ce,betaH,beta), log.likewG,control=list(iter.max=10000))
MLEEWG <- OPTEWG$par
AICEWG <- 2*OPTEWG$objective + 2*length(MLEEWG)


###################################################################################################################
# Model Comparison
###################################################################################################################

AIC <- c(AICW1,AICEW1,AICEW2,AICEW3,AICEWG)
AIC.min <- which(AIC==min(AIC))
model.names <- c("W-PH","EW-PH","EW-AH","EW-AFT","EW-General")
names(AIC) <- model.names
AIC <- as.matrix(AIC);  colnames(AIC) <- c("AIC")
print(kable(AIC),digits=4)
model.names[order(AIC)]
cat("Best Model based on AIC = ",model.names[AIC.min])

# True parameter values vs MLE of best model
print(kable(cbind(true.parc,MLEEW1),digits=4))



### Derivation of the confidence intervals

###################################################################################################################
# Confidence intervals for best model
###################################################################################################################

# Required package
library(numDeriv)

###########################################################################################
# Function to calculate the normal confidence intervals
# The parameters indicated with "index" are transformed to the real line using log()
###########################################################################################
# FUN   : minus log-likelihood function to be used to calculate the confidence intervals
# MLE   : maximum likelihood estimator of the parameters of interest
# level : confidence level
# index : position of the positive parameters under the original parameterisation

Conf.Int <- function(FUN,MLE,level=0.95,index=NULL){
  sd.int <- abs(qnorm(0.5*(1-level)))
  tempf <- function(par){
    par[index] = exp(par[index])
    return(FUN( par ))
  }
  r.MLE <- MLE
  r.MLE[index] <- log(MLE[index])
  HESS <- hessian(tempf,x=r.MLE)
  Fisher.Info <- solve(HESS)
  Sigma <- sqrt(diag(Fisher.Info))
  U<- r.MLE + sd.int*Sigma
  L<- r.MLE - sd.int*Sigma
  C.I <- cbind(L,U,r.MLE, Sigma)
  names.row <- paste0("par", seq_along(1:length(MLE)))
  names.row[index] <- paste0("log.par", seq_along(index))
  rownames(C.I)<- names.row
  colnames(C.I)<- c("Lower","Upper","Transf MLE", "Std. Error")
  return(C.I)
}

CI <- Conf.Int(log.likew1,MLEEW1,level=0.95,index=c(1,2,3))
print(kable(CI,digits=4))



### Results: Graphical comparisons and Net Survival estimates derived from the selected model


###################################################################################################################
# Fitted hazard vs True hazard
###################################################################################################################


hazard.fit <- Vectorize(function(t)  hexpw(t,MLEEW1[1],MLEEW1[2],MLEEW1[3]))
hazard.true <- Vectorize(function(t)  hexpw(t,ae,be,ce))

curve(hazard.true,0,5,lwd=3,n=10000,xlab="time",ylab="Excess Hazard",cex.axis=1.5,cex.lab=1.5,ylim=c(0,0.26), col="black")
curve(hazard.fit,0,5,lwd=3,n=10000,lty=2,col="black",add=T)
legend(3.5,0.24, c("True","Fitted"),
       text.col = c("black","black"), col = c("black","black"), lty = c(1, 2),lwd=3,
       merge = TRUE, bg = "gray90")


###################################################################################################################
# Fitted Net survival for a single patient with sex = 0, age = 70, Comorbidity = 0, 1 
###################################################################################################################

NS0 <- Vectorize(function(t)  exp(-Hexpw(t,MLEEW1[1],MLEEW1[2],MLEEW1[3])))
NS1 <- Vectorize(function(t)  exp(-Hexpw(t,MLEEW1[1],MLEEW1[2],MLEEW1[3])*exp(MLEEW1[6])))

curve(NS0,0,5,lwd=3,n=10000,xlab="time",ylab="Net Survival",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1), col="blue")
curve(NS1,0,5,lwd=3,n=10000,lty=2,col="red",add=T)
legend(3.5,0.9, c("Comorb = 0","Comorb = 1"),
       text.col = c("blue","red"), col = c("blue","red"), lty = c(1, 2),lwd=3,
       merge = TRUE, bg = "gray90")

########################################################################################################
# Confidence intervals for the net survival based on Bootstrap
########################################################################################################
# Required package
library(mvtnorm)

# t0 : time at what the confidence interval will be calculated
# level : confidence level
# n.boot : number of bootstrap iterations
# x.int : subset of the design matrix associated to the group of interest

# The variance of the normal approximation to the distribution of the MLE
tempf <- function(par){
  par[1:3] = exp(par[1:3])
  return(log.likew1( par ))
}
r.MLE <- MLEEW1
r.MLE[1:3] <- log(MLEEW1[1:3])
HESS <- hessian(tempf,x=r.MLE)
Sigma <- solve(HESS)

# The function 
conf.int.NS <- function(t0,level,n.boot,x.int){
  boot <- vector()
  S.par <- function(par) mean( exp(-Hexpw(t0,par[1],par[2],par[3])*exp(x.int%*%par[-c(1:3)])))
  
  for(i in 1:n.boot) {
    val <- rmvnorm(1,mean = r.MLE, sigma = Sigma)
    val[1:3] <- exp(val[1:3])
    boot[i] <- S.par(val)
  }
  
  L <- quantile(boot,(1-level)*0.5)
  U <- quantile(boot,(1+level)*0.5)
  
  M <- S.par(MLEEW1)
  
  return(c(L,M,U))
}

#-----------------------------------
# Net Survival for the whole Cohort
# at Comorbidity = 0 , 1 
# at years 1, 2, 3, 4, 4.9
#-----------------------------------

# times
times <- c(1,2,3,4,4.9)

CIS0 <- CIS1 <- matrix(0, ncol = 4, nrow = length(times))

for(k in 1:length(times)) CIS0[k,] <- c(times[k],conf.int.NS(times[k],0.95,10000,x[which(x[,3]==0),]))
for(k in 1:length(times)) CIS1[k,] <- c(times[k],conf.int.NS(times[k],0.95,10000,x[which(x[,3]==1),]))

colnames(CIS0) <- cbind("year","lower","net_survival","upper")
print(kable(CIS0,digits=4))

colnames(CIS1) <- cbind("year","lower","net_survival","upper")
print(kable(CIS1,digits=4))

#--------------------------------------------
# Net Survival for an age group: 60 - 70
# at Comorbidity = 0 , 1 
# sex = 1
#--------------------------------------------

CISA0 <- CISA1 <- matrix(0, ncol = 4, nrow = length(times))

for(k in 1:length(times)) CISA0[k,] <- c(times[k],conf.int.NS(times[k],0.95,10000,x[which(x[,3]==0 & x[,1]>=-10 & x[,1]<=0 & x[,2]==1),]))
for(k in 1:length(times)) CISA1[k,] <- c(times[k],conf.int.NS(times[k],0.95,10000,x[which(x[,3]==1 & x[,1]>=-10 & x[,1]<=0 & x[,2]==1),]))

colnames(CISA0) <- cbind("year","lower","net_survival","upper")
print(kable(CISA0,digits=4))

colnames(CISA1) <- cbind("year","lower","net_survival","upper")
print(kable(CISA1,digits=4))
