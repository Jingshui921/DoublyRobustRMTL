rm(list = ls())
options(scipen=15, digits=15, max.print = 10000000)
source(' Functions.R')
# Libraries
library(MASS)
library(etm)
library(survival)
library(pseudo)
library(geepack)
library(tictoc)
library(parallel)
library(doParallel)
library(foreach)
# Total Sample Size and Iteration
n <- 250
k <- 20000
# Total number of covariates
p <- 5
# Parellism Computing
numCores <- detectCores()
cl <-makeCluster(15)
registerDoParallel(cl)
# Global Parameter
# 1. Multinormal Distribution
mu <- rep(0,p)
sigma <-  matrix(c(1,0.2,0,0,0,
                   0.2,1,0,0,0,
                   0,0,1,0,0,
                   0,0,0,1,0,
                   0,0,0,0,1), nrow = p)
# 2. True Propensity Score Model (D) for generating exposure
alpha <- c(log(3), log(1.5), 0, log(0.5),0)
names(alpha) <- c('X1','X2','V1','Z1','W1')
# 3. Sub-distribution hazard model for generating survival time
gamma <- 3
rho <- -1.5

beta1 <- c(log(0.7), log(0.6), log(0.5), 0, 0,log(0.5))
names(beta1) <- c('X1','X2','V1','Z1','W1','D')

beta2 <- c(log(0.6), log(0.5), log(0.4), 0, 0,log(0.85))
names(beta2) <- c('X1','X2','V1','Z1','W1','D')
# 4. Landmark time in RMTL and the delta in equally spaced event time.
tau <- 1 
delta <- 1/200
# 5. Setting Seed and the Seed for each iteration
set.seed(123)
seed.sample <- c(sample(1000000:9999999,size=k))[-(1:10000)] 
# 1:5000 used for n=2000 
# 5001:10000 used for n=1000 

Sys.time()
tic()
sim_ALL <- foreach (i=1:5000, .combine=rbind, .packages = c('MASS','etm','geepack','pseudo')) %dopar% {
  seed <- seed.sample[i]
  set.seed(seed)
  ########################################## Data Generation ##############################################
  Z <-  mvrnorm(n, mu = mu, Sigma = sigma, empirical = FALSE)
  colnames(Z) <- c('X1','X2','V1','Z1','W1') 
  rownames(Z) <- 1:n
  
  # True Propensity Score Model (D) for generating exposure
  true.propensity <- expit(Z %*% alpha)
  D <- rbinom(n,size = 1, prob = true.propensity)
  n0 <- table(D)['0'] # 979
  n1 <- table(D)['1'] # 1021
  Z <- cbind(Z,D=D)
  
  # Assign Event Type by the Marginal Distribution of Event Type 
  # Epsilon = 1 is the primary event, Epsilon = 2 is competing risk
  # Random error comes from both sampling variability and non-deterministic counter-factuals
  temp <- exp(Z %*% beta1)*gamma/rho 
  Epsilon <- 1 + rbinom(n,size = 1, prob = exp(temp)) 
  # table(Epsilon) 
  
  # Conditional Survival time of T|Epsilon
  u <- runif(n, 0, 1)
  Time <- (2-Epsilon)*(1/rho)*cloglog(x=u*(1-exp(temp)), a=1/temp,b=1)+ # Epsilon=1
    (Epsilon-1)*(-exp(-Z %*% beta2))*log(1-u)                           # Epsilon=2
  # max(Time)
  
  # Censoring
  u <- runif(n, 0.5, 3)
  Event <- ifelse(Time<=u,1,0)*Epsilon # Event = 0 if censored
  Obs.Time <- pmin(Time,u) # table(Event, D) # sum(Event==0)/n
  
  # Order the simulated dataset by exposure groups D
  dat <- data.frame(Z,Event,Obs.Time)
  dat$D <- as.factor(dat$D)
  dat <- dat[order(dat$D,decreasing = TRUE),]
  dat$ID <- 1:n
  
  true.propensity <- true.propensity[order(D,decreasing = TRUE)]
  
  # "Censored" all events after tau
  loc <- dat$Obs.Time>tau
  dat$Event[loc] <- 0
  dat$Obs.Time[loc] <- tau
  
  dat.D1 <- dat[dat$D=='1',]
  dat.D0 <- dat[dat$D=='0',]
  
  ########################################## Estimation ##############################################
  # Unadjusted Estimator of RMTL
  diff.rmtl.unadj <- diff(rmtl_unadj(dat$Obs.Time, dat$Event, group = dat$D, eoi = '1', tau = 1)[,'rmtl'])
  ## need to add the greenwood type formula
  
  # Pseudo-Observation of Aalen-Johansen Estimator
  tmax <- sort(unique(c(dat$Obs.Time[dat$Event!=0])))
  Pseudo.AJ.D1 <- Pseudo_AJ_cif_2(dat.D1$Obs.Time, dat.D1$Event, eoi = '1', tmax)
  Pseudo.AJ.D0 <- Pseudo_AJ_cif_2(dat.D0$Obs.Time, dat.D0$Event, eoi = '1', tmax)
  
  tequal <- seq(tmax[1], 1-delta, by = delta )
  Pseudo.AJ.D1.equal <- Pseudo_AJ_cif_2(dat.D1$Obs.Time, dat.D1$Event, eoi = '1', tequal)
  Pseudo.AJ.D0.equal <- Pseudo_AJ_cif_2(dat.D0$Obs.Time, dat.D0$Event, eoi = '1', tequal)
  
  # Pseudo-Observation of RMTL
  dat.D1$ps.rmtl <- Pseudo_rmtl(Pseudo.AJ.D1$time,Pseudo.AJ.D1$pseudo.cif, tau=1)
  dat.D0$ps.rmtl <- Pseudo_rmtl(Pseudo.AJ.D0$time,Pseudo.AJ.D0$pseudo.cif, tau=1)
  dat$ps.rmtl <- c(dat.D1$ps.rmtl,dat.D0$ps.rmtl)
  
  # IPW estimator using the true propensity score and pseudo-observations of RMTL
  diff.rmtl.ipw.alpha <- rmtl_ipw(pseudo.rmtl = dat$ps.rmtl, group = dat$D , propensity = true.propensity, 
                                  dmatrix_ps = dat[,c('X1','X2','Z1')])[3,2]
  
  # IPW Estimator of RMTL
  # summary(PS) # alpha #expit(data.matrix(dat[1:5,c('X1','X2','Z1')])%*%matrix(PS$coefficients,ncol=1))
  
  PS <- glm(D ~ X1+X2+Z1-1, data=dat,family = 'binomial')  
  dat$true.propensity <- PS$fitted.values
  rmtl.ipw <- rmtl_ipw(pseudo.rmtl = dat$ps.rmtl, 
                       group = dat$D , 
                       propensity = dat$true.propensity, 
                       dmatrix_ps = dat[,c('X1','X2','Z1')])
  diff.rmtl.ipw.true <- rmtl.ipw[3,2]
  var.diff.rmtl.ipw.true <- rmtl.ipw[4,2]
  
  PS <- glm(D ~ X1+V1+W1-1, data=dat,family = 'binomial') 
  dat$f1.propensity <- PS$fitted.values
  rmtl.ipw <- rmtl_ipw(pseudo.rmtl = dat$ps.rmtl, 
                       group = dat$D , 
                       propensity = dat$f1.propensity, 
                       dmatrix_ps = dat[,c('X1','V1','W1')])
  diff.rmtl.ipw.f1 <- rmtl.ipw[3,2]
  var.diff.rmtl.ipw.f1 <- rmtl.ipw[4,2]
  
  # Model-Based Estimator of RMTL
  b.equal <- cbind(dat[rep(1:n, each=length(tequal)),c('ID','X1','X2','V1','Z1','W1','D')],
                   po.cif = c(as.vector(Pseudo.AJ.D1.equal$pseudo.cif), as.vector(Pseudo.AJ.D0.equal$pseudo.cif)),
                   time = rep(tequal,times=n))
  b.equal$time <- as.factor(b.equal$time)
  
  fit <- geese(po.cif ~ time+X1+X2+V1+D-1,
               data=b.equal,
               id=ID, # A vector which identifies the clusters. Data are assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
               jack=T, # if approximate jackknife variance estimate should be computed.
               scale.fix = TRUE, # if the scale should be fixed.
               family = gaussian,
               mean.link = "cloglog",
               corstr = 'independence') 
  
  fg.rmtl <- rmtl_mb(beta = fit$beta, X = dat[,c('X1','X2','V1')], tequal, tau)
  dat$fg.rmtl.D0.true <- fg.rmtl$fg.po.rmtl.D0
  dat$fg.rmtl.D1.true <- fg.rmtl$fg.po.rmtl.D1
  diff.rmtl.mb.true <- as.numeric(diff(apply(fg.rmtl, MARGIN = 2, mean)))
  
  fit <- geese(po.cif ~ time+X1+Z1+W1+D-1,
               data=b.equal,
               id=ID, # A vector which identifies the clusters. Data are assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
               jack=T, # if approximate jackknife variance estimate should be computed.
               scale.fix = TRUE, # if the scale should be fixed.
               family = gaussian,
               mean.link = "cloglog",
               corstr = 'independence') 
  
  fg.rmtl <- rmtl_mb(beta = fit$beta, X = dat[,c('X1','Z1','W1')], tequal, tau)
  dat$fg.rmtl.D0.f1 <- fg.rmtl$fg.po.rmtl.D0
  dat$fg.rmtl.D1.f1 <- fg.rmtl$fg.po.rmtl.D1
  diff.rmtl.mb.f1 <- as.numeric(diff(apply(fg.rmtl, MARGIN = 2, mean)))
  
  # Double Robust Estimator
  db.rmtl <- rmtl_db(group = dat$D, 
                     ps.rmtl = dat$ps.rmtl, 
                     mb.rmtl.d1 = dat$fg.rmtl.D1.true, 
                     mb.rmtl.d0 = dat$fg.rmtl.D0.true, 
                     propensity = dat$true.propensity)
  diff.rmtl.db.FGtrue.PStrue <- as.numeric(db.rmtl[3])
  var.diff.rmtl.db.FGtrue.PStrue <- as.numeric(db.rmtl[4])
  
  db.rmtl <- rmtl_db(group = dat$D, 
                     ps.rmtl = dat$ps.rmtl, 
                     mb.rmtl.d1 = dat$fg.rmtl.D1.true, 
                     mb.rmtl.d0 = dat$fg.rmtl.D0.true, 
                     propensity = dat$f1.propensity)
  diff.rmtl.db.FGtrue.PSf1 <- as.numeric(db.rmtl[3])
  var.diff.rmtl.db.FGtrue.PSf1 <- as.numeric(db.rmtl[4])
  
  db.rmtl <- rmtl_db(group = dat$D, 
                     ps.rmtl = dat$ps.rmtl, 
                     mb.rmtl.d1 = dat$fg.rmtl.D1.f1, 
                     mb.rmtl.d0 = dat$fg.rmtl.D0.f1, 
                     propensity = dat$true.propensity)
  diff.rmtl.db.FGf1.PStrue <- as.numeric(db.rmtl[3])
  var.diff.rmtl.db.FGf1.PStrue <- as.numeric(db.rmtl[4])

  db.rmtl <- rmtl_db(group = dat$D, 
                     ps.rmtl = dat$ps.rmtl, 
                     mb.rmtl.d1 = dat$fg.rmtl.D1.f1, 
                     mb.rmtl.d0 = dat$fg.rmtl.D0.f1, 
                     propensity = dat$f1.propensity)
  diff.rmtl.db.FGf1.PSf1 <- as.numeric(db.rmtl[3])
  var.diff.rmtl.db.FGf1.PSf1 <- as.numeric(db.rmtl[4])

  # Checking Difference
  dat.D1 <- dat[dat$D=='1',]
  dat.D0 <- dat[dat$D=='0',]
  
  diff.db.ipw.FGtrue.PStrue <- diff.db.ipw(as.numeric(dat$D)-1, dat$true.propensity, dat$fg.rmtl.D1.true, dat$fg.rmtl.D0.true)
  diff.db.ipw.FGtrue.PSf1   <- diff.db.ipw(as.numeric(dat$D)-1, dat$f1.propensity, dat$fg.rmtl.D1.true, dat$fg.rmtl.D0.true)
  diff.db.ipw.FGf1.PStrue   <- diff.db.ipw(as.numeric(dat$D)-1, dat$true.propensity, dat$fg.rmtl.D1.f1, dat$fg.rmtl.D0.f1)
  diff.db.mb.FGtrue.PStrue  <- diff.db.mb(dat.D1$ps.rmtl, dat.D0$ps.rmtl,dat.D1$fg.rmtl.D1.true, dat.D0$fg.rmtl.D0.true, dat.D1$true.propensity, dat.D0$true.propensity)
  diff.db.mb.FGtrue.PSf1    <- diff.db.mb(dat.D1$ps.rmtl, dat.D0$ps.rmtl,dat.D1$fg.rmtl.D1.true, dat.D0$fg.rmtl.D0.true, dat.D1$f1.propensity, dat.D0$f1.propensity)
  diff.db.mb.FGf1.PStrue    <- diff.db.mb(dat.D1$ps.rmtl, dat.D0$ps.rmtl,dat.D1$fg.rmtl.D1.f1, dat.D0$fg.rmtl.D0.f1, dat.D1$true.propensity, dat.D0$true.propensity) 
  diff.db.both.FGtrue.PStrue <- diff.db.both(dat.D1$ps.rmtl, dat.D0$ps.rmtl, dat.D1$fg.rmtl.D1.true, dat.D0$fg.rmtl.D0.true, dat.D1$true.propensity, dat.D0$true.propensity)
  diff.db.both.FGtrue.PSf1   <- diff.db.both(dat.D1$ps.rmtl, dat.D0$ps.rmtl, dat.D1$fg.rmtl.D1.true, dat.D0$fg.rmtl.D0.true, dat.D1$f1.propensity, dat.D0$f1.propensity)
  diff.db.both.FGf1.PStrue   <- diff.db.both(dat.D1$ps.rmtl, dat.D0$ps.rmtl, dat.D1$fg.rmtl.D1.f1, dat.D0$fg.rmtl.D0.f1, dat.D1$true.propensity, dat.D0$true.propensity) 
  diff.db.both.FGf1.PSf1     <- diff.db.both(dat.D1$ps.rmtl, dat.D0$ps.rmtl, dat.D1$fg.rmtl.D1.f1, dat.D0$fg.rmtl.D0.f1, dat.D1$f1.propensity, dat.D0$f1.propensity)

  return(c(i,
           seed,
           diff.rmtl.unadj,   
           diff.rmtl.ipw.alpha,
           diff.rmtl.ipw.true,
           diff.rmtl.ipw.f1,
           diff.rmtl.mb.true,
           diff.rmtl.mb.f1,
           diff.rmtl.db.FGtrue.PStrue,
           diff.rmtl.db.FGtrue.PSf1,
           diff.rmtl.db.FGf1.PStrue,
           diff.rmtl.db.FGf1.PSf1,
           var.diff.rmtl.ipw.true,
           var.diff.rmtl.ipw.f1,
           var.diff.rmtl.db.FGtrue.PStrue,
           var.diff.rmtl.db.FGtrue.PSf1,
           var.diff.rmtl.db.FGf1.PStrue,
           var.diff.rmtl.db.FGf1.PSf1,
           diff.db.ipw.FGtrue.PStrue, 
           diff.db.ipw.FGtrue.PSf1,   
           diff.db.ipw.FGf1.PStrue,   
           diff.db.mb.FGtrue.PStrue,  
           diff.db.mb.FGtrue.PSf1,    
           diff.db.mb.FGf1.PStrue,    
           diff.db.both.FGtrue.PStrue,
           diff.db.both.FGtrue.PSf1,  
           diff.db.both.FGf1.PStrue,  
           diff.db.both.FGf1.PSf1
  ))
}
toc()

colnames(sim_ALL) <- c('Iteration','Seed',
                       'diff.rmtl.unadj',   
                       'diff.rmtl.ipw.alpha',
                       'diff.rmtl.ipw.true', 
                       'diff.rmtl.ipw.f1', 
                       'diff.rmtl.mb.true', 
                       'diff.rmtl.mb.f1',
                       'diff.rmtl.db.FGtrue.PStrue',
                       'diff.rmtl.db.FGtrue.PSf1',
                       'diff.rmtl.db.FGf1.PStrue',
                       'diff.rmtl.db.FGf1.PSf1',
                       'var.diff.rmtl.ipw.true',
                       'var.diff.rmtl.ipw.f1',
                       'var.diff.rmtl.db.FGtrue.PStrue',
                       'var.diff.rmtl.db.FGtrue.PSf1',
                       'var.diff.rmtl.db.FGf1.PStrue',
                       'var.diff.rmtl.db.FGf1.PSf1',
                       'diff.db.ipw.FGtrue.PStrue', 
                       'diff.db.ipw.FGtrue.PSf1',   
                       'diff.db.ipw.FGf1.PStrue',   
                       'diff.db.mb.FGtrue.PStrue',
                       'diff.db.mb.FGtrue.PSf1',  
                       'diff.db.mb.FGf1.PStrue',  
                       'diff.db.both.FGtrue.PStrue',
                       'diff.db.both.FGtrue.PSf1',  
                       'diff.db.both.FGf1.PStrue',  
                       'diff.db.both.FGf1.PSf1'
)

write.table(data.frame(sim_ALL), 
            file = "ALLmodel_n250.csv", 
            sep = ',',
            append = TRUE,
            col.names = !file.exists("ALLmodel_n250.csv"),
            row.names = FALSE,
            quote=FALSE)

