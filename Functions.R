# Functions
expit <- function(x){1/(1+exp(-x))}
cloglog <- function(x,a=1,b=0){log(b-a*log(1-x))} # 0<x<1,  y>0
gompertz <- function(t, gamma, rho) {gamma*exp(rho*t)} # gamma>0, rho<0, t>=
int.step <- function(x,y,tau=1) sum(diff(c(x,tau))*y) # length(x)=length(y), max(x) <= tau

# Joint (T, Epsilon) CIF1, CIF2
cif1 <- function(t,rho=-1.5, gamma=1.5, Z, beta1){
  t <- matrix(t,nrow=1)
  temp <- exp(Z %*% beta1)*gamma/rho
  return(1-exp(temp%*%(1-exp(rho*t))))
}
cif2 <-function(t,rho=-1.5, gamma=1.5, Z, beta1,beta2){
  t <- matrix(t,nrow=1)
  temp <- exp(Z %*% beta1)*gamma/rho
  return(exp(temp)%*%(1-exp(-exp(Z %*% beta2)%*%t)))
}

# Conditional CIF1 given event=1 (T|Epsilon=1)
cif1.e1 <- function(t,rho, gamma, Z, beta1){  
  t <- matrix(t,nrow=1)
  temp <- exp(Z %*% beta1)*gamma/rho
  return(diag(as.vector(1/(1-exp(temp))),nrow=nrow(Z))%*%(1-exp(temp%*%(1-exp(rho*t)))))
} 

# Conditional CIF2 given event=2 (T|Epsilon=2)
cif2.e2 <- function(t,Z,beta2){
  t <- matrix(t,nrow=1)
  return(1-exp(-exp(Z %*% beta2)%*%t))
}  


cif <- function(times, event, weight=NULL, tau){
  # Author: Sarah
  # Input: event times, event indicator (1,2,0-censored), landmark tau in RMTL
  # Output: When weight=NULL, the output is the unadjusted Aalen-Johansen estimator of CIF for both events.
  if(is.null(weight)){weight <- rep(1, length(times))}
  entry=rep(0, length(times))
  
  data <-  data.frame(entry, times, event, weight)
  data$event[data$times>tau] <- 0
  data$times[data$times>tau] <- tau
  
  tj <- data$times[data$event!=0]
  tj <- unique(tj[order(tj)])
  num.tj <- length(tj)
  
  num.atrisk <- sapply(tj, function(tj) sum(data$weight[data$entry<tj & data$times>=tj]))
  num.ev1 <- sapply(tj, function(tj) sum(data$weight[data$event==1 & data$times==tj]))
  num.ev2 <- sapply(tj, function(tj) sum(data$weight[data$event==2 & data$times==tj]))
  num.ev <- num.ev1 + num.ev2 
  
  h1 <- num.ev1/num.atrisk
  h2 <- num.ev2/num.atrisk
  h <- num.ev/num.atrisk
  
  s <- cumprod(c(1, 1 - h))
  s <- s[1:length(s)-1]
  
  theta1 <- s * h1
  theta2 <- s * h2
  
  cif1 <- cumsum(theta1)
  cif2 <- cumsum(theta2)
  
  #cr1 <- 1 - cumprod(c(1, 1 - h1[1:num.tj-1])) # 1-KM of event 1
  #cr2 <- 1 - cumprod(c(1, 1 - h2[1:num.tj-1])) # 1-KM of event 2
  
  return(list(tj=tj, cif1=cif1, cif2=cif2, num.atrisk=num.atrisk))
}

rmtl1 <- function(Z, D=NULL, beta1, rho=-1.5, gamma=1.5, tau=1, delta=0.01){
  # Input: Covariates Z, coefficients beta1 in CIF1, group assignment D
  # Output: the true RMTL of event 1 of each individual assuming all subject are in group D.
  # Exposure variable D must be the last column in Design matrix D and coefficient vector beta1
  
  if(is.null(D)) {Z.D <- Z}
  else{Z.D <- cbind(matrix(Z[,1:(length(beta1)-1)],nrow = nrow(Z)),rep(D,nrow(Z)))}
  
  temp <- exp(Z.D %*% beta1)*gamma/rho
  t <- matrix(seq(0,tau-delta,delta),nrow = 1)
  rmtl.D <- tau-exp(temp)*(apply(exp(-temp%*%exp(rho*t)), MARGIN=1, FUN=sum)*delta) # dim(exp(-temp%*%exp(rho*t))) = n*length(t)
  return(rmtl.D)
}

diff.rmtl <- function(Z, beta1,seed,seed0, rho=-1.5, gamma=1.5, tau=1){
  # Input: Covariates Z and coefficients beta1 in CIF1
  # Output: The RMTL of event 1 assuming all subject are exposed(D=1) and unexposed(D=0), respectively
  
  Z.D1 <- cbind(Z[,-ncol(Z)],D=rep(1,nrow(Z)))
  Z.D0 <- cbind(Z[,-ncol(Z)],D=rep(0,nrow(Z)))
  
  set.seed(seed)
  Epsilon.D1 <- 1 + rbinom(n,size = 1, prob = exp(exp(Z.D1 %*% beta1)*gamma/rho))
  Epsilon.D0 <- 1 + rbinom(n,size = 1, prob = exp(exp(Z.D0 %*% beta1)*gamma/rho))
  
  set.seed(seed0)
  u <- runif(n, 0, 1)
  
  temp <- exp(Z.D1 %*% beta1)*gamma/rho
  Time.D1 <- (2-Epsilon.D1)*(1/rho)*cloglog(x=u*(1-exp(temp)), a=1/temp,b=1)+  
    (Epsilon.D1-1)*(-exp(-Z.D1 %*% beta2))*log(1-u)                  
  temp <- exp(Z.D0 %*% beta1)*gamma/rho
  Time.D0 <- (2-Epsilon.D0)*(1/rho)*cloglog(x=u*(1-exp(temp)), a=1/temp,b=1)+   
    (Epsilon.D0-1)*(-exp(-Z.D0 %*% beta2))*log(1-u)                  
  
  cif.d1 <- cif(Time.D1, Epsilon.D1, tau=1)
  cif.d0 <- cif(Time.D0, Epsilon.D0, tau=1)
  
  #diff.rmtl1.true <- int.step(x=cif.d1$tj,y=cif.d1$cif1,tau=tau) - int.step(x=cif.d0$tj,y=cif.d0$cif1,tau=tau) 
  return(data.frame(rmtl.d1=int.step(x=cif.d1$tj,y=cif.d1$cif1,tau=tau),
                    rmtl.d0=int.step(x=cif.d0$tj,y=cif.d0$cif1,tau=tau) ))
}

true_coef <- function(tmax,rho, gamma){
  log(gamma*(exp(rho*tmax)-1)/rho)
}


########################### Estimator #############################
AJ_cif <- function(time, event, eoi, tmax=NULL){
  
  # if (missing(time))
  #   stop('Please specify the failure time as a numeric vector')
  # if (any(is.na(time)))
  #   stop("Missingness is not allowed in the time vector")
  # if (any(time<=0))
  #   stop('Failure time can only be positive.')
  # 
  # if (missing(event))
  #   stop('Please specify the observed event types as an integer factor')
  # if (any(is.na(event)))
  #   stop("Missingness is not allowed in the event")
  # 
  # if (missing(eoi))
  #   stop('Please specify the event of interest (eoi) as an unique value in the event.')
  # if (!eoi %in% unique(event))
  #   stop('eoi (event of interest) should be observed in event.')
  # if (!is.null(tmax) & class(tmax)!='numeric')
  #   stop('Please specify the time points as a non-negative numeric vector.')
  
  event <- as.vector(event)
  data <- data.frame(time=time,event=event)
  
  evtime <- data$time[data$event!='0']
  evtime <- sort(unique(evtime))
  
  n.atrisk <- sapply(evtime, FUN = function(x) sum(data$time>=x)) 
  n.eoi <- sapply(evtime, FUN = function(x) sum(data$time==x & data$event==eoi)) # the increment of observed event of interest
  n.event <- sapply(evtime, FUN = function(x) sum(data$time==x & data$event!='0') ) # the increment of observed event
  
  nelson.aalen <- cumsum(n.eoi/n.atrisk)
  # var.nel.aal <- cumsum(n.eoi/n.atrisk/n.atrisk)
  kaplen.meier <- cumprod(c(1,1-n.event/n.atrisk))
  # var.kap.mei <- ((kaplen.meier[-1])^2)*cumsum(n.event/(n.atrisk*(n.atrisk-n.event)))
  aalen.johansen <- cumsum(kaplen.meier[-length(kaplen.meier)]*n.eoi/n.atrisk)
  # delta.hazard <- n.event/n.atrisk
  # var.aal.joh.greenwood <- 
  #   sapply(1:length(evtime), function(t){
  #     sum(
  #       ((aalen.johansen[t]-aalen.johansen[1:t])^2)*delta.hazard[1:t]/(n.atrisk[1:t]-n.event[1:t])+
  #         (kaplen.meier[1:t]^2*n.eoi[1:t]/n.atrisk[1:t]^3)*(n.atrisk[1:t]-n.eoi[1:t]-2*(n.atrisk[1:t]-n.event[1:t])*(aalen.johansen[1:t]-aalen.johansen[1:t])/kaplen.meier[(1:t)+1])
  #     )
  #   })
  
  # data.frame(time=c(0,evtime),
  #            n.atrisk=c(length(time),n.atrisk),
  #            n.eoi=c(0,n.eoi),
  #            n.event=c(0,n.event),
  #            nelson.aalen=c(0,nelson.aalen),
  #            kaplen.meier,
  #            aalen.johansen=c(0,aalen.johansen),
  #            delta.hazard = c(0,delta.hazard) )
  if (is.null(tmax)) {
    return(data.frame(evtime=evtime, cif=aalen.johansen, eoi=n.eoi))
  } else {
    tmax <- sort(tmax) # Test tmax <- c(0, evtime, 4)
    a <- sum(tmax<evtime[1])
    where <- sapply(tmax[(a+1):length(tmax)], function (x) max(which(evtime<=x)))# location of the last TRUE value
    aalen.johansen <- c(rep(0,a),aalen.johansen[where])  # AJ = 0 before the first event
    return(data.frame(time=tmax, cif=aalen.johansen))
  }
}

Pseudo_AJ_cif_1 <- function(time, event, eoi, tmax=NULL){
  
  if (missing(time))
    stop('Please specify the failure time as a numeric vector')
  if (any(is.na(time)))
    stop("Missingness is not allowed in the time vector")
  if (any(time<=0))
    stop('Failure time can only be positive.')
  
  if (missing(event))
    stop('Please specify the observed event types as an integer vector')
  if (any(is.na(event)))
    stop("Missingness is not allowed in the event vector")
  
  if (missing(eoi))
    stop('Please specify the event of interest (eoi) as an unique value in event.')
  if (!eoi %in% unique(event))
    stop('eoi (event of interest) should be observed in event.')
  
  if (!is.null(tmax) & class(tmax)!='numeric')
    stop('Please specify the time points as a numeric vector.')
  
  event <- as.vector(event)
  data <- data.frame(time=time,event=event)
  n <- length(time)
  
  evtime <- data$time[data$event!='0'] # need all event time to estimate kaplen.meier
  evtime <- sort(unique(evtime))
  
  n.atrisk <- sapply(evtime, FUN = function(x) sum(data$time>=x)) 
  n.eoi <- sapply(evtime, FUN = function(x) sum(data$time==x & data$event==eoi)) # the increment of observed event of interest
  n.event <- sapply(evtime, FUN = function(x) sum(data$time==x & data$event!='0') ) # the increment of observed event
  
  nelson.aalen <- cumsum(n.eoi/n.atrisk)
  kaplen.meier <- cumprod(c(1,1-n.event/n.atrisk)) 
  aalen.johansen <- cumsum(kaplen.meier[-length(kaplen.meier)]*n.eoi/n.atrisk)
  
  aalen.johansen.woi <- sapply(1:n, FUN = function(i) AJ_cif(time = time[-i],event = event[-i], eoi = '1', tmax = evtime)[2])
  aalen.johansen.woi <- matrix(unlist(aalen.johansen.woi),ncol=n, byrow=FALSE)
  aalen.johansen.pseudo <- sapply(1:n, FUN = function(j) n*aalen.johansen - (n-1)*aalen.johansen.woi[,j])
  colnames(aalen.johansen.pseudo) <- paste('n', 1:n, sep = '')
  
  if (is.null(tmax)) {
    rownames(aalen.johansen.pseudo) <- paste('t', 1:length(evtime), sep = '')
    return(list(time = evtime, pseudo.cif = aalen.johansen.pseudo))
  } else {
    tmax <- sort(tmax) # Test tmax <- c(0,evtime,4)
    a <- sum(tmax<evtime[1])
    where <- sapply(tmax[(a+1):length(tmax)], function (x) max(which(evtime<=x)))# location of the last TRUE value
    aalen.johansen.pseudo <- rbind(matrix(0, nrow = a, ncol = n), 
                                   aalen.johansen.pseudo[where,])
    rownames(aalen.johansen.pseudo) <- paste('t', 1:length(tmax), sep = '')
    return(list(time = tmax, pseudo.cif = aalen.johansen.pseudo))
  }
}

Pseudo_AJ_cif_2 <- function(time, event, eoi, tmax=NULL) {
  if (missing(time))
    stop('Please specify the failure time as a numeric vector')
  if (any(is.na(time)))
    stop("Missingness is not allowed in the time vector")
  if (any(time<=0))
    stop('Failure time can only be positive.')
  
  if (missing(event))
    stop('Please specify the observed event types as an integer vector')
  if (any(is.na(event)))
    stop("Missingness is not allowed in the event vector")
  if (class(event)!='numeric')
    stop('Event should be an integer vector where censoring is coded as 0.')
  
  if (missing(eoi))
    stop('Please specify the event of interest (eoi).')
  if (!eoi %in% unique(event))
    stop('eoi (event of interest) should be observed in event.')
  
  if (!is.null(tmax) & class(tmax)!='numeric')
    stop('Please specify the time points as a numeric vector.')
  
  library(pseudo)
  ps.ci <- pseudoci(time, event)
  
  evtime <- ps.ci$time
  loc <- which(ps.ci$cause==as.numeric(eoi))
  aalen.johansen.pseudo <- t(ps.ci$pseudo[[loc]])
  
  n <- length(time)
  colnames(aalen.johansen.pseudo) <- paste('n', 1:n, sep = '')
  
  if (is.null(tmax)) {
    rownames(aalen.johansen.pseudo) <- paste('t', 1:length(evtime), sep = '')
    return(list(time = evtime, pseudo.cif = aalen.johansen.pseudo))
  } else {
    tmax <- sort(tmax) # Test tmax <- c(0,evtime,4)
    a <- sum(tmax<evtime[1])
    where <- sapply(tmax[(a+1):length(tmax)], function (x) max(which(evtime<=x)))# location of the last TRUE value
    aalen.johansen.pseudo <- rbind(matrix(0, nrow = a, ncol = n), 
                                   aalen.johansen.pseudo[where,])
    rownames(aalen.johansen.pseudo) <- paste('t', 1:length(tmax), sep = '')
    return(list(time = tmax, pseudo.cif = aalen.johansen.pseudo))
  }
}

Pseudo_rmtl <- function(evtime, pseudo.cif, tau){
  if (length(evtime)!=nrow(pseudo.cif))
    stop('The i_th row of pseudo.cif is the pseudo cif at the i_th evtime.')
  if(tau<=evtime[1]) return(0)
  
  loc <- max(which(evtime<tau))
  evtime <- c(evtime[1:loc],tau)
  pseudo.cif <- pseudo.cif[1:loc,]
  
  pseudo.rmtl <- t(diff(evtime) %*% pseudo.cif)
  colnames(pseudo.rmtl) <- 'pseudo.rmtl'
  
  return(pseudo.rmtl)
}

rmtl_unadj <- function(time, event, group, eoi, tau=NULL){
  
  if (missing(time))
    stop('Please specify the failure time as a numeric vector')
  if (any(is.na(time)))
    stop("Missingness is not allowed in the time vector")
  if (any(time < 0)) 
    stop("'Negative value is not allowed in the time vector")
  
  if (missing(event))
    stop('Please specify the observed event types as an integer vector')
  if (any(is.na(event)))
    stop("Missingness is not allowed in the event vector")
  
  if (missing(eoi))
    stop('Please specify the event of interest (eoi) as a level of event.')
  if (!eoi %in% unique(event))
    stop('eoi (event of interest) should be observed in event.')
  
  event <- as.factor(event)
  data <- data.frame(time, group, event)
  if (is.null(tau)) # by default, tau is the minimum of last event time in two groups.
    tau <- min(aggregate(time ~ group, data = data , FUN = max)[,'time'])
  
  loc <- which(time>tau)
  data$time[loc] <- tau
  data$event[loc] <- 0
  
  grp <- sort(unique(group))
  
  rmtl <- sapply(grp, function(x){
    where <- which(group==x)
    cif <- AJ_cif(data$time[where], data$event[where], eoi)
    sum(diff(c(cif$evtime,tau))*cif$cif)
  })
  
  data.frame(grp=grp, rmtl=rmtl,tau=tau)
}

rmtl_ipw <- function(pseudo.rmtl, group, propensity, dmatrix_ps){
  if (any(is.na(propensity)))
    stop('Missingness is not allowed in propensity score')
  if (any(is.na(group)))
    stop('Missingness is not allowed in group')
  if (sum(!unique(group) %in% c(0,1))>0) 
    stop('Group = 1 for exposure and = 0 for non-exposure')
  if (class(group)!='factor')
    stop('Group is a factor with level 1 and 0')
  
  dmatrix_ps <-data.matrix(dmatrix_ps)
  data <- data.frame(pseudo.rmtl, group, propensity)
  n <- nrow(data)
  p <- ncol(dmatrix_ps)
  
  rmtl1 <- sum(data$pseudo.rmtl[data$group =='1']/data$propensity[data$group=='1'])/n 
  rmtl0 <- sum(data$pseudo.rmtl[data$group=='0']/(1-data$propensity[data$group=='0']))/n
  diff.rmtl.ipw <- rmtl1 - rmtl0
  
  group<- as.numeric(group)-1
  wwt <- apply(dmatrix_ps, 
               MARGIN = 1, 
               FUN = function(x) matrix(x,ncol=1)%*%matrix(x,nrow=1))
  ee <- diag(propensity*(1-propensity))
  Ebb <- matrix(apply(wwt%*%ee, 1, sum),p,p)/n
  ze <- diag(group*pseudo.rmtl*(1-propensity)/propensity+(1-group)*pseudo.rmtl*propensity/(1-propensity)) 
  Hbt <- matrix(apply(t(dmatrix_ps)%*%ze, 1, sum),nrow=1)/n
  M <- apply(dmatrix_ps, 1, function(x) Hbt%*%solve(Ebb)%*%matrix(x,ncol = 1))
  I <- group*pseudo.rmtl/propensity - 
    (1-group)*pseudo.rmtl/(1-propensity) - diff.rmtl.ipw - (group-propensity)*M
  
  var.diff.rmtl.ipw <- sum(I^2)/(n^2)
  
  return(data.frame(group=c(1,0,'1-0','Var.diff.ipw'), 
                    rmtl=c(rmtl1,rmtl0,diff.rmtl.ipw,var.diff.rmtl.ipw)))
}

rmtl_mb <- function(beta, X, tmax, tau){
  # beta is the estimated coefficient from the geese() output
  # X = Covariates + Group assignment(last column)
  
  if (max(tmax)>tau) stop('No need to compute pseudo-observations after tau when calculating RMTL.')
  
  n.tmax <- length(tmax)
  n.beta <- length(beta)
  n.covariate <- n.beta - n.tmax -1
  
  n <- nrow(X)
  
  temp <- as.matrix(X[,1:n.covariate])%*%matrix(beta[(n.tmax+1):(n.beta-1)], nrow = n.covariate)
  
  cloglog.po.cif.D1 <- apply(matrix(rep(beta[1:n.tmax],time=n),byrow = TRUE, nrow = n), 2, FUN = function(x) x+beta[n.beta]+temp) 
  cloglog.po.cif.D0 <- apply(matrix(rep(beta[1:n.tmax],time=n),byrow = TRUE, nrow = n), 2, FUN = function(x) x+temp)
  
  fg.po.cif.D1 <- 1-exp(-exp(cloglog.po.cif.D1)) 
  fg.po.cif.D0 <- 1-exp(-exp(cloglog.po.cif.D0)) 
  
  fg.po.rmtl.D1 <- fg.po.cif.D1 %*% diff(c(tmax,tau))      
  fg.po.rmtl.D0 <- fg.po.cif.D0 %*% diff(c(tmax,tau)) 
  
  return(data.frame(fg.po.rmtl.D0, fg.po.rmtl.D1))
}


rmtl_db <- function(group, ps.rmtl, mb.rmtl.d1, mb.rmtl.d0, propensity){
  
  dat <- data.frame(group, ps.rmtl, mb.rmtl.d1, mb.rmtl.d0, propensity)
  dat.D1 <- dat[dat$group=='1',]
  dat.D0 <- dat[dat$group=='0',]
  n <- nrow(dat)
  
  db.rmtl1.D1 <- (sum(dat.D1$ps.rmtl/dat.D1$propensity 
                      - dat.D1$mb.rmtl.d1*(1-dat.D1$propensity)/dat.D1$propensity)+
                    sum(dat.D0$mb.rmtl.d1))/n
  db.rmtl1.D0 <- (sum(dat.D1$mb.rmtl.d0)+
                    sum(dat.D0$ps.rmtl/(1-dat.D0$propensity)-dat.D0$mb.rmtl.d0*dat.D0$propensity/(1-dat.D0$propensity)))/n
  diff.rmtl.db <- db.rmtl1.D1-db.rmtl1.D0
  
  var.db.diff <- (1/n/n)*(
    sum(((dat.D1$ps.rmtl-(1-dat.D1$propensity)*dat.D1$mb.rmtl.d1)/dat.D1$propensity-dat.D1$mb.rmtl.d0-diff.rmtl.db)**2)+
      sum((dat.D0$mb.rmtl.d1-(dat.D0$ps.rmtl-dat.D0$propensity*dat.D0$mb.rmtl.d0)/(1-dat.D0$propensity)-diff.rmtl.db)**2)
  )
  
  return(c(db.rmtl1.D0 = db.rmtl1.D0, 
           db.rmtl1.D1 = db.rmtl1.D1, 
           diff.rmtl.db = diff.rmtl.db,
           Var = var.db.diff, 
           LowerCI = diff.rmtl.db-1.96*sqrt(var.db.diff),
           UpperCI = diff.rmtl.db+1.96*sqrt(var.db.diff)))
}

rmtl_est <- function(time, group, event, covariate,
                     var_ps, var_fg,
                     tau=1){
  # covariates is the design matrix in the fine and gray model
  
  formula_ps <- as.formula(paste('D ~ ',var_ps, '-1', sep = ''))
  formula_fg <- as.formula(paste('po.cif ~ time + ', var_fg, ' + D -1', sep = ''))
  
  n <- length(time)
  data <- data.frame(Obs.Time = time,
                     D = group,
                     Event = event,
                     covariate)
  data$D <- as.factor(data$D)
  data <- data[order(data$D,decreasing = TRUE),]
  data$ID <- 1:n
  
  data.D1 <- data[data$D=='1',]
  data.D0 <- data[data$D=='0',]
  
  # Unadjusted Estimator of RMTL
  diff.rmtl.unadj <- diff(rmtl_unadj(data$Obs.Time, data$Event, group = data$D, eoi = '1', tau = tau)[,'rmtl'])
  
  # Pseudo-Observation of Aalen-Johansen Estimator
  tmax <- sort(unique(c(data$Obs.Time[data$Event!=0])))
  Pseudo.AJ.D1 <- Pseudo_AJ_cif_2(data.D1$Obs.Time, data.D1$Event, eoi = '1', tmax)
  Pseudo.AJ.D0 <- Pseudo_AJ_cif_2(data.D0$Obs.Time, data.D0$Event, eoi = '1', tmax)
  
  # Pseudo-Observation of RMTL
  data.D1$ps.rmtl <- Pseudo_rmtl(Pseudo.AJ.D1$time,Pseudo.AJ.D1$pseudo.cif, tau=1)
  data.D0$ps.rmtl <- Pseudo_rmtl(Pseudo.AJ.D0$time,Pseudo.AJ.D0$pseudo.cif, tau=1)
  data$ps.rmtl <- c(data.D1$ps.rmtl,data.D0$ps.rmtl)
  
  # IPW Estimator of RMTL
  PS <- glm(formula_ps, data=data,family = 'binomial') # summary(PS.f1) # summary(PS.f1$fitted.values)
  data$propensity <- PS$fitted.values
  diff.rmtl.ipw <- diff(rmtl_ipw(pseudo.rmtl = data$ps.rmtl, group = data$D , propensity = data$propensity)[,'rmtl'])
  
  # Model-Based Estimator of RMTL
  b<-NULL
  for(j in 1:length(tmax)) b <- rbind(b,cbind(data[data$D=='1',],
                                              po.cif = Pseudo.AJ.D1$pseudo.cif[j,],
                                              time = tmax[j]))
  for(j in 1:length(tmax)) b <- rbind(b,cbind(data[data$D=='0',],
                                              po.cif = Pseudo.AJ.D0$pseudo.cif[j,],
                                              time = tmax[j]))
  
  b<-b[order(b$ID),]
  b$time <- as.factor(b$time)
  fit <- geese(formula_fg,
               data=b,
               id=ID, # A vector which identifies the clusters. Data are assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
               jack=T, # if approximate jackknife variance estimate should be computed.
               scale.fix = TRUE, # if the scale should be fixed.
               family = gaussian,
               mean.link = "cloglog",
               corstr = 'independence') 
  
  vars_fg <- all.vars(formula_fg)
  vars_fg <- vars_fg[c(-1,-2,-length(vars_fg))]
  
  data <- cbind(data, rmtl_mb(beta = fit$beta, X = data[,vars_fg], tmax, tau))
  diff.rmtl.mb <- as.numeric(diff(apply(rmtl_mb(beta = fit$beta, X = data[,vars_fg], tmax, tau), MARGIN = 2, mean)))
  
  # Double Robust Estimator
  data.D1 <- data[data$D=='1',]
  data.D0 <- data[data$D=='0',]
  db.rmtl.D1 <- (sum(data.D1$ps.rmtl/data.D1$propensity - data.D1$fg.po.rmtl.D1*(1-data.D1$propensity)/data.D1$propensity)+
                   sum(data.D0$fg.po.rmtl.D1))/n
  db.rmtl.D0 <- (sum(data.D1$fg.po.rmtl.D0)+
                   sum(data.D0$ps.rmtl/(1-data.D0$propensity)-data.D0$fg.po.rmtl.D0*data.D0$propensity/(1-data.D0$propensity)))/n
  diff.rmtl.db <- db.rmtl.D1-db.rmtl.D0
  
  var.db.diff <- (1/n/n)*(
    sum(((data.D1$ps.rmtl-(1-data.D1$propensity)*data.D1$fg.po.rmtl.D1)/data.D1$propensity-data.D1$fg.po.rmtl.D0-diff.rmtl.db)**2)+
      sum((data.D0$fg.po.rmtl.D1-(data.D0$ps.rmtl-data.D0$propensity*data.D0$fg.po.rmtl.D0)/(1-data.D0$propensity)-diff.rmtl.db)**2)
  )
  
  return(data.frame(diff.rmtl.unadj, diff.rmtl.ipw, diff.rmtl.mb, diff.rmtl.db, var.db.diff))
  
}

rmtl_marginal <- function(mu, sigma, beta1, beta2, tau = 1, n = 100000, seed, gamma, rho){
  # Don't consider censoring when calculating the marginal difference in RMTL
  library(etm)
  set.seed(seed)
  Z <-  mvrnorm(n, mu = mu, Sigma = sigma, empirical = FALSE)
  Z.D1 <- cbind(Z,D=rep(1,n))
  Z.D0 <- cbind(Z,D=rep(0,n))
  
  u <- runif(n, 0, 1)
  
  temp <- exp(Z.D1 %*% beta1)*gamma/rho
  Epsilon.D1 <- 1 + rbinom(n, size = 1, prob = exp(temp))
  Time.D1 <- (2-Epsilon.D1)*(1/rho)*cloglog(x=u*(1-exp(temp)), a=1/temp,b=1)+  
    (Epsilon.D1-1)*(-exp(-Z.D1 %*% beta2))*log(1-u)
  where <- which(Time.D1>tau)
  Time.D1[where] <- tau
  Epsilon.D1[where] <- 0 
  
  temp <- exp(Z.D0 %*% beta1)*gamma/rho
  Epsilon.D0 <- 1 + rbinom(n, size = 1, prob = exp(temp))
  Time.D0 <- (2-Epsilon.D0)*(1/rho)*cloglog(x=u*(1-exp(temp)), a=1/temp,b=1)+   
    (Epsilon.D0-1)*(-exp(-Z.D0 %*% beta2))*log(1-u)
  where <- which(Time.D0>tau)
  Time.D0[where] <- tau
  Epsilon.D0[where] <- 0
  
  
  ## Aalen-Johansen Estimator & Unadjusted RMTL estimator
  tra <- matrix(ncol=3,nrow=3,FALSE)
  tra[3,c(1,2)] <- TRUE
  
  # cif.D1 <- AJ_cif(Time.D1, Epsilon.D1, eoi='1')
  dat <- data.frame(id=1:n, from='Alive',to=Epsilon.D1, time=Time.D1)
  dat$to <- as.factor(dat$to) 
  levels(dat$to) <- c('Cen','1','2', 'Alive')
  tr <- etm(dat, state.names=c('1','2', 'Alive'), tra, 'Cen', 0)
  evtime.D1 <- tr$time
  cif1.D1 <- tr$est[3,1,]
  rmtl.D1 <- sum(diff(c(evtime.D1,tau))*cif1.D1)
  
  # cif.D0 <- AJ_cif(Time.D0, Epsilon.D0, eoi='1')
  dat <- data.frame(id=1:n, from='Alive',to=Epsilon.D0, time=Time.D0)
  dat$to <- as.factor(dat$to) 
  levels(dat$to) <- c('Cen','1','2', 'Alive')
  tr <- etm(dat, state.names=c('1','2', 'Alive'), tra, 'Cen', 0)
  evtime.D0 <- tr$time
  cif1.D0 <- tr$est[3,1,]
  rmtl.D0 <- sum(diff(c(evtime.D0,tau))*cif1.D0)
  
  return(data.frame(group=c('0','1', '1-0'), 
                    rmtl.marginal = c(rmtl.D0,rmtl.D1,rmtl.D1-rmtl.D0)))
}

diff.db.ipw  <- function(D, propensity, fg.rmtl.D1, fg.rmtl.D0){
  n <- length(D)
  -sum((D-propensity)*(fg.rmtl.D1/propensity + fg.rmtl.D0/(1-propensity)))/n
}

diff.db.mb   <- function(ps.rmtl.D1, ps.rmtl.D0, 
                         fg.rmtl.D1, fg.rmtl.D0, 
                         propensity.D1, propensity.D0){
  n <- length(ps.rmtl.D1)+length(ps.rmtl.D0)
  (sum((ps.rmtl.D1-fg.rmtl.D1)/propensity.D1)  -sum((ps.rmtl.D0-fg.rmtl.D0)/(1-propensity.D0)))/n
}

diff.db.both <-  function(ps.rmtl.D1, ps.rmtl.D0, 
                          fg.rmtl.D1, fg.rmtl.D0, 
                          propensity.D1, propensity.D0){
  n <- length(ps.rmtl.D1)+length(ps.rmtl.D0)
  (sum(fg.rmtl.D1/propensity.D1) - sum(fg.rmtl.D0/(1-propensity.D0)))/n
}

rmtl <- function(entry=NULL, times, event, eoi=1, group=NULL, weight=NULL, tau=NULL, alpha=.05, yaxismax=1){  
  
  if(sum(times<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=1)>0){print("Error: weights must be greater than 1.")
    }else{
      
      #--- Prep input data ---
      if(is.null(entry)){entry <- rep(0, length(times))}
      if(is.null(group)){group <- as.factor(rep(1, length(times)))}
      if(is.null(weight)){weight <- rep(1, length(times))}
      alldat <- data.frame(entry, times, event, group, weight)
      alldat <- alldat[!is.na(alldat$group) & !is.na(alldat$times),]
      alldat <- alldat[order(group),] 
      
      
      #--- If tau not specified, use minimum time from all groups. Check if provided tau is appropriate. ---
      gg <- length(levels(alldat$group))
      
      if(is.null(tau)){
        
        taui <- rep(NA, gg)
        
        for (i in (1:gg)){
          groupval <- (levels(alldat$group)[i])
          dat_group <- alldat[which(alldat$group==(groupval)),]
          taui[i] <- max(dat_group$times)
        }
        tau <- min(taui)
        
      } else {
        
        tau.error <- rep(0, gg)
        
        for (i in (1:gg)){
          groupval <- (levels(alldat$group)[i])
          dat_group <- alldat[which(alldat$group==(groupval)),]
          tau.error[i] <- ifelse(max(dat_group$times)<tau, 1, 0)
        }
      }
      
      if(sum(tau.error)>0){
        print("Error: observed times do not reach tau in each exposure group. Choose a different tau or leave unspecified for default value of tau.")
      }else{
        
        
        #--- Proceed with RMTL ----
        
        alldat$event[alldat$times>tau] <- 0
        alldat$times[alldat$times>tau] <- tau
        
        rmtl <- rep(NA, length(1:gg))
        groupval <- rep(NA, length(1:gg))
        rmtl.se <- rep(NA, length(1:gg))
        # plot(NULL, xlim=c(0, tau), ylim=c(0,yaxismax), xlab='Time', ylab='Cumulative incidence')
        
        for (g in 1:gg){
          
          #--- Derive CIF and related quantities (theta, hazards, etc) ---
          
          groupval[g] <- (levels(alldat$group)[g])
          data <- alldat[which(alldat$group==(groupval[g])),]
          
          tj <- data$times[data$event!=0]
          tj <- unique(tj[order(tj)])
          num.tj <- length(tj)
          
          num.atrisk <- sapply(tj, function(x) sum(data$weight[data$entry<x & data$times>=x]))
          num.ev1 <- sapply(tj, function(x) sum(data$weight[data$event==eoi & data$times==x]))
          num.ev2 <- sapply(tj, function(x) sum(data$weight[data$event!=eoi & data$event!=0 & data$times==x]))
          num.ev <- num.ev1 + num.ev2
          
          m <- sapply(tj, function(x){sum((data$weight[data$entry<x & data$times>=x])^2)})
          mg <- ((num.atrisk^2)/m)
          
          h1 <- num.ev1/num.atrisk
          h <- num.ev/num.atrisk
          
          s <- cumprod(c(1, 1 - h))
          s <- s[1:length(s)-1]
          
          theta <- s * h1
          cif1 <- cumsum(theta)
          # lines(c(tj, tau), c(cif1, cif1[num.tj]), type="s", col=(g+2), lwd=2)
          
          
          #---  Variance of each theta --- 
          
          a <- c(0,cumsum(num.ev/(mg * (num.atrisk - num.ev))))
          a <- a[1:num.tj]
          
          var.theta <- ((theta)^2) * (((num.atrisk - num.ev1)/(mg * num.ev1)) + a)
          var.theta[is.nan(var.theta)] <- 0
          
          #sum.var.theta <- cumsum(var.theta) 
          
          
          #---  Covariance of thetas --- 
          
          cov.theta <- matrix(NA, nrow=num.tj, ncol=num.tj)
          b <- c(0,cumsum(num.ev/(mg * (num.atrisk - num.ev))))
          
          for(j in 1:(num.tj-1)){
            for(k in (j+1):num.tj){
              #cov.theta[k,j] <- cov.theta[j,k] <- (theta[j]) * (theta[k]) * (-1/num.atrisk[j] + b[j])
              cov.theta[k,j] <- cov.theta[j,k] <- (theta[j]) * (theta[k]) * (-1/mg[j] + b[j])
            }
          }
          
          # Diagonal is variance of thetas
          diag(cov.theta) <- var.theta
          
          
          #---  Covariances of CIF --- 
          
          cov.f10 <- apply(cov.theta, 2, function(x){x[is.na(x)] <- 0;  cumsum(x)})
          cov.f1 <- apply(cov.f10, 1, function(x){x[is.na(x)] <- 0;  cumsum(x)})
          
          var.f1 <- diag(cov.f1) # not sure if this is needed, but for sanity check
          
          #etm <- etmCIF(Surv(times, event != 0) ~ 1, data=data, etype=event, failcode=1)
          #etm.se.0 <- sqrt(unname(etm[[1]]$cov[4,4,]))
          #all.equal(etm.se.0, sqrt(var.f1))
          
          
          #--- RMTL and variance ---
          
          areas <- c(tj[2:num.tj], tau)-tj
          rmtl[g] <- sum(areas*cif1)
          
          cov.weights <- outer(areas,areas)
          cov.f1.weight <- cov.weights * cov.f1
          
          rmtl.var <- sum(cov.f1.weight)
          rmtl.se[g] <- sqrt(rmtl.var) 
          
        }
        
        
        
        #--- Add legend and tau to plot ---
        # abline(v=tau, col=1, lty=3, lwd=2)
        # if(gg>1){
        #   legend('topleft', groupval, lty=rep(1, gg), lwd=rep(2, gg), col=3:(gg+2), cex=.75, bty ="n")
        # }
        
        
        #--- Compare RMTL between groups and compile output---
        
        z <- qnorm(1-alpha/2)
        res <- data.frame(groupval, rmtl, rmtl.se, cil=rmtl-(z*rmtl.se), ciu=rmtl+(z*rmtl.se))
        
        pwc <- ((gg^2)-gg)/2   #number of pairwise comparisons
        
        
        #--- Calculate difference in RMTL if group is not specified
        
        if(pwc>0){
          label.diff <- rep(NA,pwc)
          rmtl.diff <- rep(NA,pwc)
          rmtl.diff.se <- rep(NA,pwc)
          rmtl.diff.cil <- rep(NA,pwc)
          rmtl.diff.ciu <- rep(NA,pwc)
          rmtl.diff.p <- rep(NA,pwc)
          
          res.diff <- data.frame(label.diff, rmtl.diff, rmtl.diff.se, rmtl.diff.cil, rmtl.diff.ciu, rmtl.diff.p)
          l <- 1
          
          for (i in 1:(gg-1)){
            for (ii in (i+1):gg){
              
              #--- RMTL Difference ---
              res.diff[l,]$label.diff <- paste(res[ii,]$groupval, '-', res[i,]$groupval)
              res.diff[l,]$rmtl.diff <- (res[ii,]$rmtl - res[i,]$rmtl)
              res.diff[l,]$rmtl.diff.se <- sqrt(res[ii,]$rmtl.se^2 + res[i,]$rmtl.se^2)
              
              res.diff[l,]$rmtl.diff.cil <- res.diff[l,]$rmtl.diff - z*res.diff[l,]$rmtl.diff.se
              res.diff[l,]$rmtl.diff.ciu <- res.diff[l,]$rmtl.diff + z*res.diff[l,]$rmtl.diff.se
              res.diff[l,]$rmtl.diff.p <- 2*(1-pnorm(abs(res.diff[l,]$rmtl.diff)/res.diff[l,]$rmtl.diff.se))
              
              l <- l+1
            }
          }
          
        }
        
        
        if(pwc>0){
          
          allres <- list(rmtl=res, rmtl.diff=res.diff)
          
          # cat("RMTL per group, tau=", tau, "\n\n", sep="")
          # rmtl.round <- round(allres$rmtl[,c(2:5)],3)
          # colnames(rmtl.round) <- c("RMTL", "SE", "CIL", "CIU")
          # rownames(rmtl.round) <- c(levels(res[,1]))
          # print(rmtl.round)
          # cat("\n\n")
          # 
          # cat ("RMTL Differences, tau=", tau, "\n\n", sep="")
          # rmtl.diff.round <- round(allres$rmtl.diff[c(2:6)],3)
          # colnames(rmtl.diff.round) <- c("RMTL Diff.", "SE", "CIL", "CIU", "p")
          # rownames(rmtl.diff.round) <- c(allres$rmtl.diff[,1])
          # print(rmtl.diff.round)
          # cat("\n\n")
          
          return(allres)
          
        } else { # No groups, thus no pairwise comparison
          
          # cat("RMTL per group, tau=", tau, "\n\n", sep="")
          # colnames(res) <- c("Group", "RMTL", "SE", "CIL", "CIU")
          # rownames(res) <- c(paste("Group", res$Group,' '))
          # print(round(res[c(2:5)],3))
          # cat("\n\n")
          
          return(res)
          
        }
        
      }  
    }
  }
}



