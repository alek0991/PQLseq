pqlseq <- function(RawCountDataSet, Phenotypes, Covariates=NULL, RelatednessMatrix=NULL, LibSize=NULL, 
                   fit.model="BMM", fit.method = "AI.REML", fit.maxiter=500, fit.tol=1e-5, verbose=FALSE, ...) {
  start_time <- Sys.time()
  CountData   <- RawCountDataSet
  numVar <- dim(CountData)[1]
  numIDV <- dim(CountData)[2]
  
  if(is.null(Covariates)){
    numCov <- 0
  }else{
    Covariates <- as.matrix(Covariates)
    numCov     <- dim(Covariates)[2]
  }
  if (verbose) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  if(verbose) cat(paste("number of individuals: ", numIDV,"\n"))
  if(verbose) cat(paste("number of variants: ", numVar,"\n"))
  if(verbose) cat(paste("number of covariates: ", numCov,"\n"))
  if(verbose) cat(paste("tol: ", fit.tol,"\n"))
  if(verbose) cat(paste("maxiter: ", fit.maxiter,"\n"))
  
  
  
  CountData  <- as.matrix(CountData)
  Phenotypes <- as.matrix(Phenotypes)
  RelatednessMatrix <- as.matrix(RelatednessMatrix)
  RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix, doSym=T)$mat)
  RelatednessMatrix <- list(RelatednessMatrix, diag(numIDV))
  
  # ******* Binomial Mixed Model *******
  if(fit.model == "BMM"){ 
    LibSize <- as.matrix(LibSize)
    
    ratio               <- CountData/LibSize
    ratio[is.na(ratio)] <- 0
    flag                <- ratio>1.0
    sumflag             <- apply(flag,1, sum)
    idx                 <- which(sumflag>0)
    
    if (length(idx)>0){
      CountData <- CountData[-idx,]
      LibSize   <- LibSize[-idx,]
    }else{
      CountData <- CountData
      LibSize   <- LibSize
    }
    
    numVar <- dim(CountData)[1]
    numIDV <- dim(CountData)[2]	
    iVar   <- 1
    
    
    LibSize <- as.matrix(LibSize)
    
    if(numCov == 0){
      model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,])
      idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.omit)),
                      rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.pass)))
    }else{
      model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,] )
      idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.omit)),
                      rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.pass)))
    }
    
    model0$numTotal <- LibSize[iVar,idx]
    model0$numSucc  <- CountData[iVar,idx]
    
    tmpRelatednessMatrix <- RelatednessMatrix
    for(ik in seq_len(length(tmpRelatednessMatrix)) ) {
      tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]
    }
    names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")
    model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix, verbose = verbose, maxiter = fit.maxiter, tol = fit.tol))
    model1$elapsedtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    
    params_to_display <- c("elapsedtime", "iter", 'maxdiff', "converged", "intercept", 'se_intercept', "beta", 'se_beta', "pvalue", 'tau1', 'tau2', 'h2', 'sigma2')
    max_param_name_length <- max(nchar(params_to_display))
    if(verbose) cat("------------ Summary ------------\n")
    summary <- c()
    for (param_name in params_to_display) {
      summary <- c(summary, model1[param_name])
      spaces <- paste(rep(" ", max_param_name_length - nchar(param_name)), collapse = "")
      if(verbose) cat(paste0(param_name, ':', spaces), model1[[param_name]], "\n")
    }
    model1$summary <- data.frame(summary)
    return(model1)
  }
  
}# end function PQLseq


##########################################################
#           	   PQLseq FIT FUNCTION					           #
##########################################################

PQLseq.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {
  names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
  if(method.optim == "AI") {
    model1 	<- PQLseq.AI(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
  }else{
    model1 <- NULL
  }
  return(model1)
}

##########################################################
#       PQLseq FIT AVERAGE INFORMATION FUNCTION			     #
##########################################################

PQLseq.AI <- function(model0, RelatednessMatrix, tau = c(1,1,0), fixtau = c(1,0,1), maxiter = 500, tol = 1e-5, verbose = FALSE) {
  y <- model0$numSucc
  numIDV <- length(y)
  
  
  
  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  mu.eta <- family$mu.eta(eta)
  
  mu.eta <- model0$numTotal*mu.eta
  D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
  mu <- model0$numTotal*mu
  
  
  Y <- eta + (y - mu)/mu.eta	
  X <- model.matrix(model0)
  alpha <- model0$coef
  
  for (iter in seq_len(maxiter)) {
    tic = Sys.time()
    if(verbose) cat('------------ iter:', iter, '------------')
    alpha0 	<- alpha
    tau0 	<- tau
    model1 	<- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)
    
    tau <- as.numeric(model1$tau)
    alpha <- as.numeric(model1$alpha)
    if(verbose) cat('\nfixtau(t):', fixtau, '\ntau(t):', tau0, '\ntau(t+1):', tau, '\nalpha(t): ', alpha0, '\nalpha(t+1): ', alpha, '\n')
    
    cov <- as.matrix(model1$cov)
    eta <- as.numeric(model1$eta)
    
    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(family$variance(mu))
    
    mu.eta <- model0$numTotal*mu.eta
    D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
    mu <- model0$numTotal*mu
    
    Y <- eta + (y - mu)/mu.eta
    diff_alpha_tau = 2*c(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))
    maxdiff <- max(diff_alpha_tau)
    if(verbose) cat('tol:', tol, '\n')
    if(verbose) cat('maxdiff:', maxdiff, '\n')
    if(verbose) cat('diff_alpha_tau:', diff_alpha_tau, '\n')
    if(verbose) cat('(diff(alpha), diff(tau))<tol:', 1*(diff_alpha_tau<tol), '\n')
    if(verbose) cat('Time:', as.numeric(difftime(Sys.time(), tic, units = "secs")), 'secs\n')
    
    converged <- FALSE
    if(maxdiff<tol) {
      converged <- TRUE
      break
    }
    if(max(abs(alpha), tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {
      break
    }
  }
  
  res <- y - mu

  P <- model1$P
  
  
  model1$intercept   <- model1$alpha[1]
  model1$se_intercept<- sqrt(diag(model1$cov)[1] )
  
  model1$beta        <- model1$alpha[length(model1$alpha)]
  model1$se_beta     <- sqrt(diag(model1$cov)[length(model1$alpha)])
  
  model1$pvalue      <- pchisq( (model1$beta/model1$se_beta)^2, 1, lower.tail = F)
  
  model1$tau1        <- model1$tau[2]
  model1$tau2        <- model1$tau[3]
  
  model1$sigma2      <- model1$tau1+model1$tau2
  model1$h2          <- model1$tau1/(model1$sigma2)
  
  model1$maxdiff <- maxdiff
  model1$iter <- iter
  model1$mu <- mu
  model1$Y <- Y
  model1$y <- y
  model1$res <- res
  model1$converged <- converged
  model1$X <- X
  model1$K <- RelatednessMatrix[[1]]
  model1$V <- model1$sigma2*model1$K
  model1$D <- D
  model1$family = family
  model1$numTotal = model0$numTotal
  return(model1)
}# end function PQLseq.AI 

#########################################
#             CODE END                  #
#########################################
