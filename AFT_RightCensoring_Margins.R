### Background functions: 

# Useful trycatch variation:
myTryCatch <- function(expr){
  
  warn <- err <- NULL
  
  value <- withCallingHandlers(
    
    tryCatch(expr, error=function(e) {
      
      err <<- e
      
      NULL
    }), warning=function(w) {
      
      warn <<- w
      
      invokeRestart("muffleWarning")
    })
  
  list(value=value, warning=warn, error=err)
}

## helper functions to ensure smooth fitting process: 
pdffz <- function(input){
  
  pdiff <- ifelse(input > 0.999999, 0.999999, input) 
  pdiff <- ifelse(pdiff < 1e-16, 1e-16, pdiff) 
  
  return(pdiff)
  
}

dvffz <- function(derivative){
  
  deriv <- ifelse(is.na(derivative), .Machine$double.eps, derivative)
  deriv <- ifelse(deriv == Inf, 8.218407e+20, deriv)
  deriv <- ifelse(deriv == -Inf, -8.218407e+20, deriv)
  
  return(deriv)
}

weighted.sd <- function(x, w, ...) {
  if (missing(w))
    w <- rep(1, length(x))
  m <- weighted.mean(x, w, ...)
  var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
  return(sqrt(var))
}

### other stuff because these did not load for me... 
check_stabilization <- function(stabilization = c("none", "MAD", "L2")) {
  stabilization <- match.arg(stabilization)
  ## check if old stabilization interface is used and issue a warning
  if (getOption("gamboostLSS_stab_ngrad")) {
    warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)"),
            " is deprecated.\n", "Use argument ", sQuote("stabilization"),
            " in the fitting family. See ?Families for details.")
    if (stabilization == "none")
      warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
  }
  stabilization
}

stabilize_ngradient <- function(ngr, w = 1, stabilization) {
  ## set which to MAD if gamboostLSS_stab_ngrad = TRUE and which == "none"
  if (stabilization == "none" && getOption("gamboostLSS_stab_ngrad"))
    stabilization <- "MAD"
  ## stabilization using the mean absolute deviation (MAD)
  if (stabilization == "MAD") {
    div <- weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                           w = w, na.rm = TRUE)
    div <- ifelse(div < 0.0001, 0.0001, div)
    ngr <- ngr / div
  }
  if (stabilization == "L2") {
    div <- sqrt(weighted.mean(ngr^2, w =w,  na.rm = TRUE))
    div <- ifelse(div < 1e-04, 1e-04, div)
    div <- ifelse(div > 1e+04, 1e+04, div)
    ngr <- ngr / div
  }
  ngr
}


########################################################################################### WEIBULL
Custom_WeibullMu <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  # neg. log-likelihood
  loss <- function(sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dermu <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 
    
    derS_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dermu * exp(f) ) + (1-censind) * ( (1/SurvT) * derS_dermu * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( mu )) {
      temp1 <- mu
      temp1 <- log( temp1 )
      
      RET <- temp1
    }else {
      
      
      RET <- log(weighted.mean((y[,1] + weighted.mean(y[,1], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE)) 
      
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "Weibull distribution: mu (log link)")
}

Custom_WeibullSigma <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  # neg. log-likelihood
  loss <- function(mu, y, f=f){
    
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = mu, shape = param)
    
    SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, y = y, f = f))
  }
  
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = mu, shape = param)
    
    SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dersigma <- ( (1/(param)) + log(time) - log(mu) - (time/ (mu))^(param) * log(time / (mu) )  )
    
    derS_dersigma <- - ( ( exp(- (time / (mu) )^(param) ) * ( (time / (mu) )^(param) * log(time / (mu) ) ) ) )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dersigma * exp(f) ) + (1 - censind) * ( ( 1/SurvT ) * derS_dersigma * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  offset <- function(y, w){
    
    if (!is.null( sigma )) {
      
      
      temp2 <- log( sigma )
      
      RET <- temp2
      
    }else{
      
      sigma_temp <- rep(0.1, length(y[,1]))
      
      RET <- log(mean(sigma_temp))
      
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "Weibull distribution: sigma (log link)")
}

#------------------ complete gamboostLSS families
Custom_WeibullFamily <- function (mu = NULL,  sigma = NULL, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization(stabilization)
  
  Families(   mu = Custom_WeibullMu(mu = mu, sigma = sigma,  stabilization = stabilization ), 
              sigma = Custom_WeibullSigma(mu = mu, sigma = sigma, stabilization = stabilization ), 
              name = "Weibull distribution for right-censored data")
}

########################################################################################## LOG LOGISTIC
dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
  
  if (any(mu < 0)) 
    stop(paste("mu must be positive", "\n", ""))
  
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  
  z <- (x/mu)^(sigma)
  
  loglik <- log(sigma) + (sigma) * (log(x) - log(mu)) - log(x) - 2 * log(1 + z)
  
  if (log == FALSE) 
    ft <- exp(loglik)
  else ft <- loglik
  ft
}

pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
  
  if (any(mu < 0)) 
    stop(paste("mu must be positive", "\n", ""))
  
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  
  log_cdf <- log(1/(1 + (x/mu)^(-sigma)))
  
  if (log == FALSE) 
    ft <- exp(log_cdf)
  else ft <- log_cdf
  ft
}

Custom_LogLogisticMu <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  # neg. log-likelihood
  loss <- function(sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dLogLogistic(x = time, mu = param, sigma = sigma)
    
    SurvT <- 1 - pLogLogistic(x = time, mu = param, sigma = sigma)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # 
    param <- exp(f)
    
    pdfT <- dLogLogistic(x = time, mu = param, sigma = sigma)
    
    SurvT <- 1 - pLogLogistic(x = time, mu = param, sigma = sigma)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dermu <- ( -(sigma)/param + 2 * (sigma) * (time / param)^(sigma) / (param * (1 + (time / param)^(sigma))) )
    
    derS_dermu <- - ( - ( (sigma)/param * ((time)/param)^(-(sigma)) / (1 + ((time)/param)^(-(sigma)))^2 ) )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dermu * exp(f) ) + (1-censind) * ( ( 1/SurvT ) * derS_dermu * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( mu )) {
      temp1 <- mu
      temp1 <- log( temp1 )
      
      RET <- temp1
    }else {
      
      
      RET <- log(weighted.mean((y[,1] + weighted.mean(y[,1], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE)) 
      
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "LogLogistic distribution: mu (log link)")
  
}

Custom_LogLogisticSigma <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  # neg. log-likelihood
  loss <- function(mu, y, f=f){
    
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <- dLogLogistic(x = time, mu = mu, sigma = param)
    
    SurvT <- 1 - pLogLogistic(x = time, mu = mu, sigma = param)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, y = y, f = f))
  }
  
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <- dLogLogistic(x = time, mu = mu, sigma = param)
    
    SurvT <- 1 - pLogLogistic(x = time, mu = mu, sigma = param)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dersigma <- ( ( (1/param) + log(time) - log(mu) - 2/(1 + ((time) / (mu))^(param)) * ((time) / (mu))^(param) * log((time)/(mu)) ) )
    
    derS_dersigma <- - ( ((time) / (mu))^(-(param)) * log((time)/(mu)) / (1 + ((time) / (mu))^(-(param)))^2  )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dersigma * exp(f) ) + (1 - censind) * ( ( 1/SurvT ) * derS_dersigma * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  offset <- function(y, w){
    
    if (!is.null( sigma )) {
      
      
      temp2 <- log( sigma )
      
      RET <- temp2
      
    }else{
      
      sigma_temp <- rep(0.1, length(y[,1]))
      
      RET <- log(mean(sigma_temp))
      
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "LogLogistic distribution: sigma (log link)")
}

Custom_LogLogisticFamily <- function (mu = NULL,  sigma = NULL, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization(stabilization)
  
  Families(   mu = Custom_LogLogisticMu(mu = mu, sigma = sigma,  stabilization = stabilization ), 
              sigma = Custom_LogLogisticSigma(mu = mu, sigma = sigma, stabilization = stabilization ), 
              name = "LogLogistic distribution for right-censored data")
}

################################################################################################ EXPONENTIAL
Custom_ExponentialFamily <- function(){
  
  mboost::Family(ngradient = function(y, f, w = 1){
    
    theta <- exp(f)
    
    PDF <- dexp(x = y[,1], rate = theta)
    
    SurvFunc <- pexp(q = y[,1], rate = theta, lower.tail = FALSE)
    
    PDF       <- dvffz(PDF)
    SurvFunc  <- pdffz(SurvFunc)
    
    derPDF_derTheta <- exp( - theta * y[,1] ) * (1 - y[,1] * theta)
    
    derSurv_derTheta <- - y[,1] *  exp( - theta * y[,1] )
    
    ngr <- ( y[,2] * (1 / PDF) * derPDF_derTheta * exp(f) + ( 1 - y[,2] ) * (1 / SurvFunc ) * derSurv_derTheta )
    
    return(ngr)
    
  },
  
  
  loss = function(y, f){ 
    
    theta <- exp(f)
    
    PDF <- dexp(x = y[,1], rate = theta)
    
    SurvFunc <- pexp(q = y[,1], rate = theta, lower.tail = FALSE)
    
    PDF       <- dvffz(PDF)
    SurvFunc  <- pdffz(SurvFunc)
    
    loglik <-  ( y[,2] * log(PDF) + ( 1 - y[,2] ) * log( SurvFunc ) )
    
    negloglik <- -loglik
    
    return(negloglik)
    
  },
  
  offset = function(y, w = 1){
    
    temp <- weighted.mean(y[,1] * y[,1], w = w)
    
    temp <- log(temp)
    
    return(temp)
    
  },
  
  name = "Exponential distribution for right-censored time-to-event data")
}


################################################################################################ LOG NORMAL
Custom_LogNormalMu <- function(mu = NULL, sigma = NULL, stabilization){
  
  # neg. log-likelihood
  loss <- function(sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- f
    
    
    pdfT <-  dnorm(log(time), mean = param, sd = sigma)
    
    SurvT <- 1 - pnorm(log(time), mean = param, sd = sigma)
    
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- f
    
    pdfT <- dnorm(x = log(time), mean = param, sd = sigma)
    
    SurvT <- 1 - pnorm(q = log(time), mean = param, sd = sigma)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derpdf_dermu <-  exp(-( 0.5 * ((log(time) - param)^2/sigma^2))) * (log(time) - param)/(sigma^2 * sqrt(2 * (pi * sigma^2)))   
    
    derS_dermu <-  dnorm(log(time), mean = param, sd = sigma)
    
    
    
    #### negative gradient:
    ngr <- censind * ( (1 / pdfT) *  derpdf_dermu * 1 ) + (1-censind) * ( ( 1/SurvT ) * derS_dermu * 1 )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( mu )) {
      temp1 <- mu
      temp1 <- ( temp1 )
      
      RET <- temp1
    }else {
      
      
      RET <- ( weighted.mean((log(y[,1]) + weighted.mean(log(y[,1]), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) ) 
      
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) f, 
                  offset = offset,
                  name = "Log-normal distribution: mu (identity link)")
  
}

Custom_LogNormalSigma <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  # neg. log-likelihood
  loss <- function(mu, y, f=f){
    
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <-  dnorm(log(time), mean = mu, sd = param)
    
    SurvT <- 1 - pnorm(log(time), mean = mu, sd = param)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, y = y, f = f))
  }
  
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dnorm(x = log(time), mean = mu, sd = param)
    
    SurvT <- 1 - pnorm(q = log(time), mean = mu, sd = param)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    derpdf_dersigma <- -((1 - (log(time) - mu)^2/param^2) * exp(-(0.5 * ((log(time) - mu)^2/param^2)))/(param * sqrt(2 * (pi * param^2))))
    
    derS_dersigma <- - -dnorm((log(time) - mu)/param)*(log(time) - mu)/param^2
    
    
    
    #### negative gradient:
    ngr <- censind * ( ( 1 / pdfT ) * derpdf_dersigma * exp(f) ) + (1 - censind) * ( ( 1/SurvT ) * derS_dersigma * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  offset <- function(y, w){
    
    if (!is.null( sigma )) {
      
      
      temp2 <- log( sigma )
      
      RET <- temp2
      
    }else{
      
      sigma_temp <- rep(weighted.sd( (y[,1]) , w = w, na.rm = TRUE), length(y[,1])) 
      
      RET <- log(mean(sigma_temp))
      
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "Log-normal distribution: sigma (log link)")
}

Custom_LogNormalFamily <- function (mu = NULL,  sigma = NULL, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization(stabilization)
  
  Families(   mu = Custom_LogNormalMu(mu = mu, sigma = sigma,  stabilization = stabilization ), 
              sigma = Custom_LogNormalSigma(mu = mu, sigma = sigma, stabilization = stabilization ), 
              name = "Log-normal distribution for right-censored data")
}
