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




##### Try to do the sophisticated version using code from the copulas!
##### The main parameter is the CURE fraction, and you can select the margin (as well as the link for the cure)
RC_CureModel_Weibull_MU <- function(mu = NULL, sigma = NULL,  nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(sigma, nu, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
    
    SurvT_cure <- (nu) * SurvT + 1 - nu
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    
    lik <- censind * ( log(nu) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, nu = nu, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
    
    SurvT_cure <- (nu) * SurvT + 1 - nu
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    derlpdf_dermu <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 
    
    derS_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )
    
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dermu * exp(f) ) + (1-censind) * ( (1/SurvT_cure) * nu * derS_dermu * exp(f) )
    
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
                  name = "Weibull distribution with cure: mu (log link)")
}
  
  
RC_CureModel_Weibull_SIGMA <- function(mu = NULL, sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
    
    
    # neg. log-likelihood
    loss <- function(mu, nu, y, f=f){
      
      
      time <- y[,1]
      censind <- y[,2]
      
      # in here this is SIGMA or vartheta2
      param <- exp(f)
      
      pdfT <- dweibull(x = time, scale = mu, shape = param)
      
      SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = FALSE)
      
      
      ## Small check
      pdfT <- dvffz(pdfT)
      
      SurvT <- pdffz(SurvT)
      
      SurvT_cure <- (nu) * SurvT + 1 - nu
      
      lik <- censind * ( log(param) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
      
      neglik <- - lik
      
      return(neglik)
      
      
    }
    
    # risk is sum of loss
    risk <- function(y, f, w = 1) {
      sum(w * loss(mu = mu, nu = nu, y = y, f = f))
    }
    
    
    # ngradient 
    ngradient <- function(y, f, w = 1){
      
      time <- y[,1]
      censind <- y[,2]
      
      # in here this is SIGMA or vartheta2
      param <- exp(f)
      
      pdfT <- dweibull(x = time, scale = mu, shape = param)
      
      SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = FALSE)
      
      
      SurvT_cure <- (nu) * SurvT + 1 - nu
      
      ## Small check
      pdfT <- dvffz(pdfT)
      
      SurvT <- pdffz(SurvT)
      
      SurvT_cure <- pdffz(SurvT_cure)
      
      
      derlpdf_dersigma <- ( (1/(param)) + log(time) - log(mu) - (time/ (mu))^(param) * log(time / (mu) )  )
      
      derS_dersigma <- - ( ( exp(- (time / (mu) )^(param) ) * ( (time / (mu) )^(param) * log(time / (mu) ) ) ) )
      
      
      #### negative gradient:
      ngr <- censind * ( derlpdf_dersigma * exp(f) ) + (1 - censind) * ( ( 1/SurvT_cure ) * nu * derS_dersigma * exp(f) )
      
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
                    name = "Weibull distribution with cure: sigma (log link)")
  }
    
    
RC_CureModel_Weibull_NU <- function(mu = NULL, sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(mu, sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
  
    # in here this is NU or vartheta3 or pi in the cure model notation:
    param <- exp(f) / (1 + exp(f))
    
    pdfT <- dweibull(x = time, scale = mu, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = mu, shape = sigma, lower.tail = FALSE)
    
    SurvT_cure <- (param) * SurvT + 1 - param
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    lik <- censind * ( log(param) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, sigma = sigma, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f) / ( 1 + exp(f) )
    
    pdfT <- dweibull(x = time, scale = mu, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = mu, shape = sigma, lower.tail = FALSE)
    
    SurvT_cure <- (param) * SurvT + 1 - param
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- pdffz(SurvT_cure)
    
    lik <- censind * ( log(1 - param) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
    
    
    ## Derivative of nu w.r.t. eta nu: (logistic response function)
    derNu_deta_nu <- exp(f) / ( 1 + exp(f) )^2
    
    
    #### negative gradient:
    ngr <- censind * ( (1/(param)) * derNu_deta_nu ) + (1 - censind) * ( (1 / (SurvT_cure)) * ( SurvT * derNu_deta_nu - derNu_deta_nu ) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( nu )) {
      temp1 <- nu
      temp1 <- log( nu / (1 - nu) )
      
      RET <- temp1
    }else {
      
      ## Just the censoring rate as a crude starting value
      RET_tmp <- ( weighted.mean((y[,2] + weighted.mean(y[,2], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) ) 
      
      RET <- log( RET_tmp / (1 - RET_tmp) )
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f)/( 1 + exp(f) ), 
                  offset = offset,
                  name = "Weibull distribution with cure: nu (logit link)")
}
  

Weibull_Cure <- function (mu = NULL,  sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
    
    stabilization <- check_stabilization(stabilization)
    
    Families(   mu = RC_CureModel_Weibull_MU(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization), 
                sigma = RC_CureModel_Weibull_SIGMA(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization), 
                nu = RC_CureModel_Weibull_NU(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization),
                name = "Weibull distribution for right-censored data with cure fraction (mixture model)")
}



##### Try to do the sophisticated version using code from the copulas!
##### The main parameter is the CURE fraction, and you can select the margin (as well as the link for the cure)
RC_CureModelPromotion_Weibull_MU <- function(mu = NULL, sigma = NULL,  nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(sigma, nu, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = TRUE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    F_t <- pdffz(SurvT)
    
    
    lik <- censind * ( log(nu) + log(pdfT) - nu * F_t ) + (1 - censind) * ( - nu * F_t )
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, nu = nu, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = TRUE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dermu <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 
    
    derS_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )
    
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dermu * exp(f) - nu * (-derS_dermu) * exp(f) ) + (1-censind) * (  - nu * (-derS_dermu) * exp(f) )
    
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
                  name = "Weibull distribution with cure: mu (log link)")
}


RC_CureModelPromotion_Weibull_SIGMA <- function(mu = NULL, sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(mu, nu, y, f=f){
    
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = mu, shape = param)
    
    SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = TRUE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    F_t <- pdffz(SurvT)
    
    
    lik <- censind * ( log(param) + log(pdfT) - nu * F_t) + (1 - censind) * ( -nu * F_t )
    
    neglik <- - lik
    
    return(neglik)
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, nu = nu, y = y, f = f))
  }
  
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = mu, shape = param)
    
    SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = TRUE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    F_t <- pdffz(SurvT)
    
    
    derlpdf_dersigma <- ( (1/(param)) + log(time) - log(mu) - (time/ (mu))^(param) * log(time / (mu) )  )
    
    derS_dersigma <- - ( ( exp(- (time / (mu) )^(param) ) * ( (time / (mu) )^(param) * log(time / (mu) ) ) ) )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dersigma * exp(f) - nu * (-derS_dersigma) * exp(f) ) + (1 - censind) * ( - nu * (-derS_dersigma) * exp(f) ) 
    
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
                  name = "Weibull distribution with cure: sigma (log link)")
}


RC_CureModelPromotion_Weibull_NU <- function(mu = NULL, sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(mu, sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    
    # in here this is NU or vartheta3 in the cure model notation:
    param <- exp(f) 
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = TRUE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    F_t <- pdffz(SurvT)
    
    
    
    lik <- censind * ( log(param) + log(pdfT) - param * F_t ) + (1 - censind) * ( - param * F_t )
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, sigma = sigma, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f) 
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = TRUE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    F_t <- pdffz(SurvT)
    
    
    ## Derivative of nu w.r.t. eta nu: (logistic response function)
    derNu_deta_nu <- exp(f) 
    
    
    #### negative gradient:
    ngr <- censind * ( (1/param) * derNu_deta_nu - derNu_deta_nu * F_t ) + (1 - censind) * ( - derNu_deta_nu * F_t )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( nu )) {
      temp1 <- nu
      temp1 <- log( nu )
      
      RET <- temp1
    }else {
      
      ## Just the censoring rate as a crude starting value
      RET_tmp <- ( weighted.mean((y[,2] + weighted.mean(y[,2], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) ) 
      
      RET <- log( RET_tmp )
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "Weibull distribution with cure: nu (log link)")
}


Weibull_Cure_Promotion <- function (mu = NULL,  sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization(stabilization)
  
  Families(   mu = RC_CureModelPromotion_Weibull_MU(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization), 
              sigma = RC_CureModelPromotion_Weibull_SIGMA(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization), 
              nu = RC_CureModelPromotion_Weibull_NU(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization),
              name = "Weibull distribution for right-censored data with cure fraction (promotion model)")
}




###################################################################################################################
### Log Normal distribution: 

RC_CureModel_LogNormal_MU <- function(mu = NULL, sigma = NULL,  nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(sigma, nu, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    # in here this is MU or vartheta1
    param <- f
    
    
    pdfT <-  dnorm(log(time), mean = param, sd = sigma)
    
    SurvT <- 1 - pnorm(log(time), mean = param, sd = sigma)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- (nu) * SurvT + 1 - nu
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    lik <- censind * ( log(nu) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, nu = nu, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    # in here this is MU or vartheta1
    param <- f
    
    
    pdfT <-  dnorm(log(time), mean = param, sd = sigma)
    
    SurvT <- 1 - pnorm(log(time), mean = param, sd = sigma)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    
    SurvT_cure <- (nu) * SurvT + 1 - nu
    SurvT_cure <- pdffz(SurvT_cure)
    
    derpdf_dermu <-  exp(-( 0.5 * ((log(time) - param)^2/sigma^2))) * (log(time) - param)/(sigma^2 * sqrt(2 * (pi * sigma^2)))   
    
    derS_dermu <-  dnorm(log(time), mean = param, sd = sigma)
    
    der_uncens <- (1/pdfT) * derpdf_dermu
    
    der_cens <- (1/SurvT_cure) * nu * derS_dermu

    
    #### negative gradient:
    ngr <- censind * der_uncens + (1-censind) * der_cens
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( mu )) {
      temp1 <- mu
      temp1 <- ( temp1 )
      
      RET <- temp1
    }else {
      
      
      RET <- (weighted.mean((y[,1] + weighted.mean(y[,1], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE)) 
      
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) f, 
                  offset = offset,
                  name = "Log-Normal distribution with cure: mu (identity link)")
}


RC_CureModel_LogNormal_SIGMA <- function(mu = NULL, sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(mu, nu, y, f=f){
    
    
    time <- y[,1]
    censind <- y[,2]
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    
    pdfT <-  dnorm(log(time), mean = mu, sd = param)
    
    SurvT <- 1 - pnorm(log(time), mean = mu, sd = param)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- (nu) * SurvT + 1 - nu
    
    lik <- censind * ( log(param) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
    
    neglik <- - lik
    
    return(neglik)
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, nu = nu, y = y, f = f))
  }
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    # in here this is SIGMA or vartheta2
    param <- exp(f)
    
    pdfT <-  dnorm(log(time), mean = mu, sd = param)
    
    SurvT <- 1 - pnorm(log(time), mean = mu, sd = param)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- (nu) * SurvT + 1 - nu
    
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    derpdf_dersigma <- -((1 - (log(time) - mu)^2/param^2) * exp(-(0.5 * ((log(time) - mu)^2/param^2)))/(param * sqrt(2 * (pi * param^2))))
    
    derS_dersigma <- - -dnorm((log(time) - mu)/param)*(log(time) - mu)/param^2
    
    der_uncens <- (1/pdfT) * derpdf_dersigma * param
    
    der_cens <- (1/SurvT_cure) * nu * derS_dersigma * param
    
    #### negative gradient:
    ngr <- censind * der_uncens + (1 - censind) * der_cens
    
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
                  name = "Log-Normal distribution with cure: sigma (log link)")
}


RC_CureModel_LogNormal_NU <- function(mu = NULL, sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  
  # neg. log-likelihood
  loss <- function(mu, sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is NU or vartheta3 or pi in the cure model notation:
    param <- exp(f) / (1 + exp(f))
    
    pdfT <-  dnorm(log(time), mean = mu, sd = sigma)
    
    SurvT <- 1 - pnorm(log(time), mean = mu, sd = sigma)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- (param) * SurvT + 1 - param
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    lik <- censind * ( log(param) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu=mu, sigma=sigma, y=y, f=f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    param <- exp(f) / ( 1 + exp(f) )
    
    # 
    pdfT <-  dnorm(log(time), mean = mu, sd = sigma)
    
    SurvT <- 1 - pnorm(log(time), mean = mu, sd = sigma)
    
    ## Small check
    pdfT <- dvffz(pdfT)
    SurvT <- pdffz(SurvT)
    
    SurvT_cure <- (param) * SurvT + 1 - param
    SurvT_cure <- pdffz(SurvT_cure)
    
    
    ## Derivative of nu w.r.t. eta nu: (logistic response function)
    derNu_deta_nu <- exp(f) / ( 1 + exp(f) )^2
    
    der_uncens <-  (1/(param)) * derNu_deta_nu 
    
    der_cens <- (1 / SurvT_cure ) * ( SurvT  - 1 ) * derNu_deta_nu
    
    #### negative gradient:
    ngr <- censind * der_uncens + (1 - censind) * der_cens
    
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( nu )) {
      temp1 <- nu
      temp1 <- log( temp1 / (1 - temp1) )
      
      RET <- temp1
    }else {
      
      ## Just the censoring rate as a crude starting value
      RET_tmp <- ( weighted.mean((y[,2] + weighted.mean(y[,2], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) ) 
      
      RET <- log( RET_tmp / (1 - RET_tmp) )
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f)/( 1 + exp(f) ), 
                  offset = offset,
                  name = "Log-Normal distribution with cure: nu (logit link)")
}


LogNormal_Cure <- function (mu = NULL,  sigma = NULL, nu = NULL, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization(stabilization)
  
  Families(   mu = RC_CureModel_LogNormal_MU(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization), 
              sigma = RC_CureModel_LogNormal_SIGMA(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization), 
              nu = RC_CureModel_LogNormal_NU(mu = mu, sigma = sigma, nu = nu, stabilization = stabilization),
              name = "Log-Normal distribution for right-censored data with cure fraction (mixture model)")
}





##### Only Cure fraction loss: 
CureFraction_Loss <- function(){
  mboost::Family(ngradient = function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    param <- exp(f)/(1 + exp(f))
    
    pdfT <-  dvffz(y[,3])
    
    SurvT <- pdffz(y[,4])
    
    SurvT_cure <- (param) * SurvT + 1 - param
    SurvT_cure <- pdffz(SurvT_cure)
    
  
    ## Derivative of nu w.r.t. eta nu: (logistic response function)
    derNu_deta_nu <- exp(f) / ( 1 + exp(f) )^2
    
    der_uncens <-  (1/(param)) * derNu_deta_nu 
    
    der_cens <- (1 / SurvT_cure ) * ( SurvT  - 1 ) * derNu_deta_nu
    
    #### negative gradient:
    ngr <- censind * der_uncens + (1 - censind) * der_cens
    

    return(ngr)
    
  },
         loss = function(y, f){
           
           time <- y[,1]
           censind <- y[,2]
           
           param <- exp(f)/(1 + exp(f))
           
           pdfT <-  dvffz(y[,3])
           
           SurvT <- pdffz(y[,4])
           
           SurvT_cure <- (param) * SurvT + 1 - param
           SurvT_cure <- pdffz(SurvT_cure)
           
           
           lik <- censind * ( log(param) + log(pdfT) ) + (1 - censind) * log(SurvT_cure)
           
           neglik <- - lik
           
           return(neglik)
         },
  offset = function(y, w = 1){
    
    temp <- weighted.mean(y[,2] * y[,2], w = w)
    
    temp <- log(temp / (1-temp))
    
    return(temp)
    
  },
         name = "Cure fraction with arbitrary margins")
}



























