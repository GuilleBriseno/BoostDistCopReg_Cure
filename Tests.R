library(gamboostLSS)
library(survival)
library(mvtnorm)


## Simulate: 
n.train <- n.mstop <- n.test <- 1000
seed <- 1
p <- 10
corr_toe <- 0.5

n <- n.train + n.mstop

weights.mstop <- c(rep(1, times = n.train), rep(0, times = n.mstop)) 

set.seed(seed)

# generate survival times: 
# #### sample design matrix for train:
BigXMatrix_Train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr_toe^x)))
BigXMatrix_Train <- apply(BigXMatrix_Train, 2, pnorm)

x.train <- BigXMatrix_Train

#### sample design matrix for test / evaluation:
BigXMatrix_Test <- rmvnorm(n = n.test, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr_toe^x)))
BigXMatrix_Test <- apply(BigXMatrix_Test, 2, pnorm) 

x.test <- BigXMatrix_Test

colnames(x.train) <- paste0("X", 1:p)
colnames(x.test) <- paste0("X", 1:p)

beta11     <- c( 2, 0, 0, 0, 0, 0, rep(0, p-6)) 

#beta12    <- c( 0, +1, 0, 1.5, 0, 0, rep(0, p-6)) 
beta12    <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 

betarho   <- c( 0, -2, 0, +2, +0, 0, rep(0, p-6)) 


########################### 30% censoring
eta_11 <- x.train %*% beta11
eta_12 <- -2 + x.train %*% beta12

eta_rho <- x.train %*% betarho

### for the cutting / right-censoring of uniform times:
timecutM1 <- 9.5

mu1     <- (eta_11)
sigma1  <- exp(eta_12)

rho     <- exp(eta_rho) / (1 + exp(eta_rho))

plot(rho)



#sigma1 <- rep(exp(0), length(mu1))


# transform into non-uniform random variables:
#t_true <- rweibull(n = length(mu1), scale = mu1, shape = sigma1)
t_true <- rlnorm(n = length(mu1), meanlog = mu1, sdlog = sigma1)


plot(t_true)


### Binary variable (observation is cured, yes / no)
v_bin <- rbinom(n = length(mu1), size = 1, prob = rho)

table(v_bin)

t_true[which(v_bin == 0)] <- Inf

# Sample censoring times: 
censoring_time <- runif(length(mu1), min = 0, max = timecutM1)

status_m1 <- 1 * (t_true <= censoring_time)

# # censoring rate margin 1: 
table(status_m1)/length(mu1)

# Replace TRUE event times with observed times: 
y <- pmin(t_true, censoring_time)

plot(y)

dat <- data.frame(y, status_m1, x.train)

#### 
# library(mixcure)
# 
# cure_model <- mixcure(Surv(y, status_m1) ~ X1, ~ X2 + X4, 
#                       data = dat, 
#                       lmodel = list(fun = "survreg", 
#                                     dist = "weibull"), 
#                       savedata = F)





################# Check quantities: pdf, Surv, SURV_CURE, pi
# deriv_pdf, deriv_Surv, deriv_pi, 

## parameter true values
eta_11
eta_12

### Cure fraction true values: 
eta_rho

#nu <- exp(eta_rho) / (1 + exp(eta_rho))
#mu1
#sigma1



#### 
time <- dat$y
censind <- dat$status_m1

SurvT <- pweibull(q = time, scale = mu1, shape = sigma1, lower.tail = FALSE)

pdfT <- dweibull(x = time, scale = mu1, shape = sigma1)

plot(SurvT)

## Small check
pdfT <- dvffz(pdfT)
SurvT <- pdffz(SurvT)


SurvT_Cure <- nu * SurvT + (1 - nu)


par(mfrow = c(1, 2))
plot(SurvT) 
plot(SurvT_Cure)
par(mfrow = c(1, 1))


derMu_dereta <- exp(eta_11)
derSigma_dereta <- exp(eta_12)
derNu_dereta <- exp(eta_rho) / (1+exp(eta_rho))^2


### derivatives wrt mu
param <- mu1
sigma <- sigma1
derlpdf_dermu <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 

derS_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )

par(mfrow = c(1, 2))
plot(derlpdf_dermu * derMu_dereta) 
plot(derS_dermu * derMu_dereta )
par(mfrow = c(1, 1))


### Derivatives wrt sigma
param <- sigma1
mu <- mu1
derlpdf_dersigma <- ( (1/(param)) + log(time) - log(mu) - (time/ (mu))^(param) * log(time / (mu) )  )

derS_dersigma <- - ( ( exp(- (time / (mu) )^(param) ) * ( (time / (mu) )^(param) * log(time / (mu) ) ) ) )

par(mfrow = c(1, 2))
plot(derlpdf_dersigma * derSigma_dereta , ylim = c(-2, 2)) 
plot(derS_dersigma * derSigma_dereta )
par(mfrow = c(1, 1))



# Derivatives wrt NU
derNu_dereta
plot(derNu_dereta)

delta <- dat$status_m1


pdfT <- dweibull(x = time, scale = param, shape = sigma, log = FALSE)

SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = TRUE)

SurvT <- pdffz(SurvT)

par(mfrow = c(1,2))
Surv_Cure <- nu * SurvT + (1 - nu)
plot(Surv_Cure)

Surv_Cure <- pdffz(Surv_Cure)
plot(Surv_Cure)
par(mfrow = c(1,1))


# derlogpdf_dermu
derlpdf_derparam <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 

# derSurvivalfunc_dermu
derSt_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )

deriv_uncens <- derlpdf_derparam * param

deriv_cens <- (1/Surv_Cure) * nu * derSt_dermu * param

derivs <- delta * deriv_uncens + (1 - delta) * deriv_cens

plot(derivs)




##############################
param <- exp(eta_12)

# derlogpdf_dersigma
derlpdf_dersigma <- ( (1/(param)) + log(time) - log(mu) - (time/ (mu))^(param) * log(time / (mu) )  )

# derSurvivalfunc_dersigma
derSt_dersigma <- - ( ( exp(- (time / (mu) )^(param) ) * ( (time / (mu) )^(param) * log(time / (mu) ) ) ) )

deriv_uncens_sigma <- derlpdf_dersigma * sigma

deriv_cens_sigma <- (1/Surv_Cure) * nu * derSt_dersigma * sigma

derivs_sigma <- delta * deriv_uncens_sigma + (1 - delta) * deriv_cens_sigma

plot(derivs_sigma)
plot(derivs)











### check with numerical derivatives: x is ETA_PARAM
derlpdf_dermu_func <- function(x, sigma, time){
  
  param <- exp(x)
  
  dermu_deret <- exp(x)
  
  derlpdf_derparam <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 
  
  
  derivs <- derlpdf_derparam * dermu_deret
  
  return(as.numeric(derivs))
  
}

logpdf_mu_func <- function(x, sigma, time){
  
  param <- exp(x)
  
  ret <- dweibull(x = time, scale = param, shape = sigma, log = TRUE)
  
  return(as.numeric(ret))
  
}

derS_dermu_func <- function(x, sigma){}

logpdf_func <- function(x, sigma){}

derlpdf_dersigma_func <- function(x, mu){}

logpdf_func <- function(x, sigma){}

derS_dersigma_func <- function(x, mu){}

logpdf_func <- function(x, sigma){}






derLoss_dermu_func <- function(x, sigma, nu, time, delta){
  
  param <- exp(x)
  
  pdfT <- dweibull(x = time, scale = param, shape = sigma, log = FALSE)
  
  SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = TRUE)
  
  SurvT <- pdffz(SurvT)
  
  Surv_Cure <- nu * SurvT + (1 - nu)
  
  Surv_Cure <- pdffz(Surv_Cure)
  
  
  # derlogpdf_dermu
  derlpdf_derparam <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 
  
  # derSurvivalfunc_dermu
  derSt_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )
  
  deriv_uncens <- derlpdf_derparam * exp(x)
  
  deriv_cens <- (1/Surv_Cure) * nu * derSt_dermu * exp(x)
  
  
  derivs <- delta * deriv_uncens + (1 - delta) * deriv_cens
  
  return(derivs)
}

derLoss_dersigma_func <- function(x, sigma, nu, time){}

derLoss_dernu_func <- function(x, sigma, mu, time){}



Loss_mu_func <- function(x, sigma, nu, time, delta){
  
  
  
  param <- exp(x)

  nu_int <- nu
  
  SurvT_int <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
  
  pdfT_int <- dweibull(x = time, scale = param, shape = sigma)
  
  SurvT_Cure <- nu_int * SurvT_int + (1 - nu_int)
  
  ## Small check
  pdfT_int <- dvffz(pdfT_int)
  #SurvT_int <- pdffz(SurvT_int)
  SurvT_Cure <- pdffz(SurvT_Cure)
  
  lss <- delta * ( log(nu) + log(pdfT_int) ) + (1 - delta) * log(SurvT_Cure)
  
  return(lss)
}

Loss_sigma_func <- function(x, mu, nu, time){
  
  return(lss)
}

Loss_nu_func <- function(x, sigma, mu, time){
  
  return(lss)
}


#derlpdf_dermu_func(x = eta_11, sigma = sigma1)

library(numDeriv)
all.equal(grad(logpdf_mu_func, x = eta_11, sigma = sigma1, time = dat$y), 
          derlpdf_dermu_func(x = eta_11, sigma = sigma1, time = dat$y) 
          )




all.equal(grad(Loss_mu_func, x = eta_11, sigma = sigma1, nu = rho, time = dat$y, delta = dat$status_m1, 
               method = "simple"),
          derLoss_dermu_func(x = eta_11, sigma = sigma1, nu = rho, time = dat$y, delta = dat$status_m1)
          )









TryMSTOP <- 2000

boost.nu.steplength <- 0.1

Cure_formula <- cbind(y, status_m1) ~ .
  
#weights.mstop <- sample(c(0,1), size = nrow(x.train), replace = TRUE, prob = c(0.3, 0.7))

GLM1TRY <- glmboostLSS(formula = Cure_formula, 
                       data = dat, 
                       families = LogNormal_Cure(stabilization = "L2"), 
                       weights = weights.mstop, 
                       method = "noncyclic",
                       control = boost_control(mstop = TryMSTOP, 
                                               nu = boost.nu.steplength, 
                                               risk = "oobag", 
                                               trace = TRUE)) 
plot(risk(GLM1TRY, merge = T))

mstop(GLM1TRY)

coef(GLM1TRY)


#### compute things:
mu_hat <- predict(GLM1TRY$mu, type = "link")

sigma_hat <- exp(predict(GLM1TRY$sigma, type = "link"))

marginal_pdf <- dlnorm(x = dat$y, meanlog = mu_hat, sdlog = sigma_hat)

marginal_Surv <- 1-plnorm(q = dat$y, meanlog = mu_hat, sdlog = sigma_hat)


dat$pdf <- marginal_pdf

dat$Surv <- marginal_Surv

Cure_formula <- cbind(y, status_m1, pdf, Surv) ~ .

CureTRY <- mboost::glmboost(Cure_formula,
                       data = dat, 
                       family = CureFraction_Loss(), 
                       weights = weights.mstop, 
                       control = boost_control(mstop = 5000, 
                                               nu = 0.01, 
                                               risk = "oobag", 
                                               trace = TRUE)) 
coef(CureTRY)





