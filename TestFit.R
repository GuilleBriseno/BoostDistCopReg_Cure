#### Try a two stage approach: 


TryMSTOP <- 2000

boost.nu.steplength <- 0.1

Cure_formula <- cbind(y, status_m1) ~ .

#weights.mstop <- sample(c(0,1), size = nrow(x.train), replace = TRUE, prob = c(0.3, 0.7))

GLM1TRY <- glmboostLSS(formula = Cure_formula, 
                       data = dat, 
                       families = LogNormal_Cure(stabilization = "none"), 
                       weights = weights.mstop, 
                       method = "noncyclic",
                       control = boost_control(mstop = TryMSTOP, 
                                               nu = boost.nu.steplength, 
                                               risk = "oobag", 
                                               trace = TRUE)) 
plot(risk(GLM1TRY, merge = T))

mstop(GLM1TRY)

coef(GLM1TRY)
