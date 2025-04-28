Performance.Dyad.Model.12.lag = function(N0.dyad,N1.dyad,T.obs,  
c0,c1,rho.Y0,rho.Y1,a0,a1,a02,a12,p0,p1,p02,p12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,
sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML){ 

# Generate data  
  
data = Sim.Dyad.Model.12.lag(N0.dyad,N1.dyad,T.obs,  
c0,c1,rho.Y0,rho,Y1,a0,a1,a02,a12,p0,p1,p02,p12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,
sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)

# Fit multilevel model

if(is.center.X==TRUE){
# Person-centered the predictors
data <- data %>% 
group_by(subject.ID,dyad.ID) %>% 
mutate(X.Actor = X.Actor - mean(X.Actor),
       X.Partner = X.Partner - mean(X.Partner))
}

# Compute interaction with lagged variable
data$Y.lag.Z = data$Y.lag*data$Z

# Compute quadratic effects
data$X.Actor.2 = I(data$X.Actor^2)
data$X.Partner.2 = I(data$X.Partner^2)

# Compute interactions
data$X.Actor.Z = data$X.Actor*data$Z
data$X.Partner.Z = data$X.Partner*data$Z

data$X.Actor.2.Z = data$X.Actor.2*data$Z
data$X.Partner.2.Z = data$X.Partner.2*data$Z


# Model estimation

if (is.REML==TRUE){
fit.Model.12 = lme(fixed = Y ~ 1 + Z + Y.lag + Y.lag.Z +
                  X.Actor + X.Actor.2 + X.Actor.Z + X.Actor.2.Z + 
                  X.Partner + X.Partner.2 + X.Partner.Z + X.Partner.2.Z,
                  random = list(dyad.ID = pdCompSymm(~ Gender -1)), 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'REML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  maxIter=500, msMaxIter=500, msVerbose = FALSE))
}

if (is.REML==FALSE){
fit.Model.12 = lme(fixed = Y ~ 1 + Z + Y.lag + Y.lag.Z +
                  X.Actor + X.Actor.2 + X.Actor.Z + X.Actor.2.Z + 
                  X.Partner + X.Partner.2 + X.Partner.Z + X.Partner.2.Z,
                  random = list(dyad.ID = pdCompSymm(~ Gender -1)), 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'ML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  maxIter=500, msMaxIter=500, msVerbose = FALSE))
}

# Performance measures

# Fixed Effects
# Estimated values
coef = coef(summary(fit.Model.12))[,1]

# Bias
bias = coef(summary(fit.Model.12))[,1] - c(c0,c1,rho.Y0,rho.Y1,a0,a02,a1,a12,p0,p02,p1,p12)

# Power
power = coef(summary(fit.Model.12))[,5]<alpha

# Coverage rate
CI = intervals(fit.Model.12, level = 1-alpha, which = "fixed")$fixed
CI.width = CI[,3] - CI[,1]
CR = CI[,1] < c(c0,c1,rho.Y0,rho.Y1,a0,a02,a1,a12,p0,p02,p1,p12) & CI[,3] > c(c0,c1,rho.Y0,rho.Y1,a0,a02,a1,a12,p0,p02,p1,p12)

summary.fixed.effect = list(coef=coef,bias=bias,power=power,CI.width=CI.width,CR=CR)

# Variance components
Sigma.hat = VarCorr(fit.Model.12)
Sigma.weight.hat = 1/unique(varWeights(fit.Model.12$modelStruct$varStruct))
sigma.eps.F.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[1]
sigma.eps.M.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[2]
rho.eps.FM.hat = as.numeric(coef(fit.Model.12$modelStruct$corStruct,unconstrained=FALSE)) 
sigma.nu.hat = as.numeric(Sigma.hat['GenderF',2])

sigma.eps.F.bias = sigma.eps.F.hat - sigma.eps.F
sigma.eps.M.bias = sigma.eps.M.hat - sigma.eps.M
rho.eps.FM.bias = rho.eps.FM.hat - rho.eps.FM 
sigma.nu.bias = sigma.nu.hat - sigma.nu

summary.var.hat = c(sigma.eps.F.hat=sigma.eps.F.hat,sigma.eps.M.hat=sigma.eps.M.hat,
rho.eps.FM.hat=rho.eps.FM.hat,sigma.nu.hat=sigma.nu.hat)

summary.var.bias = c(sigma.eps.F.bias=sigma.eps.F.bias,
sigma.eps.M.bias=sigma.eps.M.bias,rho.eps.FM.bias=rho.eps.FM.bias,
sigma.nu.bias=sigma.nu.bias)

return(list(summary.fixed.effect=summary.fixed.effect,summary.var.hat=summary.var.hat,summary.var.bias=summary.var.bias))
}