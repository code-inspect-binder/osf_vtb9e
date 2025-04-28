Performance.Dyad.Model.15.lag = function(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML){ 

# Generate data  
  
data = Sim.Dyad.Model.15.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)

# Fit multilevel model

if(is.center.X==TRUE){
# Person-centered the predictors
data <- data %>% 
group_by(subject.ID,dyad.ID) %>% 
mutate(X.Actor = X.Actor - mean(X.Actor),
       X.Partner = X.Partner - mean(X.Partner))
}

# Compute quadratic variables
data$X.Actor.2 = I(data$X.Actor^2)
data$X.Partner.2 = I(data$X.Partner^2)

# Compute interactions
data$X.Actor.D = data$X.Actor*data$D
data$X.Partner.D = data$X.Partner*data$D
data$X.Actor.D.2 = data$X.Actor.2*data$D
data$X.Partner.D.2 = data$X.Partner.2*data$D

# Model estimation

if (is.REML==TRUE){
fit.Model.15 = lme(fixed = Y ~ -1 + Female + Female:Y.lag + 
                  Female:X.Actor + Female:X.Actor.2 +
                  Female:X.Partner + Female:X.Partner.2 + Female:D +
                  Female:X.Actor.D + Female:X.Actor.D.2 + 
                  Female:X.Partner.D + Female:X.Partner.D.2 + 
                  Male + Male:Y.lag + 
                  Male:X.Actor + Male:X.Actor.2 + Male:X.Partner +
                  Male:X.Partner.2 + Male:D + 
                  Male:X.Actor.D + Male:X.Actor.D.2 + Male:X.Partner.D + 
                  Male:X.Partner.D.2,
                  random = ~ -1 + Female + Male|dyad.ID, 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'REML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  maxIter=500, msMaxIter=500, msVerbose = FALSE))
}

if (is.REML==FALSE){
fit.Model.15 = lme(fixed = Y ~ -1 + Female + Female:Y.lag + 
                  Female:X.Actor + Female:X.Actor.2 +
                  Female:X.Partner + Female:X.Partner.2 + Female:D +
                  Female:X.Actor.D + Female:X.Actor.D.2 + 
                  Female:X.Partner.D + Female:X.Partner.D.2 + 
                  Male + Male:Y.lag + 
                  Male:X.Actor + Male:X.Actor.2 + Male:X.Partner +
                  Male:X.Partner.2 + Male:D + 
                  Male:X.Actor.D + Male:X.Actor.D.2 + Male:X.Partner.D + 
                  Male:X.Partner.D.2,
                  random = ~ -1 + Female + Male|dyad.ID, 
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
coef = coef(summary(fit.Model.15))[,1]

# Bias
bias = coef(summary(fit.Model.15))[,1] - c(c.F,c.M,rho.YF,a.FF,a.FF2,p.MF,p.MF2,d.F,d.FF,d.FF2,d.MF,d.MF2,rho.YM,a.MM,a.MM2,p.FM,p.FM2,d.M,d.MM,d.MM2,d.FM,d.FM2)

# Power
power = coef(summary(fit.Model.15))[,5]<alpha

# Coverage rate
CI = intervals(fit.Model.15, level = 1-alpha, which = "fixed")$fixed
CI.width = CI[,3] - CI[,1]
CR = CI[,1] < c(c.F,c.M,rho.YF,a.FF,a.FF2,p.MF,p.MF2,d.F,d.FF,d.FF2,d.MF,d.MF2,rho.YM,a.MM,a.MM2,p.FM,p.FM2,d.M,d.MM,d.MM2,d.FM,d.FM2) & CI[,3] > c(c.F,c.M,rho.YF,a.FF,a.FF2,p.MF,p.MF2,d.F,d.FF,d.FF2,d.MF,d.MF2,rho.YM,a.MM,a.MM2,p.FM,p.FM2,d.M,d.MM,d.MM2,d.FM,d.FM2)

summary.fixed.effect = list(coef=coef,bias=bias,power=power,CI.width=CI.width,CR=CR)

# Variance components
Sigma.hat = VarCorr(fit.Model.15)
Sigma.weight.hat = 1/unique(varWeights(fit.Model.15$modelStruct$varStruct))
sigma.eps.F.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[1]
sigma.eps.M.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[2]
rho.eps.FM.hat = as.numeric(coef(fit.Model.15$modelStruct$corStruct,unconstrained=FALSE)) 
sigma.nu.F.hat = as.numeric(Sigma.hat['Female',2])
sigma.nu.M.hat = as.numeric(Sigma.hat['Male',2])
rho.nu.F.M.hat = as.numeric(Sigma.hat['Male',3])

sigma.eps.F.bias = sigma.eps.F.hat - sigma.eps.F
sigma.eps.M.bias = sigma.eps.M.hat - sigma.eps.M
rho.eps.FM.bias = rho.eps.FM.hat - rho.eps.FM 
sigma.nu.F.bias = sigma.nu.F.hat - sigma.nu.F
sigma.nu.M.bias = sigma.nu.M.hat - sigma.nu.M
rho.nu.F.M.bias = rho.nu.F.M.hat - rho.nu.F.M

summary.var.hat = c(sigma.eps.F.hat=sigma.eps.F.hat,sigma.eps.M.hat=sigma.eps.M.hat,
rho.eps.FM.hat=rho.eps.FM.hat,sigma.nu.F.hat=sigma.nu.F.hat,sigma.nu.M.hat=sigma.nu.M.hat,rho.nu.F.M.hat=rho.nu.F.M.hat)

summary.var.bias = c(sigma.eps.F.bias=sigma.eps.F.bias,
sigma.eps.M.bias=sigma.eps.M.bias,rho.eps.FM.bias=rho.eps.FM.bias,
sigma.nu.F.bias=sigma.nu.F.bias,sigma.nu.M.bias=sigma.nu.M.bias,
rho.nu.F.M.bias=rho.nu.F.M.bias)

return(list(summary.fixed.effect=summary.fixed.effect,summary.var.hat=summary.var.hat,summary.var.bias=summary.var.bias))
}