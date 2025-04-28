Performance.Dyad.Model.2.lag = function(N.dyad,T.obs,  
c,rho.Y,a,p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML){ 

# Generate data  
  
data = Sim.Dyad.Model.2.lag(N.dyad,T.obs,  
c,rho.Y,a,p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)

# Fit multilevel model

if(is.center.X==TRUE){
# Person-centered the predictors
data <- data %>% 
group_by(subject.ID,dyad.ID) %>% 
mutate(X.Actor = X.Actor - mean(X.Actor),
       X.Partner = X.Partner - mean(X.Partner))
}

# Model estimation

if (is.REML==TRUE){
fit.Model.2 = lme(fixed = Y ~ 1 + Y.lag + X.Actor + X.Partner,
                  random = list(dyad.ID = pdCompSymm(~ Gender -1)), 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'REML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  msMaxIter = 100, msVerbose = FALSE))
}

if (is.REML==FALSE){
fit.Model.2 = lme(fixed = Y ~ 1 + Y.lag + X.Actor + X.Partner,
                  random = list(dyad.ID = pdCompSymm(~ Gender -1)), 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'ML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  msMaxIter = 100, msVerbose = FALSE))
}

# Performance measures 

# Fixed Effects
# Estimated values
coef = coef(summary(fit.Model.2))[,1]

# Bias
bias = coef(summary(fit.Model.2))[,1] - c(c,rho.Y,a,p)

# Power
power = coef(summary(fit.Model.2))[,5]<alpha

# Coverage rate
CI = intervals(fit.Model.2, level = 1-alpha, which = "fixed")$fixed
CI.width = CI[,3] - CI[,1]
CR = CI[,1] < c(c,rho.Y,a,p) & CI[,3] > c(c,rho.Y,a,p)

summary.fixed.effect = list(coef=coef,bias=bias,power=power,CI.width=CI.width,CR=CR)

# Two-part mixed model: variance components
Sigma.hat = VarCorr(fit.Model.2)
Sigma.weight.hat = 1/unique(varWeights(fit.Model.2$modelStruct$varStruct))
sigma.eps.F.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[1]
sigma.eps.M.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[2]
rho.eps.FM.hat = as.numeric(coef(fit.Model.2$modelStruct$corStruct,unconstrained=FALSE)) 
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