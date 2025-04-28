Performance.Dyad.Model.3.lag = function(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,
sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML){ 

# Generate data  
  
data = Sim.Dyad.Model.3.lag(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
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

# Compute interactions
data$Female.Z = data$Female*data$Z
data$Male.Z = data$Male*data$Z
data$Y.lag.Z = data$Y.lag*data$Z
data$X.Actor.Z = data$X.Actor*data$Z
data$X.Partner.Z = data$X.Partner*data$Z

# Model estimation

if (is.REML==TRUE){
fit.Model.3 = lme(fixed = Y ~ -1 + Female + Female.Z + Female:Y.lag + Female:Y.lag.Z +
                  Female:X.Actor + Female:X.Actor.Z + Female:X.Partner + Female:X.Partner.Z + 
                  Male + Male.Z +  Male:Y.lag + Male:Y.lag.Z +
                  Male:X.Actor + Male:X.Actor.Z + Male:X.Partner + Male:X.Partner.Z,
                  random = ~ -1 + Female + Male |dyad.ID, 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'REML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  msMaxIter = 100, msVerbose = FALSE))
}

if (is.REML==FALSE){
fit.Model.3 = lme(fixed = Y ~ -1 + Female + Female.Z + Female:Y.lag + Female:Y.lag.Z +
                  Female:X.Actor + Female:X.Actor.Z + Female:X.Partner + Female:X.Partner.Z + 
                  Male + Male.Z +  Male:Y.lag + Male:Y.lag.Z +
                  Male:X.Actor + Male:X.Actor.Z + Male:X.Partner + Male:X.Partner.Z,
                  random = ~ -1 + Female + Male |dyad.ID, 
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
coef = coef(summary(fit.Model.3))[,1]

# Bias
bias = coef(summary(fit.Model.3))[,1] - c(c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,a.FF0,a.FF1,p.MF0,p.MF1,rho.YM0,rho.YM1,a.MM0,a.MM1,p.FM0,p.FM1)

# Power
power = coef(summary(fit.Model.3))[,5]<alpha

# Coverage rate
CI = intervals(fit.Model.3, level = 1-alpha, which = "fixed")$fixed
CI.width = CI[,3] - CI[,1]
CR = CI[,1] < c(c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,a.FF0,a.FF1,p.MF0,p.MF1,rho.YM0,rho.YM1,a.MM0,a.MM1,p.FM0,p.FM1) & CI[,3] > c(c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,a.FF0,a.FF1,p.MF0,p.MF1,rho.YM0,rho.YM1,a.MM0,a.MM1,p.FM0,p.FM1)

summary.fixed.effect = list(coef=coef,bias=bias,power=power,CI.width=CI.width,CR=CR)

# Two-part mixed model: variance components
Sigma.hat = VarCorr(fit.Model.3)
Sigma.weight.hat = 1/unique(varWeights(fit.Model.3$modelStruct$varStruct))
sigma.eps.F.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[1]
sigma.eps.M.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[2]
rho.eps.FM.hat = as.numeric(coef(fit.Model.3$modelStruct$corStruct,unconstrained=FALSE)) 
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
rho.eps.FM.hat=rho.eps.FM.hat,sigma.nu.F.hat=sigma.nu.F.hat,sigma.nu.M.hat=sigma.nu.M.hat,
rho.nu.F.M.hat=rho.nu.F.M.hat)

summary.var.bias = c(sigma.eps.F.bias=sigma.eps.F.bias,
sigma.eps.M.bias=sigma.eps.M.bias,rho.eps.FM.bias=rho.eps.FM.bias,
sigma.nu.F.bias=sigma.nu.F.bias,sigma.nu.M.bias=sigma.nu.M.bias,
rho.nu.F.M.bias=rho.nu.F.M.bias)

return(list(summary.fixed.effect=summary.fixed.effect,summary.var.hat=summary.var.hat,summary.var.bias=summary.var.bias))
}