Performance.Dyad.Model.5.lag = function(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML){ 

# Generate data  
  
data = Sim.Dyad.Model.5.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)

# Fit multilevel model

if(is.center.X==TRUE){
# Person-centered the predictors
data <- data %>% 
group_by(subject.ID,dyad.ID) %>% 
mutate(X.Actor = X.Actor - mean(X.Actor),
       X.Partner = X.Partner - mean(X.Partner))
}

if(is.center.W==TRUE){
# Person-centered the predictors
data <- data %>% 
group_by(subject.ID,dyad.ID) %>% 
mutate(W = W - mean(W))
}

# Compute interactions
data$X.Actor.W = data$X.Actor*data$W
data$X.Partner.W = data$X.Partner*data$W

# Model estimation

if (is.REML==TRUE){
fit.Model.5 = lme(fixed = Y ~ -1 + Female + Female:Y.lag +
                  Female:X.Actor + Female:X.Partner +
                  Female:W + Female:X.Actor.W + Female:X.Partner.W + 
                  Male + Male:Y.lag +  
                  Male:X.Actor + Male:X.Partner + Male:W + 
                  Male:X.Actor.W + Male:X.Partner.W,
                  random = ~ -1 + Female + Male|dyad.ID, 
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'REML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  maxIter=500, msMaxIter=500, msVerbose = FALSE)) 
}

if (is.REML==FALSE){
fit.Model.5 = lme(fixed = Y ~ -1 + Female + Female:Y.lag +
                  Female:X.Actor + Female:X.Partner +
                  Female:W + Female:X.Actor.W + Female:X.Partner.W + 
                  Male + Male:Y.lag +  
                  Male:X.Actor + Male:X.Partner + Male:W + 
                  Male:X.Actor.W + Male:X.Partner.W,
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
coef = coef(summary(fit.Model.5))[,1]

# Bias
bias = coef(summary(fit.Model.5))[,1] - c(c.F,c.M,rho.YF,a.FF,p.MF,b.F,b.FF,b.MF,rho.YM,a.MM,p.FM,b.M,b.MM,b.FM)

# Power
power = coef(summary(fit.Model.5))[,5]<alpha

# Coverage rate
CI = intervals(fit.Model.5, level = 1-alpha, which = "fixed")$fixed
CI.width = CI[,3] - CI[,1]
CR = CI[,1] < c(c.F,c.M,rho.YF,a.FF,p.MF,b.F,b.FF,b.MF,rho.YM,a.MM,p.FM,b.M,b.MM,b.FM) & CI[,3] > c(c.F,c.M,rho.YF,a.FF,p.MF,b.F,b.FF,b.MF,rho.YM,a.MM,p.FM,b.M,b.MM,b.FM)

summary.fixed.effect = list(coef=coef,bias=bias,power=power,CI.width=CI.width,CR=CR)

# Variance components
Sigma.hat = VarCorr(fit.Model.5)
Sigma.weight.hat = 1/unique(varWeights(fit.Model.5$modelStruct$varStruct))
sigma.eps.F.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[1]
sigma.eps.M.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[2]
rho.eps.FM.hat = as.numeric(coef(fit.Model.5$modelStruct$corStruct,unconstrained=FALSE)) 
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