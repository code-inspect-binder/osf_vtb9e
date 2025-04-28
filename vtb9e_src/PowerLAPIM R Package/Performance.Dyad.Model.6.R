Performance.Dyad.Model.6 = function(N.dyad,T.obs,  
c,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML){ 

# Generate data  
  
data = Sim.Dyad.Model.6(N.dyad,T.obs,  
c,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
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
fit.Model.6 = lme(fixed = Y ~ 1 + X.Actor + X.Partner + W +                        
                  X.Actor.W + X.Partner.W,
                  random = list(dyad.ID = pdCompSymm(~ Gender -1)),
                  correlation = corCompSymm(form = ~1|dyad.ID/Obs),
                  weights = varIdent(form = ~1|Gender),
                  data = data, na.action=na.omit, 
                  method = 'REML',
                  control = lmeControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), 
                  maxIter=500, msMaxIter=500, msVerbose = FALSE))  
}

if (is.REML==FALSE){
fit.Model.6 = lme(fixed = Y ~ 1 + X.Actor + X.Partner + W +                        
                  X.Actor.W + X.Partner.W,
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
coef = coef(summary(fit.Model.6))[,1]

# Bias
bias = coef(summary(fit.Model.6))[,1] - c(c,a,p,b,b.a,b.p)

# Power
power = coef(summary(fit.Model.6))[,5]<alpha

# Coverage rate
CI = intervals(fit.Model.6, level = 1-alpha, which = "fixed")$fixed
CI.width = CI[,3] - CI[,1]
CR = CI[,1] < c(c,a,p,b,b.a,b.p) & CI[,3] > c(c,a,p,b,b.a,b.p)

summary.fixed.effect = list(coef=coef,bias=bias,power=power,CI.width=CI.width,CR=CR)

# Variance components
Sigma.hat = VarCorr(fit.Model.6)
Sigma.weight.hat = 1/unique(varWeights(fit.Model.6$modelStruct$varStruct))
sigma.eps.F.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[1]
sigma.eps.M.hat = as.numeric(Sigma.hat['Residual',2])*Sigma.weight.hat[2]
rho.eps.FM.hat = as.numeric(coef(fit.Model.6$modelStruct$corStruct,unconstrained=FALSE)) 
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