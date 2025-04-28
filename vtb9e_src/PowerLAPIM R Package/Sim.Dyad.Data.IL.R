###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to simulate an IL data set

Sim.Dyad.Data.IL = function(Model,N.dyad,N0.dyad,N1.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.FF2,p.MF2,a.MM,p.FM,a.MM2,p.FM2,
c,a,a.2,p,p.2,
c.F0,c.F1,c.M0,c.M1,a.FF0,a.FF1,a.FF02,a.FF12,p.MF0,p.MF1,p.MF02,p.MF12,
a.MM0,a.MM1,a.MM02,a.MM12,p.FM0,p.FM1,p.FM02,p.FM12,
c0,c1,a0,a1,a02,a12,p0,p1,p02,p12,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,b.FF2,b.MF2,b.MM2,b.FM2,
d.F,d.M,d.FF,d.MF,d.MM,d.FM,d.FF2,d.MF2,d.MM2,d.FM2,
b,b.a,b.a2,b.p,b.p2,
d,d.a,d.a2,d.p,d.p2,
rho.YF,rho.YM,rho.Y,rho.YF0,rho.YF1,rho.YM0,rho.YM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,
rho.X0,rho.X1,
mu.W,sigma.W,prob.D,
is.center.X,is.center.W){

# Check parameters & and output messages

if (length(T.obs) == 0) {stop('Set the number of time points')}
if (sigma.eps.F<=0) {stop('The standard deviation of the Level 1 partner A error must be positive')}
if (length(sigma.eps.F)==0) {stop('The standard deviation of the Level 1 partner A error must be positive')}
if (sigma.eps.M<=0) {stop('The standard deviation of the Level 1 partner B error must be positive')}
if (length(sigma.eps.M)==0) {stop('The standard deviation of the Level 1 partner B error must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 1 partner A and M errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 1 partner A and M errors must be included in the interval [-1,1]')}

########################################################################################

# Simulate data from APIM model

if (Model == 1 | Model == 17){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation between the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation between the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(c.F) == 0) {stop('Set the value of the fixed intercept for partner A')}
if (length(c.M) == 0) {stop('Set the value of the fixed intercept for partner B')}
if (length(a.FF) == 0) {stop('Set the value of the fixed actor effect for partner A')}
if (length(p.MF) == 0) {stop('Set the value of the fixed partner effect for partner A')}
if (length(a.MM) == 0) {stop('Set the value of the fixed actor effect for partner B')}
if (length(p.FM) == 0) {stop('Set the value of the fixed partner effect for partner B')}

if (Model == 1){
data = Sim.Dyad.Model.1(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}

if (Model == 17){
if (abs(rho.YF)>1) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (length(rho.YF)==0) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (abs(rho.YM)>1) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}
if (length(rho.YM)==0) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.1.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}}

########################################################################################

# Simulate data from APIM model with indistinguishable dyads 

if (Model == 2 | Model == 18){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed partner effect')}

if (Model == 2){
data = Sim.Dyad.Model.2(N.dyad,T.obs,  
c,a,p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}

if (Model == 18){
if (abs(rho.Y)>1) {stop('The autoregressive effect must be included in the interval [-1,1]')}
if (length(rho.Y)==0) {stop('The autoregressive effect must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.2.lag(N.dyad,T.obs,  
c,rho.Y,a,p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}}

########################################################################################

# Simulate data from APIM model: group differences in actor partner effects

if (Model == 3 | Model == 19){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF0) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 0')}
if (length(mu.XM0) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 0')}
if (length(mu.XF1) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 1')}
if (length(mu.XM1) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 1')}
if (sigma.XF0<=0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XF0)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (sigma.XM0<=0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XM0)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (sigma.XF0<=1) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XF1)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (sigma.XM0<=1) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XM1)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (abs(rho.X0)>1) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (length(rho.X0)==0) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.X1)>1) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(rho.X1)==0) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(c.F0) == 0) {stop('Set the value of the fixed intercept for partner A in Group 0')}
if (length(c.M0) == 0) {stop('Set the value of the fixed intercept for partner B in Group 0')}
if (length(a.FF0) == 0) {stop('Set the value of the fixed actor effect for partner A in Group 0')}
if (length(p.MF0) == 0) {stop('Set the value of the fixed partner effect for partner A in Group 0')}
if (length(a.MM0) == 0) {stop('Set the value of the fixed actor effect for partner B in Group 0')}
if (length(p.FM0) == 0) {stop('Set the value of the fixed partner effect for partner B in Group 0')}
if (length(c.F1) == 0) {stop('Set the value of the fixed intercept for partner A in Group 1')}
if (length(c.M1) == 0) {stop('Set the value of the fixed intercept for partner B in Group 1')}
if (length(a.FF1) == 0) {stop('Set the value of the fixed actor effect for partner A in Group 1')}
if (length(p.MF1) == 0) {stop('Set the value of the fixed partner effect for partner A in Group 1')}
if (length(a.MM1) == 0) {stop('Set the value of the fixed actor effect for partner B in Group 1')}
if (length(p.FM1) == 0) {stop('Set the value of the fixed partner effect for partner B in Group 1')}

if (Model == 3){
data = Sim.Dyad.Model.3(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}

if (Model == 19){
if (abs(rho.YF0)>1) {stop('The autoregressive effect for partner A in Group 0 must be included in the interval [-1,1]')}
if (length(rho.YF0)==0) {stop('The autoregressive effect for partner A in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.YF0+rho.YF1)>1) {stop('The autoregressive effect for partner A in Group 1 must be included in the interval [-1,1]')}
if (length(rho.YF1)==0) {stop('The autoregressive effect for partner A in Group 1 must be included in the interval [-1,1]')}
if (abs(rho.YM0)>1) {stop('The autoregressive effect for partner B in Group 0 must be included in the interval [-1,1]')}
if (length(rho.YM0)==0) {stop('The autoregressive effect for partner B in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.YM0+rho.YM1)>1) {stop('The autoregressive effect for partner B in Group 1 must be included in the interval [-1,1]')}
if (length(rho.YM1)==0) {stop('The autoregressive effect for partner B in Group 1 must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.3.lag(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}}

########################################################################################

# Simulate data from APIM model: group differences in actor partner effects
# with indistinguishable dyads 

if (Model == 4 | Model == 20){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF0) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 0')}
if (length(mu.XM0) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 0')}
if (length(mu.XF1) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 1')}
if (length(mu.XM1) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 1')}
if (sigma.XF0<=0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XF0)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (sigma.XM0<=0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XM0)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (sigma.XF0<=1) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XF1)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (sigma.XM0<=1) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XM1)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (abs(rho.X0)>1) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (length(rho.X0)==0) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.X1)>1) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(rho.X1)==0) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(c0) == 0) {stop('Set the value of the fixed intercept in Group 0')}
if (length(a0) == 0) {stop('Set the value of the fixed actor effect in Group 0')}
if (length(p0) == 0) {stop('Set the value of the fixed partner effect in Group 0')}
if (length(c1) == 0) {stop('Set the value of the fixed intercept in Group 1')}
if (length(a1) == 0) {stop('Set the value of the fixed actor effect in Group 1')}
if (length(p1) == 0) {stop('Set the value of the fixed partner effect in Group 1')}

if (Model == 4){
data = Sim.Dyad.Model.4(N0.dyad,N1.dyad,T.obs,  
c0,c1,a0,a1,p0,p1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}

if (Model == 20){
if (abs(rho.Y0)>1) {stop('The autoregressive effect in Group 0 must be included in the interval [-1,1]')}
if (length(rho.Y0)==0) {stop('The autoregressive effect in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.Y0+rho.Y1)>1) {stop('The autoregressive effect in Group 1 must be included in the interval [-1,1]')}
if (length(rho.Y1)==0) {stop('The autoregressive effect in Group 1 must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.4.lag(N0.dyad,N1.dyad,T.obs,  
c0,c1,rho.Y0,rho.Y1,a0,a1,p0,p1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}}

########################################################################################

# Simulate data from APIM model with a continuos time-varying dyad moderator

if (Model == 5 | Model == 21){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(mu.W) == 0) {stop('Set the value of the mean of the Level 1 continuous dyad moderator')}
if (sigma.W==0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(sigma.W)==0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(c.F) == 0) {stop('Set the value of the fixed intercept for partner A')}
if (length(c.M) == 0) {stop('Set the value of the fixed intercept for partner B')}
if (length(a.FF) == 0) {stop('Set the value of the fixed actor effect for partner A')}
if (length(p.MF) == 0) {stop('Set the value of the fixed partner effect for partner A')}
if (length(a.MM) == 0) {stop('Set the value of the fixed actor effect for partner B')}
if (length(p.FM) == 0) {stop('Set the value of the fixed partner effect for partner B')}
if (length(b.F) == 0) {stop('Set the value of the fixed effect of the Level 1 continuous dyad moderator for partner A')}
if (length(b.M) == 0) {stop('Set the value of the fixed effect of the Level 1 continuous dyad moderator for partner B')}
if (length(b.FF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor effect for partner A')}
if (length(b.MF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner effect for partner A')}
if (length(b.MM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor effect for partner B')}
if (length(b.FM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner effect for partner B')}

if (Model == 5){
data = Sim.Dyad.Model.5(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}

if (Model == 21){
if (abs(rho.YF)>1) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (length(rho.YF)==0) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (abs(rho.YM)>1) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}
if (length(rho.YM)==0) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.5.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}}

########################################################################################

# Simulate data from APIM model with a continuos time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 6 | Model == 22){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(mu.W) == 0) {stop('Set the value of the mean of the Level 1 continuous dyad moderator')}
if (sigma.W == 0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(sigma.W)==0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed partner effect')}
if (length(b) == 0) {stop('Set the value of the fixed effect of the Level 1 continuous dyad moderator')}
if (length(b.a) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor effect')}
if (length(b.p) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner effect')}

if (Model == 6){
data = Sim.Dyad.Model.6(N.dyad,T.obs,  
c,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}

if (Model == 22){
if (abs(rho.Y)>1) {stop('The autoregressive effect must be included in the interval [-1,1]')}
if (length(rho.Y)==0) {stop('The autoregressive effect must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.6.lag(N.dyad,T.obs,  
c,rho.Y,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}}

########################################################################################

# Simulate data from APIM model with a dichotonomous time-varying dyad moderator

if (Model == 7 | Model == 23){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (abs(prob.D)>1) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(prob.D)==0) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(c.F) == 0) {stop('Set the value of the fixed intercept for partner A')}
if (length(c.M) == 0) {stop('Set the value of the fixed intercept for partner B')}
if (length(a.FF) == 0) {stop('Set the value of the fixed actor effect for partner A')}
if (length(p.MF) == 0) {stop('Set the value of the fixed partner effect for partner A')}
if (length(a.MM) == 0) {stop('Set the value of the fixed actor effect for partner B')}
if (length(p.FM) == 0) {stop('Set the value of the fixed partner effect for partner B')}
if (length(d.F) == 0) {stop('Set the value of the fixed effect of the Level 1 dichotomous dyad moderator for partner A')}
if (length(d.M) == 0) {stop('Set the value of the fixed effect of the Level 1 dichotomous dyad moderator for partner B')}
if (length(d.FF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor effect for partner A')}
if (length(d.MF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner effect for partner A')}
if (length(d.MM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor effect for partner B')}
if (length(d.FM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner effect for partner B')}

if (Model == 7){
data = Sim.Dyad.Model.7(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
d.F,d.M,d.FF,d.MF,d.MM,d.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}

if (Model == 23){
if (abs(rho.YF)>1) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (length(rho.YF)==0) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (abs(rho.YM)>1) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}
if (length(rho.YM)==0) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.7.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
d.F,d.M,d.FF,d.MF,d.MM,d.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}}

########################################################################################

# Simulate data from APIM model with a dichotonomous time-varying dyad moderator
# with indistinguishable dyads

if (Model == 8 | Model == 24){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed partner effect')}
if (abs(prob.D)>1) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(prob.D)==0) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed partner effect')}
if (length(d) == 0) {stop('Set the value of the fixed effect of the Level 1 dichotomous dyad moderator')}
if (length(d.a) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor effect')}
if (length(d.p) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner effect')}

if (Model == 8){
data = Sim.Dyad.Model.8(N.dyad,T.obs,  
c,a,p,d,d.a,d.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}

if (Model == 24){
if (abs(rho.Y)>1) {stop('The autoregressive effect must be included in the interval [-1,1]')}
if (length(rho.Y)==0) {stop('The autoregressive effect must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.8.lag(N.dyad,T.obs,  
c,rho.Y,a,p,d,d.a,d.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model 

if (Model == 9 | Model == 25){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(c.F) == 0) {stop('Set the value of the fixed intercept for partner A')}
if (length(c.M) == 0) {stop('Set the value of the fixed intercept for partner B')}
if (length(a.FF) == 0) {stop('Set the value of the fixed linear actor effect for partner A')}
if (length(p.MF) == 0) {stop('Set the value of the fixed linear partner effect for partner A')}
if (length(a.MM) == 0) {stop('Set the value of the fixed linear actor effect for partner B')}
if (length(p.FM) == 0) {stop('Set the value of the fixed linear partner effect for partner B')}
if (length(a.FF2) == 0) {stop('Set the value of the fixed quadratic actor effect for partner A')}
if (length(p.MF2) == 0) {stop('Set the value of the fixed quadratic partner effect for partner A')}
if (length(a.MM2) == 0) {stop('Set the value of the fixed quadratic actor effect for partner B')}
if (length(p.FM2) == 0) {stop('Set the value of the fixed quadratic partner effect for partner B')}

if (Model == 9){
data = Sim.Dyad.Model.9(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.FF2,p.MF2,a.MM,p.FM,a.MM2,p.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}

if (Model == 25){
if (abs(rho.YF)>1) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (length(rho.YF)==0) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (abs(rho.YM)>1) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}
if (length(rho.YM)==0) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.9.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.FF2,p.MF2,a.MM,p.FM,a.MM2,p.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with indistinguishable dyads 

if (Model == 10 | Model == 26){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed linear actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed linear partner effect')}
if (length(a.2) == 0) {stop('Set the value of the fixed quadratic actor effect')}
if (length(p.2) == 0) {stop('Set the value of the fixed quadratic partner effect')}

if (Model == 10){
data = Sim.Dyad.Model.10(N.dyad,T.obs,  
c,a,a.2,p,p.2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}

if (Model == 26){
if (abs(rho.Y)>1) {stop('The autoregressive effect must be included in the interval [-1,1]')}
if (length(rho.Y)==0) {stop('The autoregressive effect must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.10.lag(N.dyad,T.obs,  
c,rho.Y,a,a.2,p,p.2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model: group differences in actor partner effects

if (Model == 11 | Model == 27){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF0) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 0')}
if (length(mu.XM0) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 0')}
if (length(mu.XF1) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 1')}
if (length(mu.XM1) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 1')}
if (sigma.XF0<=0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XF0)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (sigma.XM0<=0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XM0)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (sigma.XF0<=1) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XF1)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (sigma.XM0<=1) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XM1)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (abs(rho.X0)>1) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (length(rho.X0)==0) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.X1)>1) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(rho.X1)==0) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(c.F0) == 0) {stop('Set the value of the fixed intercept for partner A in Group 0')}
if (length(c.M0) == 0) {stop('Set the value of the fixed intercept for partner B in Group 0')}
if (length(a.FF0) == 0) {stop('Set the value of the fixed actor linear effect for partner A in Group 0')}
if (length(p.MF0) == 0) {stop('Set the value of the fixed partner linear effect for partner A in Group 0')}
if (length(a.MM0) == 0) {stop('Set the value of the fixed actor linear effect for partner B in Group 0')}
if (length(p.FM0) == 0) {stop('Set the value of the fixed partner linear effect for partner B in Group 0')}
if (length(a.FF02) == 0) {stop('Set the value of the fixed actor quadratic effect for partner A in Group 0')}
if (length(p.MF02) == 0) {stop('Set the value of the fixed partner quadratic effect for partner A in Group 0')}
if (length(a.MM02) == 0) {stop('Set the value of the fixed actor quadratic effect for partner B in Group 0')}
if (length(p.FM02) == 0) {stop('Set the value of the fixed partner quadratic effect for partner B in Group 0')}
if (length(c.F1) == 0) {stop('Set the value of the fixed intercept for partner A in Group 1')}
if (length(c.M1) == 0) {stop('Set the value of the fixed intercept for partner B in Group 1')}
if (length(a.FF1) == 0) {stop('Set the value of the fixed actor linear effect for partner A in Group 1')}
if (length(p.MF1) == 0) {stop('Set the value of the fixed partner linear effect for partner A in Group 1')}
if (length(a.MM1) == 0) {stop('Set the value of the fixed actor linear effect for partner B in Group 1')}
if (length(p.FM1) == 0) {stop('Set the value of the fixed partner linear effect for partner B in Group 1')}
if (length(a.FF12) == 0) {stop('Set the value of the fixed actor quadratic effect for partner A in Group 1')}
if (length(p.MF12) == 0) {stop('Set the value of the fixed partner quadratic effect for partner A in Group 1')}
if (length(a.MM12) == 0) {stop('Set the value of the fixed actor quadratic effect for partner B in Group 1')}
if (length(p.FM12) == 0) {stop('Set the value of the fixed partner quadratic effect for partner B in Group 1')}

if (Model == 11){
data = Sim.Dyad.Model.11(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,a.FF0,a.FF1,a.FF02,a.FF12,p.MF0,p.MF1,p.MF02,p.MF12,
a.MM0,a.MM1,a.MM02,a.MM12,p.FM0,p.FM1,p.FM02,p.FM12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}

if (Model == 27){
if (abs(rho.YF0)>1) {stop('The autoregressive effect for partner A in Group 0 must be included in the interval [-1,1]')}
if (length(rho.YF0)==0) {stop('The autoregressive effect for partner A in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.YF0+rho.YF1)>1) {stop('The autoregressive effect for partner A in Group 1 must be included in the interval [-1,1]')}
if (length(rho.YF1)==0) {stop('The autoregressive effect for partner A in Group 1 must be included in the interval [-1,1]')}
if (abs(rho.YM0)>1) {stop('The autoregressive effect for partner B in Group 0 must be included in the interval [-1,1]')}
if (length(rho.YM0)==0) {stop('The autoregressive effect for partner B in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.YM0+rho.YM1)>1) {stop('The autoregressive effect for partner B in Group 1 must be included in the interval [-1,1]')}
if (length(rho.YM1)==0) {stop('The autoregressive effect for partner B in Group 1 must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.11.lag(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,a.FF0,a.FF1,a.FF02,a.FF12,p.MF0,p.MF1,p.MF02,p.MF12,
a.MM0,a.MM1,a.MM02,a.MM12,p.FM0,p.FM1,p.FM02,p.FM12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}}


########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model: group differences in actor partner effects
# with indistinguishable dyads 

if (Model == 12 | Model == 28){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF0) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 0')}
if (length(mu.XM0) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 0')}
if (length(mu.XF1) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor in Group 1')}
if (length(mu.XM1) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor in Group 1')}
if (sigma.XF0<=0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XF0)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 0 must be positive')}
if (sigma.XM0<=0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (length(sigma.XM0)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 0 must be positive')}
if (sigma.XF0<=1) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XF1)==0) {stop('The standard deviation of the partner A Level 1 predictor in Group 1 must be positive')}
if (sigma.XM0<=1) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (length(sigma.XM1)==0) {stop('The standard deviation of the partner B Level 1 predictor in Group 1 must be positive')}
if (abs(rho.X0)>1) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (length(rho.X0)==0) {stop('The correlation between the parter F and M predictors in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.X1)>1) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(rho.X1)==0) {stop('The correlation between the parter F and M predictors in Group 1 must be included in the interval [-1,1]')}
if (length(c0) == 0) {stop('Set the value of the fixed intercept in Group 0')}
if (length(c1) == 0) {stop('Set the value of the fixed intercept in Group 1')}
if (length(a0) == 0) {stop('Set the value of the fixed actor linear effect in Group 0')}
if (length(p0) == 0) {stop('Set the value of the fixed partner linear effect in Group 0')}
if (length(a1) == 0) {stop('Set the value of the fixed actor linear effect in Group 1')}
if (length(p1) == 0) {stop('Set the value of the fixed partner linear effect in Group 1')}
if (length(a02) == 0) {stop('Set the value of the fixed actor quadratic effect in Group 0')}
if (length(p02) == 0) {stop('Set the value of the fixed partner quadratic effect in Group 0')}
if (length(a12) == 0) {stop('Set the value of the fixed actor quadratic effect in Group 1')}
if (length(p12) == 0) {stop('Set the value of the fixed partner quadratic effect in Group 1')}

if (Model == 12){
data = Sim.Dyad.Model.12(N0.dyad,N1.dyad,T.obs,  
c0,c1,a0,a1,a02,a12,p0,p1,p02,p12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}

if (Model == 28){
if (abs(rho.Y0)>1) {stop('The autoregressive effect in Group 0 must be included in the interval [-1,1]')}
if (length(rho.Y0)==0) {stop('The autoregressive effect in Group 0 must be included in the interval [-1,1]')}
if (abs(rho.Y0+rho.Y1)>1) {stop('The autoregressive effect in Group 1 must be included in the interval [-1,1]')}
if (length(rho.Y1)==0) {stop('The autoregressive effect in Group 1 must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.12.lag(N0.dyad,N1.dyad,T.obs,  
c0,c1,rho.Y0,rho,Y1,a0,a1,a02,a12,p0,p1,p02,p12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,
sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a continuos time-varying dyad moderator

if (Model == 13 | Model == 29){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(mu.W) == 0) {stop('Set the value of the mean of the Level 1 continuous dyad moderator')}
if (sigma.W == 0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(sigma.W)==0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(c.F) == 0) {stop('Set the value of the fixed intercept for partner A')}
if (length(c.M) == 0) {stop('Set the value of the fixed intercept for partner B')}
if (length(a.FF) == 0) {stop('Set the value of the fixed linear actor effect for partner A')}
if (length(p.MF) == 0) {stop('Set the value of the fixed linear partner effect for partner A')}
if (length(a.MM) == 0) {stop('Set the value of the fixed linear actor effect for partner B')}
if (length(p.FM) == 0) {stop('Set the value of the fixed linear partner effect for partner B')}
if (length(a.FF2) == 0) {stop('Set the value of the fixed quadratic actor effect for partner A')}
if (length(p.MF2) == 0) {stop('Set the value of the fixed quadratic partner effect for partner A')}
if (length(a.MM2) == 0) {stop('Set the value of the fixed quadratic actor effect for partner B')}
if (length(p.FM2) == 0) {stop('Set the value of the fixed quadratic partner effect for partner B')}
if (length(b.F) == 0) {stop('Set the value of the fixed effect of the Level 1 continuous dyad moderator for partner A')}
if (length(b.M) == 0) {stop('Set the value of the fixed effect of the Level 1 continuous dyad moderator for partner B')}
if (length(b.FF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor linear effect for partner A')}
if (length(b.MF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner linear effect for partner A')}
if (length(b.MM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor linear effect for partner B')}
if (length(b.FM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner linear effect for partner B')}
if (length(b.FF2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor quadratic effect for partner A')}
if (length(b.MF2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner quadratic effect for partner A')}
if (length(b.MM2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the actor quadratic effect for partner B')}
if (length(b.FM2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the partner quadratic effect for partner B')}

if (Model == 13){
data = Sim.Dyad.Model.13(N.dyad,T.obs,  
c.F,c.M,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
b.F,b.M,b.FF,b.FF2,b.MF,b.MF2,b.MM,b.MM2,b.FM,b.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}

if (Model == 29){
if (abs(rho.YF)>1) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (length(rho.YF)==0) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (abs(rho.YM)>1) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}
if (length(rho.YM)==0) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.13.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
b.F,b.M,b.FF,b.FF2,b.MF,b.MF2,b.MM,b.MM2,b.FM,b.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a continuos time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 14 | Model == 30){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(mu.W) == 0) {stop('Set the value of the mean of the Level 1 continuous dyad moderator')}
if (sigma.W == 0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(sigma.W)==0) {stop('The standard deviation of the Level 1 continuous dyad moderator')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed linear actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed linear partner effect')}
if (length(a.2) == 0) {stop('Set the value of the fixed quadratic actor effect')}
if (length(p.2) == 0) {stop('Set the value of the fixed quadratic partner effect')}
if (length(b) == 0) {stop('Set the value of the fixed effect of the Level 1 continuous dyad moderator')}
if (length(b.a) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the linear actor effect')}
if (length(b.p) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the linear partner effect')}
if (length(b.a2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the quadratic actor effect')}
if (length(b.p2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 continuous dyad moderator and the quadratic partner effect')}

if (Model == 14){
data = Sim.Dyad.Model.14(N.dyad,T.obs,  
c,a,a.2,p,p.2,b,b.a,b.a2,b.p,b.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}

if (Model == 30){
if (abs(rho.Y)>1) {stop('The autoregressive effect must be included in the interval [-1,1]')}
if (length(rho.Y)==0) {stop('The autoregressive effect must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.14.lag(N.dyad,T.obs,  
c,rho.Y,a,a.2,p,p.2,b,b.a,b.a2,b.p,b.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a dichotonomous time-varying dyad moderator

if (Model == 15 | Model == 31){
if (sigma.nu.F<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu.F)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (sigma.nu.M<=0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (length(sigma.nu.M)==0) {stop('The variance of the Level 2 partner B random intercept must be positive')}
if (abs(rho.eps.FM)>1) {stop('The correlation betwee the Level 2 partner A and M random intercepts errors must be included in the interval [-1,1]')}
if (length(rho.eps.FM)==0) {stop('The correlation betwee the Level 2 partner A and M random intercepts must be included in the interval [-1,1]')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (abs(prob.D)>1) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(prob.D)==0) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(c.F) == 0) {stop('Set the value of the fixed intercept for partner A')}
if (length(c.M) == 0) {stop('Set the value of the fixed intercept for partner B')}
if (length(a.FF) == 0) {stop('Set the value of the fixed linear actor effect for partner A')}
if (length(p.MF) == 0) {stop('Set the value of the fixed linear partner effect for partner A')}
if (length(a.MM) == 0) {stop('Set the value of the fixed linear actor effect for partner B')}
if (length(p.FM) == 0) {stop('Set the value of the fixed linear partner effect for partner B')}
if (length(a.FF2) == 0) {stop('Set the value of the fixed quadratic actor effect for partner A')}
if (length(p.MF2) == 0) {stop('Set the value of the fixed quadratic partner effect for partner A')}
if (length(a.MM2) == 0) {stop('Set the value of the fixed quadratic actor effect for partner B')}
if (length(p.FM2) == 0) {stop('Set the value of the fixed quadratic partner effect for partner B')}
if (length(d.F) == 0) {stop('Set the value of the fixed effect of the Level 1 dichotomous dyad moderator for partner A')}
if (length(d.M) == 0) {stop('Set the value of the fixed effect of the Level 1 dichotomous dyad moderator for partner B')}
if (length(d.FF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor linear effect for partner A')}
if (length(d.MF) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner linear effect for partner A')}
if (length(d.MM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor linear effect for partner B')}
if (length(d.FM) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner linear effect for partner B')}
if (length(d.FF2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor quadratic effect for partner A')}
if (length(d.MF2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner quadratic effect for partner A')}
if (length(d.MM2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the actor quadratic effect for partner B')}
if (length(d.FM2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the partner quadratic effect for partner B')}
if (sigma.W<= 0) {stop('The variance of the time-invarying predictor must be positive')}

if (Model == 15){
data = Sim.Dyad.Model.15(N.dyad,T.obs,  
c.F,c.M,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}

if (Model == 31){
if (abs(rho.YF)>1) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (length(rho.YF)==0) {stop('The autoregressive effect for partner A must be included in the interval [-1,1]')}
if (abs(rho.YM)>1) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}
if (length(rho.YM)==0) {stop('The autoregressive effect for partner B must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.15.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a dichotonomous time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 16 | Model == 32){
if (sigma.nu<=0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(sigma.nu)==0) {stop('The variance of the Level 2 partner A random intercept must be positive')}
if (length(mu.XF) == 0) {stop('Set the value of the mean of partner A of the Level 1 predictor')}
if (length(mu.XM) == 0) {stop('Set the value of the mean of partner B of the Level 1 predictor')}
if (sigma.XF<=0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (length(sigma.XF)==0) {stop('The standard deviation of the partner A Level 1 predictor must be positive')}
if (sigma.XM<=0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (length(sigma.XM)==0) {stop('The standard deviation of the partner B Level 1 predictor must be positive')}
if (abs(rho.X)>1) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (length(rho.X)==0) {stop('The correlation between the parter F and M predictors must be included in the interval [-1,1]')}
if (abs(prob.D)>1) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(prob.D)==0) {stop('The probability of the Level 1 dichotomous dyad moderator must be included in the interval [0,1]')}
if (length(c) == 0) {stop('Set the value of the fixed intercept')}
if (length(a) == 0) {stop('Set the value of the fixed linear actor effect')}
if (length(p) == 0) {stop('Set the value of the fixed linear partner effect')}
if (length(a.2) == 0) {stop('Set the value of the fixed quadratic actor effect')}
if (length(p.2) == 0) {stop('Set the value of the fixed quadratic partner effect')}
if (length(d) == 0) {stop('Set the value of the fixed effect of the Level 1 dichotomous dyad moderator')}
if (length(d.a) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the linear actor effect')}
if (length(d.p) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the linear partner effect')}
if (length(d.a2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the quadratic actor effect')}
if (length(d.p2) == 0) {stop('Set the value of the fixed interaction effect between the Level 1 dichotomous dyad moderator and the quadratic partner effect')}

if (Model == 16){
data = Sim.Dyad.Model.16(N.dyad,T.obs,  
c,a,a.2,p,p.2,d,d.a,d.a2,d.p,d.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}

if (Model == 32){
if (abs(rho.Y)>1) {stop('The autoregressive effect must be included in the interval [-1,1]')}
if (length(rho.Y)==0) {stop('The autoregressive effect must be included in the interval [-1,1]')}

data = Sim.Dyad.Model.16.lag(N.dyad,T.obs,  
c,rho.Y,a,a.2,p,p.2,d,d.a,d.a2,d.p,d.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X)
}}

# End of function ---> Return simulated IL dyad data set

return(data)}