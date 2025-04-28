###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for dyadic IL studies as a function of the number of participants

Summary.model.Dyad.IL = function(Model,N.dyad,N0.dyad,N1.dyad,T.obs,  
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
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,
mu.W,sigma.W,prob.D,
is.center.X,is.center.W,R,alpha,is.REML){

########################################################################################
########################################################################################
########################################################################################

# Simulate data from APIM model

if (Model == 1){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.1(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with indistinguishable dyads 

if (Model == 2){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.2(N.dyad,T.obs,  
c,a,p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model: group differences in actor partner effects

if (Model == 3){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.3(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}

########################################################################################

# Simulate data from APIM model: group differences in actor partner effects
# with indistinguishable dyads 

if (Model == 4){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.4(N0.dyad,N1.dyad,T.obs,  
c0,c1,a0,a1,p0,p1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}


########################################################################################

# Simulate data from APIM model with a continuos time-varying dyad moderator

if (Model == 5){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.5(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with a continuos time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 6){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.6(N.dyad,T.obs,  
c,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with a dichotonomous time-varying dyad moderator

if (Model == 7){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.7(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
d.F,d.M,d.FF,d.MF,d.MM,d.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with a dichotonomous time-varying dyad moderator
# with indistinguishable dyads

if (Model == 8){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.8(N.dyad,T.obs,  
c,a,p,d,d.a,d.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model 

if (Model == 9){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.9(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.FF2,p.MF2,a.MM,p.FM,a.MM2,p.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with indistinguishable dyads 

if (Model == 10){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.10(N.dyad,T.obs,  
c,a,a.2,p,p.2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model: group differences in actor partner effects

if (Model == 11){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.11(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,a.FF0,a.FF1,a.FF02,a.FF12,p.MF0,p.MF1,p.MF02,p.MF12,
a.MM0,a.MM1,a.MM02,a.MM12,p.FM0,p.FM1,p.FM02,p.FM12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}


########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model: group differences in actor partner effects
# with indistinguishable dyads 

if (Model == 12){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.12(N0.dyad,N1.dyad,T.obs,  
c0,c1,a0,a1,a02,a12,p0,p1,p02,p12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}


########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a continuos time-varying dyad moderator

if (Model == 13){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.13(N.dyad,T.obs,  
c.F,c.M,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
b.F,b.M,b.FF,b.FF2,b.MF,b.MF2,b.MM,b.MM2,b.FM,b.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a continuos time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 14){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.14(N.dyad,T.obs,  
c,a,a.2,p,p.2,b,b.a,b.a2,b.p,b.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a dichotonomous time-varying dyad moderator

if (Model == 15){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.15(N.dyad,T.obs,  
c.F,c.M,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a dichotonomous time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 16){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.16(N.dyad,T.obs,  
c,a,a.2,p,p.2,d,d.a,d.a2,d.p,d.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model

if (Model == 17){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.1.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with indistinguishable dyads 

if (Model == 18){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.2.lag(N.dyad,T.obs,  
c,rho.Y,a,p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model: group differences in actor partner effects

if (Model == 19){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.3.lag(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}

########################################################################################

# Simulate data from APIM model: group differences in actor partner effects
# with indistinguishable dyads 

if (Model == 20){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.4.lag(N0.dyad,N1.dyad,T.obs,  
c0,c1,rho.Y0,rho.Y1,a0,a1,p0,p1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}


########################################################################################

# Simulate data from APIM model with a continuos time-varying dyad moderator

if (Model == 21){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.5.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with a continuos time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 22){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.6.lag(N.dyad,T.obs,  
c,rho.Y,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with a dichotonomous time-varying dyad moderator

if (Model == 23){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.7.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.MM,p.FM,
d.F,d.M,d.FF,d.MF,d.MM,d.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Simulate data from APIM model with a dichotonomous time-varying dyad moderator
# with indistinguishable dyads

if (Model == 24){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.8.lag(N.dyad,T.obs,  
c,rho.Y,a,p,d,d.a,d.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model 

if (Model == 25){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.9.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,p.MF,a.FF2,p.MF2,a.MM,p.FM,a.MM2,p.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with indistinguishable dyads 

if (Model == 26){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.10.lag(N.dyad,T.obs,  
c,rho.Y,a,a.2,p,p.2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model: group differences in actor partner effects

if (Model == 27){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.11.lag(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,a.FF0,a.FF1,a.FF02,a.FF12,p.MF0,p.MF1,p.MF02,p.MF12,
a.MM0,a.MM1,a.MM02,a.MM12,p.FM0,p.FM1,p.FM02,p.FM12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}


########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model: group differences in actor partner effects
# with indistinguishable dyads 

if (Model == 28){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.12.lag(N0.dyad,N1.dyad,T.obs,  
c0,c1,rho.Y0,rho.Y1,a0,a1,a02,a12,p0,p1,p02,p12,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads in Group 0 larger to',
N0.dyad,'or the number of dyads in Group 1 larger to',N1.dyad))}}
}


########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a continuos time-varying dyad moderator

if (Model == 29){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.13.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
b.F,b.M,b.FF,b.FF2,b.MF,b.MF2,b.MM,b.MM2,b.FM,b.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a continuos time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 30){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.14.lag(N.dyad,T.obs,  
c,rho.Y,a,a.2,p,p.2,b,b.a,b.a2,b.p,b.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a dichotonomous time-varying dyad moderator

if (Model == 31){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.15.lag(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

########################################################################################

# Curvilinear actor and partner effects
# Simulate data from APIM model with a dichotonomous time-varying dyad moderator
# with indistinguishable dyads 

if (Model == 32){
plan(multisession)
fit = future_lapply(1:R, function(r)
try(Performance.Dyad.Model.16.lag(N.dyad,T.obs,  
c,rho.Y,a,a.2,p,p.2,d,d.a,d.a2,d.p,d.p2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X,alpha,is.REML),
silent = FALSE),future.seed = 0xBEEF)

errors = rep(0,R)
for (r in 1:R){errors[r] = length(fit[[r]])}
errors = sum(ifelse(errors==1,1,0)) 

if (errors>0){
if (errors<(R-1)){message(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}
if (errors>=(R-1)){stop(paste(errors, 'replications produce convergence errors. 
Check the value of the parameters or set the number of dyads larger to',
N.dyad))}}
}

# End of function ---> Return estimated model

return(fit)}

#####################################################################################