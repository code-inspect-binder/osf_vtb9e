###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for IL studies as a function of the number of dyads

Sim.model.Dyad.IL = function(Model,N.dyad,N0.dyad,N1.dyad,T.obs,  
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
  
  ###############################################################
  ###############################################################
  ###############################################################
  
  message('Initializing...')
  
  # Check parameters & and output messages
  
  if (length(alpha) == 0) {stop('The Type I error must be between 0 and 1')}
  if (alpha > 1) {stop('The Type I error must be between 0 and 1')}
  if (alpha == 0) {stop('The Type I error must be between 0 and 1')}
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  # Simulate data from the L-APIM model
  
  if (Model == 1 | Model == 2 | Model == 5 | Model == 6 | Model == 7 | Model == 8 | 
      Model == 9 | Model == 10 | Model == 13 | Model == 14 | Model == 15 | Model == 16 | 
      Model == 17 | Model == 18 | Model == 21 | Model == 22 | Model == 23 | Model == 24 | 
      Model == 25 | Model == 26 | Model == 29 | Model == 30 | Model == 31 | Model == 32){
    
    N.dyad = as.numeric(unlist(strsplit(N.dyad,",")))
    
    # Monte Carlo replicates
    
    plan(multisession)
    fit = future_lapply(1:length(N.dyad), function(i){
      message(paste('Estimating Power for Number of dyads =',N.dyad[i]))
      Summary.model.Dyad.IL(Model,N.dyad[i],N0.dyad,N1.dyad,T.obs,  
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
                            is.center.X,is.center.W,R,alpha,is.REML)},future.seed = 0xBEEF)
    
    # Distribution of the parameters: fixed effects 
    
    Rr = lapply(1:length(N.dyad), function(i) 
      which(lapply(1:R, function(r) length(fit[[i]][[r]])==3)==TRUE))
    
    Coef.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$coef)),nrow=length(Rr[[i]]),byrow=T))
    Coef.bias.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$bias)),nrow=length(Rr[[i]]),byrow=T))
    Coef.power.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$power)),nrow=length(Rr[[i]]),byrow=T))
    Coef.CI.width.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$CI.width)),nrow=length(Rr[[i]]),byrow=T))
    Coef.CR.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$CR)),nrow=length(Rr[[i]]),byrow=T))
    
    Coef.hat = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Coef.list[[i]]))),nrow=length(N.dyad),byrow=T)
    Coef.se = matrix(unlist(lapply(1:length(N.dyad), function(i) apply(Coef.list[[i]],2,sd)/sqrt(length(Rr[[i]])))),nrow=length(N.dyad),byrow=T)
    Coef.bias = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Coef.bias.list[[i]]))),nrow=length(N.dyad),byrow=T)
    Coef.CI.Width = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Coef.CI.width.list[[i]]))),nrow=length(N.dyad),byrow=T)
    Coef.CR = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Coef.CR.list[[i]]))),nrow=length(N.dyad),byrow=T)
    Coef.power = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Coef.power.list[[i]]))),nrow=length(N.dyad),byrow=T)
    Coef.power.se = matrix(unlist(lapply(1:length(N.dyad), function(i) apply(Coef.power.list[[i]],2,sd)/sqrt(length(Rr[[i]])))),nrow=length(N.dyad),byrow=T)
    
    # Distribution of the parameters: random effects
    
    Var.hat.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.var.hat)),nrow=length(Rr[[i]]),byrow=T))
    Var.bias.list = lapply(1:length(N.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.var.bias)),nrow=length(Rr[[i]]),byrow=T))
    
    Var.hat = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Var.hat.list[[i]]))),nrow=length(N.dyad),byrow=T)
    Var.hat.se = matrix(unlist(lapply(1:length(N.dyad), function(i) apply(Var.hat.list[[i]],2,sd)/sqrt(length(Rr[[i]])))),nrow=length(N.dyad),byrow=T)
    Var.bias = matrix(unlist(lapply(1:length(N.dyad), function(i) colMeans(Var.bias.list[[i]]))),nrow=length(N.dyad),byrow=T)
    
    # Number of replicates
    
    replicates = unlist(lapply(1:length(N.dyad), function(i) paste('(',length(Rr[[i]]),', Number of Dyad=',N.dyad[i],')')))
  }
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  if (Model == 3 | Model == 4 | Model == 11 | Model == 12 | 
      Model == 19 | Model == 20 | Model == 27 | Model == 28){
    
    N0.dyad = as.numeric(unlist(strsplit(N0.dyad,",")))
    N1.dyad = as.numeric(unlist(strsplit(N1.dyad,",")))
    
    # Monte Carlo replicates
    
    plan(multisession)
    fit = future_lapply(1:length(N0.dyad), function(i){
      message(paste('Estimating power for number of dyads in Group 0 =',N0.dyad[i],'and number of dyads in Group 1=',N0.dyad[i]))
      Summary.model.Dyad.IL(Model,N.dyad,N0.dyad[i],N1.dyad[i],T.obs,  
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
                            is.center.X,is.center.W,R,alpha,is.REML)},future.seed = 0xBEEF)
    
    # Distribution of the parameters: fixed effects 
    
    Rr = lapply(1:length(N0.dyad), function(i) 
      which(lapply(1:R, function(r) length(fit[[i]][[r]])==3)==TRUE))
    
    Coef.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$coef)),nrow=length(Rr[[i]]),byrow=T))
    Coef.bias.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$bias)),nrow=length(Rr[[i]]),byrow=T))
    Coef.power.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$power)),nrow=length(Rr[[i]]),byrow=T))
    Coef.CI.width.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$CI.width)),nrow=length(Rr[[i]]),byrow=T))
    Coef.CR.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.fixed.effect$CR)),nrow=length(Rr[[i]]),byrow=T))
    
    Coef.hat = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Coef.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    Coef.se = matrix(unlist(lapply(1:length(N0.dyad), function(i) apply(Coef.list[[i]],2,sd)/sqrt(length(Rr[[i]])))),nrow=length(N0.dyad),byrow=T)
    Coef.bias = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Coef.bias.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    Coef.CI.Width = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Coef.CI.width.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    Coef.CR = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Coef.CR.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    Coef.power = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Coef.power.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    Coef.power.se = matrix(unlist(lapply(1:length(N0.dyad), function(i) apply(Coef.power.list[[i]],2,sd)/sqrt(length(Rr[[i]])))),nrow=length(N0.dyad),byrow=T)
    
    # Distribution of the parameters: random effects
    
    Var.hat.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.var.hat)),nrow=length(Rr[[i]]),byrow=T))
    Var.bias.list = lapply(1:length(N0.dyad), function(i)
      matrix(unlist(lapply(Rr[[i]], function(r) fit[[i]][[r]]$summary.var.bias)),nrow=length(Rr[[i]]),byrow=T))
    
    Var.hat = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Var.hat.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    Var.hat.se = matrix(unlist(lapply(1:length(N0.dyad), function(i) apply(Var.hat.list[[i]],2,sd)/sqrt(length(Rr[[i]])))),nrow=length(N0.dyad),byrow=T)
    Var.bias = matrix(unlist(lapply(1:length(N0.dyad), function(i) colMeans(Var.bias.list[[i]]))),nrow=length(N0.dyad),byrow=T)
    
    # Number of replicates
    
    replicates = unlist(lapply(1:length(N0.dyad), function(i) paste('(',length(Rr[[i]]),', Number of Dyad=',N0.dyad[i],N1.dyad[i],')')))
  }
  
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  
  # Simulate data from the L-APIM model 
  
  if (Model == 1){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),rep(a.FF,length(N.dyad)),rep(p.MF,length(N.dyad)),rep(a.MM,length(N.dyad)),rep(p.FM,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],Coef.hat[,6])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],Coef.bias[,6])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],Coef.se[,6])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],Coef.CI.Width[,6])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],Coef.CR[,6])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],Coef.power[,6])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],Coef.power.se[,5],Coef.power.se[,6])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','a.FF','p.MF','a.MM','p.FM',
                              'c.F.se','c.M.se','a.FF.se','p.MF.se','a.MM.se','p.FM.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from the L-APIM model with indistinguishable dyads
  
  if (Model == 2){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(a,length(N.dyad)),rep(p,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','a','p',
                              'c.se','a.se','p.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from L-APIM model: group differences in actor partner effects
  
  if (Model == 3){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F0,length(N0.dyad)),rep(c.F1,length(N0.dyad)),
                           rep(c.M,length(N0.dyad)),rep(c.M0,length(N0.dyad)),
                           rep(a.FF0,length(N0.dyad)),rep(a.FF1,length(N0.dyad)),
                           rep(p.MF0,length(N0.dyad)),rep(p.MF1,length(N0.dyad)),
                           rep(a.MM0,length(N0.dyad)),rep(a.MM1,length(N0.dyad)),
                           rep(p.FM0,length(N0.dyad)),rep(p.FM1,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Actor effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in actor effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Partner effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in partner effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Actor effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in actor effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Partner effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in partner effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu.F,length(N0.dyad)),
                     rep(sigma.nu.M,length(N0.dyad)),rep(rho.nu.F.M,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between random intercepts for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F1','c.M0','c.M1','a.FF0',
                              'a.FF1','p.MF0','p.MF1','a.MM0','a.MM1','p.FM0','p.FM1',
                              'c.F0.se','c.F1.se','c.M0.se','c.M1.se','a.FF0.se',
                              'a.FF1.se','p.MF0.se','p.MF1.se','a.MM0.se','a.MM1.se','p.FM0.se','p.FM1.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from L-APIM model: group differences in actor partner effects
  # with indistinguishable dyads 
  
  if (Model == 4){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c0,length(N0.dyad)),rep(c1,length(N0.dyad)),
                           rep(a0,length(N0.dyad)),rep(a1,length(N0.dyad)),
                           rep(p0,length(N0.dyad)),rep(p1,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Actor effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in actor effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Partner effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in partner effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c1','a0',
                              'a1','p0','p1',
                              'c0.se','c1.se','a0.se',
                              'a1.se','p0.se','p1.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from the L-APIM model with a continuous time-varying dyad moderator
  
  if (Model == 5){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(p.MF,length(N.dyad)),
                           rep(b.F,length(N.dyad)),rep(b.FF,length(N.dyad)),rep(b.MF,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(p.FM,length(N.dyad)),
                           rep(b.M,length(N.dyad)),rep(b.MM,length(N.dyad)),rep(b.FM,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the Partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','a.FF','p.MF','b.F','b.FF','b.MF','a.MM','p.FM','b.M','b.MM','b.FM',
                              'c.F.se','c.M.se','a.FF.se','p.MF.se','b.F.se','b.FF.se','b.MF.se','a.MM.se','p.FM.se','b.M.se','b.MM.se','b.FM.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from the L-APIM model with a continuous time-varying dyad moderator
  # with indistinguishable dyads 
  
  if (Model == 6){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(p,length(N.dyad)),
                           rep(b,length(N.dyad)),rep(b.a,length(N.dyad)),rep(b.p,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','a','p','b','b.a','b.p',
                              'c.se','a.se','p.se','b.se','b.a.se','b.p.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from the L-APIM model with a dichotomous time-varying dyad moderator
  
  if (Model == 7){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(p.MF,length(N.dyad)),
                           rep(d.F,length(N.dyad)),rep(d.FF,length(N.dyad)),rep(d.MF,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(p.FM,length(N.dyad)),
                           rep(d.M,length(N.dyad)),rep(d.MM,length(N.dyad)),rep(d.FM,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','a.FF','p.MF','d.F','d.FF','d.MF','a.MM','p.FM','d.M','d.MM','d.FM',
                              'c.F.se','c.M.se','a.FF.se','p.MF.se','d.F.se','d.FF.se','d.MF.se','a.MM.se','p.FM.se','d.M.se','d.MM.se','d.FM.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from the L-APIM model with a dichotomous time-varying dyad moderator
  # with indistinguishable dyads
  
  if (Model == 8){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(p,length(N.dyad)),
                           rep(d,length(N.dyad)),rep(d.a,length(N.dyad)),rep(d.p,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the Partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','a','p','d','d.a','d.p',
                              'c.se','a.se','p.se','d.se','d.a.se','d.p.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L-APIM model 
  
  if (Model == 9){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(a.FF2,length(N.dyad)),
                           rep(p.MF,length(N.dyad)),rep(p.MF2,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(a.MM2,length(N.dyad)),
                           rep(p.FM,length(N.dyad)),rep(p.FM2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],Coef.power.se[,5],
                               Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],Coef.power.se[,9],Coef.power.se[,10])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','a.FF','a.FF2','p.MF','p.MF2','a.MM','a.MM2','p.FM','p.FM2',
                              'c.F.se','c.M.se','a.FF.se','a.FF2.se','p.MF.se','p.MF2.se','a.MM.se','a.MM2.se','p.FM.se','p.FM2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model with indistinguishable dyads 
  
  if (Model == 10){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(a.2,length(N.dyad)),
                           rep(p,length(N.dyad)),rep(p.2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],Coef.power.se[,5])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linar partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','a','a.2','p','p.2',
                              'c.se','a.se','a.2.se','p.se','p.2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from L-APIM model: group differences in actor partner effects
  
  if (Model == 11){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F0,length(N0.dyad)),rep(c.F1,length(N0.dyad)),
                           rep(c.M,length(N0.dyad)),rep(c.M0,length(N0.dyad)),
                           rep(a.FF0,length(N0.dyad)),rep(a.FF02,length(N0.dyad)),
                           rep(a.FF1,length(N0.dyad)),rep(a.FF12,length(N0.dyad)),
                           rep(p.MF0,length(N0.dyad)),rep(p.MF02,length(N0.dyad)),
                           rep(p.MF1,length(N0.dyad)),rep(p.MF12,length(N0.dyad)),
                           rep(a.MM0,length(N0.dyad)),rep(a.MM02,length(N0.dyad)),
                           rep(a.MM1,length(N0.dyad)),rep(a.MM12,length(N0.dyad)),
                           rep(p.FM0,length(N0.dyad)),rep(p.FM02,length(N0.dyad)),
                           rep(p.FM1,length(N0.dyad)),rep(p.FM12,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16],
                          Coef.hat[,17],Coef.hat[,18],Coef.hat[,19],Coef.hat[,20])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16],
                           Coef.bias[,17],Coef.bias[,18],Coef.bias[,19],Coef.bias[,20]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16],
                         Coef.se[,17],Coef.se[,18],Coef.se[,19],Coef.se[,20])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16],
                               Coef.CI.Width[,17],Coef.CI.Width[,18],Coef.CI.Width[,19],Coef.CI.Width[,20]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16],
                         Coef.CR[,17],Coef.CR[,18],Coef.CR[,19],Coef.CR[,20])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16],
                            Coef.power[,17],Coef.power[,18],Coef.power[,19],Coef.power[,20])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16],
                               Coef.power.se[,17],Coef.power.se[,18],Coef.power.se[,19],Coef.power.se[,20])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner B in Group 0','N Dyad','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear actor effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic actor effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear actor effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic actor effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear partner effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic partner effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear partner effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic partner effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear actor effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic actor effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear actor effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic actor effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear partner effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic partner effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear partner effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic partner effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu.F,length(N0.dyad)),
                     rep(sigma.nu.M,length(N0.dyad)),rep(rho.nu.F.M,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between random intercepts for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F1','c.M0','c.M1','a.FF0','a.FF02','a.FF1',
                              'a.FF12','p.MF0','p.MF02','p.MF1','p.MF12','a.MM0','a.MM02','a.MM1','a.MM12',
                              'p.FM0','p.FM02','p.FM1','p.FM12',
                              'c.F0.se','c.F1.se','c.M0.se','c.M1.se','a.FF0.se','a.FF02.se','a.FF1.se',
                              'a.FF12.se','p.MF0.se','p.MF02.se','p.MF1.se','p.MF12.se','a.MM0.se','a.MM02.se','a.MM1.se','a.MM12.se',
                              'p.FM0.se','p.FM02.se','p.FM1.se','p.FM12.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model: group differences in actor partner effects
  # with indistinguishable dyads 
  
  if (Model == 12){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c0,length(N.dyad)),rep(c1,length(N0.dyad)),
                           rep(a0,length(N0.dyad)),rep(a02,length(N0.dyad)),
                           rep(a1,length(N0.dyad)),rep(a12,length(N0.dyad)),
                           rep(p0,length(N0.dyad)),rep(p02,length(N0.dyad)),
                           rep(p1,length(N0.dyad)),rep(p12,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear actor effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic actor effect in Group 0','N Dyad','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear actor effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic actor effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear partner effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic partner effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear partner effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic partner effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c1','a0','a02','a1','a12','p0','p02','p1','p12',
                              'c0.se','c1.se','a0.se','a02.se','a1.se','a12.se','p0.se','p02.se','p1.se','p12.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L- APIM model with a continuous time-varying dyad moderator
  
  if (Model == 13){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(a.FF2,length(N.dyad)),
                           rep(p.MF,length(N.dyad)),rep(p.MF2,length(N.dyad)),
                           rep(b.F,length(N.dyad)),
                           rep(b.FF,length(N.dyad)),rep(b.FF2,length(N.dyad)),
                           rep(b.MF,length(N.dyad)),rep(b.MF2,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(a.MM2,length(N.dyad)),
                           rep(p.FM,length(N.dyad)),rep(p.FM2,length(N.dyad)),
                           rep(b.M,length(N.dyad)),
                           rep(b.MM,length(N.dyad)),rep(b.MM2,length(N.dyad)),
                           rep(b.FM,length(N.dyad)),rep(b.FM2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16],
                          Coef.hat[,17],Coef.hat[,18],Coef.hat[,19],Coef.hat[,20])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16],
                           Coef.bias[,17],Coef.bias[,18],Coef.bias[,19],Coef.bias[,20])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16],
                         Coef.se[,17],Coef.se[,18],Coef.se[,19],Coef.se[,20])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16],
                               Coef.CI.Width[,17],Coef.CI.Width[,18],Coef.CI.Width[,19],Coef.CI.Width[,20])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16],
                         Coef.CR[,17],Coef.CR[,18],Coef.CR[,19],Coef.CR[,20])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16],
                            Coef.power[,17],Coef.power[,18],Coef.power[,19],Coef.power[,20])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16],
                               Coef.power.se[,17],Coef.power.se[,18],Coef.power.se[,19],Coef.power.se[,20])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear aActor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','a.FF','a.FF2','p.MF','p.MF2','b.F','b.FF','b.FF2','b.MF','b.MF2',
                              'a.MM','a.MM2','p.FM','p.FM2','b.M','b.MM','b.MM2','b.FM','b.FM2',
                              'c.F.se','c.M.se','a.FF.se','a.FF2.se','p.MF.se','p.MF2.se','b.F.se','b.FF.se','b.FF2.se','b.MF.se','b.MF2.se',
                              'a.MM.se','a.MM2.se','p.FM.se','p.FM2.se','b.M.se','b.MM.se','b.MM2.se','b.FM.se','b.FM2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L-APIM model with a continuous time-varying dyad moderator
  # with indistinguishable dyads 
  
  if (Model == 14){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(a.2,length(N.dyad)),
                           rep(p,length(N.dyad)),rep(p.2,length(N.dyad)),
                           rep(b,length(N.dyad)),
                           rep(b.a,length(N.dyad)),rep(b.a2,length(N.dyad)),
                           rep(b.p,length(N.dyad)),rep(b.p2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','a','a.2','p','p.2','b','b.a','b.a2','b.p','b.p2',
                              'c.se','a.se','a.2.se','p.se','p.2.se','b.se','b.a.se','b.a2.se','b.p.se','b.p2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L-APIM model with a dichotomous time-varying dyad moderator
  
  if (Model == 15){
    
    # Table fixed predictors  
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(a.FF2,length(N.dyad)),
                           rep(p.MF,length(N.dyad)),rep(p.MF2,length(N.dyad)),
                           rep(d.F,length(N.dyad)),
                           rep(d.FF,length(N.dyad)),rep(d.FF2,length(N.dyad)),
                           rep(d.MF,length(N.dyad)),rep(d.MF2,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(a.MM2,length(N.dyad)),
                           rep(p.FM,length(N.dyad)),rep(p.FM2,length(N.dyad)),
                           rep(d.M,length(N.dyad)),
                           rep(d.MM,length(N.dyad)),rep(d.MM2,length(N.dyad)),
                           rep(d.FM,length(N.dyad)),rep(d.FM2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16],
                          Coef.hat[,17],Coef.hat[,18],Coef.hat[,19],Coef.hat[,20])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16],
                           Coef.bias[,17],Coef.bias[,18],Coef.bias[,19],Coef.bias[,20])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16],
                         Coef.se[,17],Coef.se[,18],Coef.se[,19],Coef.se[,20])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16],
                               Coef.CI.Width[,17],Coef.CI.Width[,18],Coef.CI.Width[,19],Coef.CI.Width[,20])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16],
                         Coef.CR[,17],Coef.CR[,18],Coef.CR[,19],Coef.CR[,20])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16],
                            Coef.power[,17],Coef.power[,18],Coef.power[,19],Coef.power[,20])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16],
                               Coef.power.se[,17],Coef.power.se[,18],Coef.power.se[,19],Coef.power.se[,20])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','a.FF','a.FF2','p.MF','p.MF2','d.F','d.FF','d.FF2','d.MF','d.MF2',
                              'a.MM','a.MM2','p.FM','p.FM2','d.M','d.MM','d.MM2','d.FM','d.FM2',
                              'c.F.se','c.M.se','a.FF.se','a.FF2.se','p.MF.se','p.MF2.se','d.F.se','d.FF.se','d.FF2.se','d.MF.se','d.MF2.se',
                              'a.MM.se','a.MM2.se','p.FM.se','p.FM2.se','d.M.se','d.MM.se','d.MM2.se','d.FM.se','d.FM2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L-APIM model with a dichotomous time-varying dyad moderator
  # with indistinguishable dyads 
  
  if (Model == 16){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(a.2,length(N.dyad)),
                           rep(p,length(N.dyad)),rep(p.2,length(N.dyad)),
                           rep(d,length(N.dyad)),
                           rep(d.a,length(N.dyad)),rep(d.a2,length(N.dyad)),
                           rep(d.p,length(N.dyad)),rep(d.p2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','a','a.2','p','p.2','d','d.a','d.a2','d.p','d.p2',
                              'c.se','a.se','a.2.se','p.se','p.2.se','d.se','d.a.se','d.a2.se','d.p.se','d.p2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  #####################################################################################
  
  # Simulate data from the L-APIM model - AR(1)
  
  if (Model == 17){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(rho.YF,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(p.MF,length(N.dyad)),
                           rep(rho.YM,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(p.FM,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','rho.YF','a.FF','p.MF','rho.YM','a.MM','p.FM',
                              'c.F.se','c.M.se','rho.YF.se','a.FF.se','p.MF.se','rho.YM.se','a.MM.se','p.FM.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model with indistinguishable dyads  - AR(1)
  
  if (Model == 18){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(rho.Y,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(p,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','rho.Y','a','p',
                              'c.se','rho.Y.se','a.se','p.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model: group differences in actor partner effects - AR(1)
  
  if (Model == 19){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F0,length(N0.dyad)),rep(c.F1,length(N0.dyad)),
                           rep(c.M,length(N0.dyad)),rep(c.M0,length(N0.dyad)),
                           rep(rho.YF0,length(N0.dyad)),rep(rho.YF1,length(N0.dyad)),
                           rep(a.FF0,length(N0.dyad)),rep(a.FF1,length(N0.dyad)),
                           rep(p.MF0,length(N0.dyad)),rep(p.MF1,length(N0.dyad)),
                           rep(rho.YM0,length(N0.dyad)),rep(rho.YM1,length(N0.dyad)),
                           rep(a.MM0,length(N0.dyad)),rep(a.MM1,length(N0.dyad)),
                           rep(p.FM0,length(N0.dyad)),rep(p.FM1,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Autoregressive effect for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the Autoregressive effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Actor effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in actor effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Partner effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in partner effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Autoregressive effect for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the Autoregressive effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Actor effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in actor effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Partner effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in partner effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu.F,length(N0.dyad)),
                     rep(sigma.nu.M,length(N0.dyad)),rep(rho.nu.F.M,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between random intercepts for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F1','c.M0','c.M1',
                              'rho.YF0','rho.YF1','a.FF0','a.FF1','p.MF0','p.MF1','rho.YM0','rho.YM1','a.MM0','a.MM1','p.FM0','p.FM1',
                              'c.F0.se','c.F1.se','c.M0.se','c.M1.se',
                              'rho.YF0.se','rho.YF1.se','a.FF0.se','a.FF1.se','p.MF0.se','p.MF1.se','rho.YM0.se','rho.YM1.se','a.MM0.se','a.MM1.se','p.FM0.se','p.FM1.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model: group differences in actor partner effects - AR(1)
  # with indistinguishable dyads 
  
  if (Model == 20){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c0,length(N0.dyad)),rep(c1,length(N0.dyad)),
                           rep(rho.Y0,length(N0.dyad)),rep(rho.Y1,length(N0.dyad)),
                           rep(a0,length(N0.dyad)),rep(a1,length(N0.dyad)),
                           rep(p0,length(N0.dyad)),rep(p1,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Autoregressive effect for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the autoregressive effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Actor effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in actor effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Partner effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in partner effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c1','rho.Y0','rho.Y1',
                              'a0','a1','p0','p1',
                              'c0.se','c1.se','rho.Y0.se','rho.Y1.se',
                              'a0.se','a1.se','p0.se','p1.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model with a continuos time-varying dyad moderator - AR(1)
  
  if (Model == 21){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(rho.YF,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(p.MF,length(N.dyad)),
                           rep(b.F,length(N.dyad)),rep(b.FF,length(N.dyad)),rep(b.MF,length(N.dyad)),
                           rep(rho.YM,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(p.FM,length(N.dyad)),
                           rep(b.M,length(N.dyad)),rep(b.MM,length(N.dyad)),rep(b.FM,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner P','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the Partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','rho.YF','a.FF','p.MF','b.F','b.FF','b.MF','rho.YM','a.MM','p.FM','b.M','b.MM','b.FM',
                              'c.F.se','c.M.se','rho.YF.se','a.FF.se','p.MF.se','b.F.se','b.FF.se','b.MF.se','rho.YM.se','a.MM.se','p.FM.se','b.M.se','b.MM.se','b.FM.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model with a continuos time-varying dyad moderator - AR(1)
  # with indistinguishable dyads 
  
  if (Model == 22){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(rho.Y,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(p,length(N.dyad)),
                           rep(b,length(N.dyad)),rep(b.a,length(N.dyad)),rep(b.p,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','rho.Y','a','p','b','b.a','b.p',
                              'c.se','rho.Y.se','a.se','p.se','b.se','b.a.se','b.p.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model with a dichotomous time-varying dyad moderator - AR(1)
  
  if (Model == 23){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(rho.YF,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(p.MF,length(N.dyad)),
                           rep(d.F,length(N.dyad)),rep(d.FF,length(N.dyad)),rep(d.MF,length(N.dyad)),
                           rep(rho.YM,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(p.FM,length(N.dyad)),
                           rep(d.M,length(N.dyad)),rep(d.MM,length(N.dyad)),rep(d.FM,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','rho.YF','a.FF','p.MF','d.F','d.FF','d.MF','rho.YM','a.MM','p.FM','d.M','d.MM','d.FM',
                              'c.F.se','c.M.se','rho.YF.se','a.FF.se','p.MF.se','d.F.se','d.FF.se','d.MF.se','rho.YM.se','a.MM.se','p.FM.se','d.M.se','d.MM.se','d.FM.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Simulate data from APIM model with a dichotomous time-varying dyad moderator - AR(1)
  # with indistinguishable dyads
  
  if (Model == 24){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(rho.Y,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(p,length(N.dyad)),
                           rep(d,length(N.dyad)),rep(d.a,length(N.dyad)),rep(d.p,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the Partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','rho.Y','a','p','d','d.a','d.p',
                              'c.se','rho.Y.se','a.se','p.se','d.se','d.a.se','d.p.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model - AR(1) 
  
  if (Model == 25){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(rho.YF,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(a.FF2,length(N.dyad)),
                           rep(p.MF,length(N.dyad)),rep(p.MF2,length(N.dyad)),
                           rep(rho.YM,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(a.MM2,length(N.dyad)),
                           rep(p.FM,length(N.dyad)),rep(p.FM2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','rho.YF','a.FF','a.FF2','p.MF','p.MF2','rho.YM','a.MM','a.MM2','p.FM','p.FM2',
                              'c.F.se','c.M.se','rho.YF.se','a.FF.se','a.FF2.se','p.MF.se','p.MF2.se','rho.YM.se','a.MM.se','a.MM2.se','p.FM.se','p.FM2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model with indistinguishable dyads - AR(1) 
  
  if (Model == 26){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(rho.Y,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(a.2,length(N.dyad)),
                           rep(p,length(N.dyad)),rep(p.2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],Coef.hat[,6])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],Coef.bias[,6])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],Coef.se[,6])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],Coef.CI.Width[,6])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],Coef.CR[,6])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],Coef.power[,6])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linar partner Effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner Effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','rho.Y','a','a.2','p','p.2',
                              'c.se','rho.Y.se','a.se','a.2.se','p.se','p.2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L-APIM model: group differences in actor partner effects - AR(1)
  
  if (Model == 27){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F0,length(N0.dyad)),rep(c.F1,length(N0.dyad)),
                           rep(c.M,length(N0.dyad)),rep(c.M0,length(N0.dyad)),
                           rep(rho.YF0,length(N0.dyad)),rep(rho.YF1,length(N0.dyad)),
                           rep(a.FF0,length(N0.dyad)),rep(a.FF02,length(N0.dyad)),
                           rep(a.FF1,length(N0.dyad)),rep(a.FF12,length(N0.dyad)),
                           rep(p.MF0,length(N0.dyad)),rep(p.MF02,length(N0.dyad)),
                           rep(p.MF1,length(N0.dyad)),rep(p.MF12,length(N0.dyad)),
                           rep(rho.YM0,length(N0.dyad)),rep(rho.YM1,length(N0.dyad)),
                           rep(a.MM0,length(N0.dyad)),rep(a.MM02,length(N0.dyad)),
                           rep(a.MM1,length(N0.dyad)),rep(a.MM12,length(N0.dyad)),
                           rep(p.FM0,length(N0.dyad)),rep(p.FM02,length(N0.dyad)),
                           rep(p.FM1,length(N0.dyad)),rep(p.FM12,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16],
                          Coef.hat[,17],Coef.hat[,18],Coef.hat[,19],Coef.hat[,20],
                          Coef.hat[,21],Coef.hat[,22],Coef.hat[,23],Coef.hat[,24])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16],
                           Coef.bias[,17],Coef.bias[,18],Coef.bias[,19],Coef.bias[,20],
                           Coef.bias[,21],Coef.bias[,22],Coef.bias[,23],Coef.bias[,24]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16],
                         Coef.se[,17],Coef.se[,18],Coef.se[,19],Coef.se[,20],
                         Coef.se[,21],Coef.se[,22],Coef.se[,23],Coef.se[,24])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16],
                               Coef.CI.Width[,17],Coef.CI.Width[,18],Coef.CI.Width[,19],Coef.CI.Width[,20],
                               Coef.CI.Width[,21],Coef.CI.Width[,22],Coef.CI.Width[,23],Coef.CI.Width[,24]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16],
                         Coef.CR[,17],Coef.CR[,18],Coef.CR[,19],Coef.CR[,20],
                         Coef.CR[,21],Coef.CR[,22],Coef.CR[,23],Coef.CR[,24])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16],
                            Coef.power[,17],Coef.power[,18],Coef.power[,19],Coef.power[,20],
                            Coef.power[,21],Coef.power[,22],Coef.power[,23],Coef.power[,24])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16],
                               Coef.power.se[,17],Coef.power.se[,18],Coef.power.se[,19],Coef.power.se[,20],
                               Coef.power.se[,21],Coef.power.se[,22],Coef.power.se[,23],Coef.power.se[,24])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Autoregressive effect for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the autoregressive effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear actor effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic actor effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear actor effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic actor effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear partner effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic partner effect for partner A in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear partner effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic partner effect between Group 0 and 1 for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Autoregressive effect for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the autoregressive effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear actor effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic actor effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear actor effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic actor effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear partner effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic partner effect for partner B in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear partner effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic partner effect between Group 0 and 1 for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu.F,length(N0.dyad)),
                     rep(sigma.nu.M,length(N0.dyad)),rep(rho.nu.F.M,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between random intercepts for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F1','rho.YF0','rho.YF1','c.M0','c.M1','a.FF0','a.FF02','a.FF1',
                              'a.FF12','p.MF0','p.MF02','p.MF1','p.MF12','rho.YM0','rho.YM1','a.MM0','a.MM02','a.MM1','a.MM12',
                              'p.FM0','p.FM02','p.FM1','p.FM12',
                              'c.F0.se','c.F1.se','rho.YF0.se','rho.YF1.se','c.M0.se','c.M1.se','a.FF0.se','a.FF02.se','a.FF1.se',
                              'a.FF12.se','p.MF0.se','p.MF02.se','p.MF1.se','p.MF12.se','rho.YM0.se','rho.YM1.se','a.MM0.se','a.MM02.se','a.MM1.se','a.MM12.se',
                              'p.FM0.se','p.FM02.se','p.FM1.se','p.FM12.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model: group differences in actor partner effects - AR(1)
  # with indistinguishable dyads 
  
  if (Model == 28){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c0,length(N0.dyad)),rep(c1,length(N0.dyad)),
                           rep(rho.Y0,length(N0.dyad)),rep(rho.Y1,length(N0.dyad)),
                           rep(a0,length(N0.dyad)),rep(a02,length(N0.dyad)),
                           rep(a1,length(N0.dyad)),rep(a12,length(N0.dyad)),
                           rep(p0,length(N0.dyad)),rep(p02,length(N0.dyad)),
                           rep(p1,length(N0.dyad)),rep(p12,length(N0.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12]) 
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12]) 
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Intercept in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in intercept between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Autoregressive effect','N Dyad','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the autoregressive effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear actor effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic actor effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear actor effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic actor effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Linear partner effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Quadratic partner effect in Group 0','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the linear partner effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Difference in the quadratic partner effect between Group 0 and 1','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N0.dyad)),rep(sigma.eps.M,length(N0.dyad)),
                     rep(rho.eps.FM,length(N0.dyad)),rep(sigma.nu,length(N0.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner A','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev Level 1 error for partner B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i]))),
      unlist(lapply(1:length(N0.dyad), function(i) paste('Std.dev random intercept','N0 Dyad',N0.dyad[i],'N1 Dyad',N1.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N0.dyad,N1.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c1','rho.Y0','rho.Y1','a0','a02','a1','a12','p0','p02','p1','p12',
                              'c0.se','c1.se','rho.Y0.se','rho.Y1.se','a0.se','a02.se','a1.se','a12.se','p0.se','p02.se','p1.se','p12.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model with a continuos time-varying dyad moderator - AR(1)
  
  if (Model == 29){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(rho.YF,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(a.FF2,length(N.dyad)),
                           rep(p.MF,length(N.dyad)),rep(p.MF2,length(N.dyad)),
                           rep(b.F,length(N.dyad)),
                           rep(b.FF,length(N.dyad)),rep(b.FF2,length(N.dyad)),
                           rep(b.MF,length(N.dyad)),rep(b.MF2,length(N.dyad)),
                           rep(rho.YM,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(a.MM2,length(N.dyad)),
                           rep(p.FM,length(N.dyad)),rep(p.FM2,length(N.dyad)),
                           rep(b.M,length(N.dyad)),
                           rep(b.MM,length(N.dyad)),rep(b.MM2,length(N.dyad)),
                           rep(b.FM,length(N.dyad)),rep(b.FM2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16],
                          Coef.hat[,17],Coef.hat[,18],Coef.hat[,19],Coef.hat[,20],Coef.hat[,21],Coef.hat[,22])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16],
                           Coef.bias[,17],Coef.bias[,18],Coef.bias[,19],Coef.bias[,20],Coef.bias[,21],Coef.bias[,22])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16],
                         Coef.se[,17],Coef.se[,18],Coef.se[,19],Coef.se[,20],Coef.se[,21],Coef.se[,22])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16],
                               Coef.CI.Width[,17],Coef.CI.Width[,18],Coef.CI.Width[,19],Coef.CI.Width[,20],Coef.CI.Width[,21],Coef.CI.Width[,22])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16],
                         Coef.CR[,17],Coef.CR[,18],Coef.CR[,19],Coef.CR[,20],Coef.CR[,21],Coef.CR[,22])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16],
                            Coef.power[,17],Coef.power[,18],Coef.power[,19],Coef.power[,20],Coef.power[,21],Coef.power[,22])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16],
                               Coef.power.se[,17],Coef.power.se[,18],Coef.power.se[,19],Coef.power.se[,20],
                               Coef.power.se[,21],Coef.power.se[,22])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear aActor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','rho.YF','a.FF','a.FF2','p.MF','p.MF2','b.F','b.FF','b.FF2','b.MF','b.MF2','rho.YM',
                              'a.MM','a.MM2','p.FM','p.FM2','b.M','b.MM','b.MM2','b.FM','b.FM2',
                              'c.F.se','c.M.se','rho.YF.se','a.FF.se','a.FF2.se','p.MF.se','p.MF2.se','b.F.se','b.FF.se','b.FF2.se','b.MF.se','b.MF2.se','rho.YM.se',
                              'a.MM.se','a.MM2.se','p.FM.se','p.FM2.se','b.M.se','b.MM.se','b.MM2.se','b.FM.se','b.FM2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from the L-APIM model with a continuous time-varying dyad moderator - AR(1)
  # with indistinguishable dyads 
  
  if (Model == 30){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(rho.Y,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(a.2,length(N.dyad)),
                           rep(p,length(N.dyad)),rep(p.2,length(N.dyad)),
                           rep(b,length(N.dyad)),
                           rep(b.a,length(N.dyad)),rep(b.a2,length(N.dyad)),
                           rep(b.p,length(N.dyad)),rep(b.p2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10],Coef.hat[,11])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10],Coef.bias[,11])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10],Coef.se[,11])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10],Coef.CR[,11])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10],Coef.power[,11])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the continuous time-varying variable on the quadratic partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','rho.Y','a','a.2','p','p.2','b','b.a','b.a2','b.p','b.p2',
                              'c.se','rho.Y.se','a.se','a.2.se','p.se','p.2.se','b.se','b.a.se','b.a2.se','b.p.se','b.p2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model with a dichotomous time-varying dyad moderator - AR(1)
  
  if (Model == 31){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c.F,length(N.dyad)),rep(c.M,length(N.dyad)),
                           rep(rho.YF,length(N.dyad)),
                           rep(a.FF,length(N.dyad)),rep(a.FF2,length(N.dyad)),
                           rep(p.MF,length(N.dyad)),rep(p.MF2,length(N.dyad)),
                           rep(d.F,length(N.dyad)),
                           rep(d.FF,length(N.dyad)),rep(d.FF2,length(N.dyad)),
                           rep(d.MF,length(N.dyad)),rep(d.MF2,length(N.dyad)),
                           rep(rho.YM,length(N.dyad)),
                           rep(a.MM,length(N.dyad)),rep(a.MM2,length(N.dyad)),
                           rep(p.FM,length(N.dyad)),rep(p.FM2,length(N.dyad)),
                           rep(d.M,length(N.dyad)),
                           rep(d.MM,length(N.dyad)),rep(d.MM2,length(N.dyad)),
                           rep(d.FM,length(N.dyad)),rep(d.FM2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],
                          Coef.hat[,5],Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],
                          Coef.hat[,9],Coef.hat[,10],Coef.hat[,11],Coef.hat[,12],
                          Coef.hat[,13],Coef.hat[,14],Coef.hat[,15],Coef.hat[,16],
                          Coef.hat[,17],Coef.hat[,18],Coef.hat[,19],Coef.hat[,20],Coef.hat[,21],Coef.hat[,22])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],
                           Coef.bias[,5],Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],
                           Coef.bias[,9],Coef.bias[,10],Coef.bias[,11],Coef.bias[,12],
                           Coef.bias[,13],Coef.bias[,14],Coef.bias[,15],Coef.bias[,16],
                           Coef.bias[,17],Coef.bias[,18],Coef.bias[,19],Coef.bias[,20],Coef.bias[,21],Coef.bias[,22])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],
                         Coef.se[,5],Coef.se[,6],Coef.se[,7],Coef.se[,8],
                         Coef.se[,9],Coef.se[,10],Coef.se[,11],Coef.se[,12],
                         Coef.se[,13],Coef.se[,14],Coef.se[,15],Coef.se[,16],
                         Coef.se[,17],Coef.se[,18],Coef.se[,19],Coef.se[,20],Coef.se[,21],Coef.se[,22])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],
                               Coef.CI.Width[,5],Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],
                               Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11],Coef.CI.Width[,12],
                               Coef.CI.Width[,13],Coef.CI.Width[,14],Coef.CI.Width[,15],Coef.CI.Width[,16],
                               Coef.CI.Width[,17],Coef.CI.Width[,18],Coef.CI.Width[,19],Coef.CI.Width[,20],Coef.CI.Width[,21],Coef.CI.Width[,22])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],
                         Coef.CR[,5],Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],
                         Coef.CR[,9],Coef.CR[,10],Coef.CR[,11],Coef.CR[,12],
                         Coef.CR[,13],Coef.CR[,14],Coef.CR[,15],Coef.CR[,16],
                         Coef.CR[,17],Coef.CR[,18],Coef.CR[,19],Coef.CR[,20],Coef.CR[,21],Coef.CR[,22])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],
                            Coef.power[,5],Coef.power[,6],Coef.power[,7],Coef.power[,8],
                            Coef.power[,9],Coef.power[,10],Coef.power[,11],Coef.power[,12],
                            Coef.power[,13],Coef.power[,14],Coef.power[,15],Coef.power[,16],
                            Coef.power[,17],Coef.power[,18],Coef.power[,19],Coef.power[,20],Coef.power[,21],Coef.power[,22])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],
                               Coef.power.se[,5],Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],
                               Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11],Coef.power.se[,12],
                               Coef.power.se[,13],Coef.power.se[,14],Coef.power.se[,15],Coef.power.se[,16],
                               Coef.power.se[,17],Coef.power.se[,18],Coef.power.se[,19],Coef.power.se[,20],Coef.power.se[,21],Coef.power.se[,22])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic actor effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic partner effect for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic actor effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear partner effect for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic partner effect for partner B','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu.F,length(N.dyad)),
                     rep(sigma.nu.M,length(N.dyad)),rep(rho.nu.F.M,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4],Var.hat[,5],Var.hat[,6])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4],Var.bias[,5],Var.bias[,6])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4],Var.hat.se[,5],Var.hat.se[,6])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between random intercepts for partners A and B','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c.F','c.M','rho.YF','a.FF','a.FF2','p.MF','p.MF2','d.F','d.FF','d.FF2','d.MF','d.MF2',
                              'rho.YM','a.MM','a.MM2','p.FM','p.FM2','d.M','d.MM','d.MM2','d.FM','d.FM2',
                              'c.F.se','c.M.se','rho.YF.se','a.FF.se','a.FF2.se','p.MF.se','p.MF2.se','d.F.se','d.FF.se','d.FF2.se','d.MF.se','d.MF2.se',
                              'rho.YM.se','a.MM.se','a.MM2.se','p.FM.se','p.FM2.se','d.M.se','d.MM.se','d.MM2.se','d.FM.se','d.FM2.se')
  }
  
  ########################################################################################
  ########################################################################################
  
  # Curvilinear actor and partner effects
  # Simulate data from APIM model with a dichotomous time-varying dyad moderator - AR(1)
  # with indistinguishable dyads 
  
  if (Model == 32){
    
    # Table fixed predictors 
    
    coef.sim.True.list = c(rep(c,length(N.dyad)),rep(rho.Y,length(N.dyad)),
                           rep(a,length(N.dyad)),rep(a.2,length(N.dyad)),
                           rep(p,length(N.dyad)),rep(p.2,length(N.dyad)),
                           rep(d,length(N.dyad)),
                           rep(d.a,length(N.dyad)),rep(d.a2,length(N.dyad)),
                           rep(d.p,length(N.dyad)),rep(d.p2,length(N.dyad)))
    coef.sim.coef.hat = c(Coef.hat[,1],Coef.hat[,2],Coef.hat[,3],Coef.hat[,4],Coef.hat[,5],
                          Coef.hat[,6],Coef.hat[,7],Coef.hat[,8],Coef.hat[,9],Coef.hat[,10],Coef.hat[,11])
    coef.sim.coef.bias = c(Coef.bias[,1],Coef.bias[,2],Coef.bias[,3],Coef.bias[,4],Coef.bias[,5],
                           Coef.bias[,6],Coef.bias[,7],Coef.bias[,8],Coef.bias[,9],Coef.bias[,10],Coef.bias[,11])
    coef.sim.coef.se = c(Coef.se[,1],Coef.se[,2],Coef.se[,3],Coef.se[,4],Coef.se[,5],
                         Coef.se[,6],Coef.se[,7],Coef.se[,8],Coef.se[,9],Coef.se[,10],Coef.se[,11])
    coef.sim.coef.CI.Width = c(Coef.CI.Width[,1],Coef.CI.Width[,2],Coef.CI.Width[,3],Coef.CI.Width[,4],Coef.CI.Width[,5],
                               Coef.CI.Width[,6],Coef.CI.Width[,7],Coef.CI.Width[,8],Coef.CI.Width[,9],Coef.CI.Width[,10],Coef.CI.Width[,11])
    coef.sim.coef.CR = c(Coef.CR[,1],Coef.CR[,2],Coef.CR[,3],Coef.CR[,4],Coef.CR[,5],
                         Coef.CR[,6],Coef.CR[,7],Coef.CR[,8],Coef.CR[,9],Coef.CR[,10],Coef.CR[,11])
    coef.sim.coef.power = c(Coef.power[,1],Coef.power[,2],Coef.power[,3],Coef.power[,4],Coef.power[,5],
                            Coef.power[,6],Coef.power[,7],Coef.power[,8],Coef.power[,9],Coef.power[,10],Coef.power[,11])
    coef.sim.coef.power.se = c(Coef.power.se[,1],Coef.power.se[,2],Coef.power.se[,3],Coef.power.se[,4],Coef.power.se[,5],
                               Coef.power.se[,6],Coef.power.se[,7],Coef.power.se[,8],Coef.power.se[,9],Coef.power.se[,10],Coef.power.se[,11])
    
    coef.sim = cbind(coef.sim.True.list,coef.sim.coef.hat,coef.sim.coef.bias,coef.sim.coef.se,
                     coef.sim.coef.CI.Width,coef.sim.coef.CR,coef.sim.coef.power,coef.sim.coef.power.se)
    coef.sim = data.frame(round(coef.sim,digits=4))
    
    rownames(coef.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Intercept','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Autoregressive effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Quadratic partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic actor effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the linear partner effect','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Effect of the dichotomous time-varying variable on the quadratic partner effect','N Dyad',N.dyad[i]))))
    
    colnames(coef.sim) = c('True','Mean','Bias','Std.Error','CI width','(1-alpha)% Coverage Rate','Power','Std.Error Power')
    
    # Table random effects 
    
    cov.sim.True = c(rep(sigma.eps.F,length(N.dyad)),rep(sigma.eps.M,length(N.dyad)),
                     rep(rho.eps.FM,length(N.dyad)),rep(sigma.nu,length(N.dyad)))
    cov.sim.hat = c(Var.hat[,1],Var.hat[,2],Var.hat[,3],Var.hat[,4])
    cov.sim.bias = c(Var.bias[,1],Var.bias[,2],Var.bias[,3],Var.bias[,4])
    cov.sim.se = c(Var.hat.se[,1],Var.hat.se[,2],Var.hat.se[,3],Var.hat.se[,4])
    
    cov.sim = cbind(cov.sim.True,cov.sim.hat,cov.sim.bias,cov.sim.se)
    cov.sim = data.frame(round(cov.sim,digits=4))
    
    rownames(cov.sim) = c(
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner A','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev Level 1 error for partner B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Cor between Level 1 error for partners A and B','N Dyad',N.dyad[i]))),
      unlist(lapply(1:length(N.dyad), function(i) paste('Std.dev random intercept','N Dyad',N.dyad[i])))) 
    
    colnames(cov.sim) = c('True','Mean','Bias','Std.Error')
    
    # Power (for the power curve!)
    
    power.curve = cbind(N.dyad,Coef.power,Coef.power.se)
    colnames(power.curve) = c('N dyad','c','rho.Y','a','a.2','p','p.2','d','d.a','d.a2','d.p','d.p2',
                              'c.se','rho.Y.se','a.se','a.2.se','p.se','p.2.se','d.se','d.a.se','d.a2.se','d.p.se','d.p2.se')
  }
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  return(list(replicates=replicates,power.curve=power.curve,coef.sim=coef.sim,cov.sim=cov.sim))}

#####################################################################################

