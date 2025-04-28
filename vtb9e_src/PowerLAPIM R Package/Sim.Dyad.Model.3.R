Sim.Dyad.Model.3 = function(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,a.FF0,a.FF1,p.MF0,p.MF1,a.MM0,a.MM1,p.FM0,p.FM1,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF0,mu.XF1,sigma.XF0,sigma.XF1,mu.XM0,mu.XM1,
sigma.XM0,sigma.XM1,rho.X0,rho.X1,is.center.X){

# Number of dyads
N.dyad = N0.dyad + N1.dyad  
    
# Number of subject in Group 0
N0.subject = 2*N0.dyad
# Number of subject in Group 1
N1.subject = 2*N1.dyad

N.subject = N0.subject + N1.subject

# Create variables observations and subjects
data.Model = cbind(expand.grid(Obs=1:T.obs,subject.ID=1:N.subject),
Z=c(rep(0,T.obs*N0.subject),rep(1,T.obs*N1.subject)),dyad.ID=c(unlist(lapply(1:N0.dyad, function(i) rep(i,2*T.obs))),unlist(lapply((N0.dyad+1):N.dyad, function(i) rep(i,2*T.obs)))),Gender=c(rep(c(rep('F',T.obs),rep('M',T.obs)),N0.dyad),rep(c(rep('F',T.obs),rep('M',T.obs)),N1.dyad)),Female=c(rep(c(rep(1,T.obs),rep(0,T.obs)),N0.dyad),rep(c(rep(1,T.obs),rep(0,T.obs)),N1.dyad)),Male=c(rep(c(rep(0,T.obs),rep(1,T.obs)),N0.dyad),rep(c(rep(0,T.obs),rep(1,T.obs)),N1.dyad)))

# Simulate error within-person errors

Sigma.eps = diag(2)
Sigma.eps[lower.tri(Sigma.eps, diag=FALSE)] = rho.eps.FM
Sigma.eps = pmax(Sigma.eps, t(Sigma.eps), na.rm=TRUE)
Sigma.eps = diag(c(sigma.eps.F,sigma.eps.M))%*%Sigma.eps%*%diag(c(sigma.eps.F,sigma.eps.M))
E = mvrnorm(T.obs*N.dyad, mu=c(0,0), Sigma.eps)
colnames(E) = c('E.F','E.M')

# Simulate error level-2
# Simulate between-subject random effect
var.diag.nu = c(sigma.nu.F,sigma.nu.M)
Sigma.nu = diag(length(var.diag.nu))
Sigma.nu[lower.tri(Sigma.nu, diag=FALSE)] = rho.nu.F.M
Sigma.nu = pmax(Sigma.nu, t(Sigma.nu), na.rm=TRUE)
Sigma.nu = diag(var.diag.nu)%*%Sigma.nu%*%diag(var.diag.nu)
V.j = mvrnorm(N.dyad,rep(0,ncol(Sigma.nu)),Sigma.nu)
colnames(V.j) = c('V.F','V.M')

V = NULL
for (j in 1:N.dyad){
V = rbind(V,matrix(unlist(lapply(1:length(var.diag.nu), function(p) rep(V.j[j,p],T.obs))), ncol=length(var.diag.nu), byrow=F))
}
colnames(V) = c('V.F','V.M')

# Simulate time varying variable X
data.X = expand.grid(Obs=1:T.obs,ID=1:N.dyad)
var.diag.X0 = c(sigma.XF0,sigma.XM0)
Sigma.X0 = diag(length(var.diag.X0))
Sigma.X0[lower.tri(Sigma.X0, diag=FALSE)] = rho.X0
Sigma.X0 = pmax(Sigma.X0, t(Sigma.X0), na.rm=TRUE)
Sigma.X0 = diag(var.diag.X0)%*%Sigma.X0%*%diag(var.diag.X0)
data.X0 = mvrnorm(N0.dyad*T.obs,c(mu.XF0,mu.XM0),Sigma.X0)
var.diag.X1 = c(sigma.XF1,sigma.XM1)
Sigma.X1 = diag(length(var.diag.X1))
Sigma.X1[lower.tri(Sigma.X1, diag=FALSE)] = rho.X1
Sigma.X1 = pmax(Sigma.X1, t(Sigma.X1), na.rm=TRUE)
Sigma.X1 = diag(var.diag.X1)%*%Sigma.X1%*%diag(var.diag.X1)
data.X1 = mvrnorm(N1.dyad*T.obs,c(mu.XF1,mu.XM1),Sigma.X1)
X.F0 = data.X0[,1]
X.F1 = data.X1[,1]
X.M0 = data.X0[,2]
X.M1 = data.X1[,2]
X.F=c(X.F0,X.F1)
X.M=c(X.M0,X.M1)
Z=c(rep(0,N0.dyad*T.obs),rep(1,N0.dyad*T.obs))
data.X = cbind(data.X,Z,X.F,X.M)

# If is.center.X is equal to TRUE, person-mean centered the predictors
if(is.center.X==TRUE){
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M))
# Simulate Dependent Variables
Y.F = c.F0 + c.F1*data.X[,'Z'] + a.FF0*data.X[,'X.F.c'] + a.FF1*data.X[,'X.F.c']*data.X[,'Z'] + p.MF0*data.X[,'X.M.c'] + p.MF1*data.X[,'X.M.c']*data.X[,'Z'] + V[,'V.F'] + E[,'E.F']

Y.M = c.M0 + c.M1*data.X[,'Z'] + a.MM0*data.X[,'X.M.c'] + a.MM1*data.X[,'X.M.c']*data.X[,'Z'] + p.FM0*data.X[,'X.F.c'] + p.FM1*data.X[,'X.F.c']*data.X[,'Z'] + V[,'V.M'] + E[,'E.M']
}

if(is.center.X==FALSE){
# Simulate Dependent Variables
Y.F = c.F0 + c.F1*data.X[,'Z'] + a.FF0*data.X[,'X.F'] + a.FF1*data.X[,'X.F']*data.X[,'Z'] + p.MF0*data.X[,'X.M'] + p.MF1*data.X[,'X.M']*data.X[,'Z'] + V[,'V.F'] + E[,'E.F']

Y.M = c.M0 + c.M1*data.X[,'Z'] + a.MM0*data.X[,'X.M'] + a.MM1*data.X[,'X.M']*data.X[,'Z'] + p.FM0*data.X[,'X.F'] + p.FM1*data.X[,'X.F']*data.X[,'Z'] + V[,'V.M'] + E[,'E.M']
}

# Create a data frame
X.Actor = rep(0,nrow(data.Model))
X.Partner = rep(0,nrow(data.Model))
Y = rep(0,nrow(data.Model))

chunk = split(1:(N.dyad*T.obs), factor(sort(rank(1:(N.dyad*T.obs))%%N.dyad)))

for (i in 1:N.dyad){
N.id = which(data.Model$dyad.ID==i)
Y[N.id[which(data.Model[N.id,]$Gender=='F')]] = Y.F[chunk[[i]],]
Y[N.id[which(data.Model[N.id,]$Gender=='M')]] = Y.M[chunk[[i]],]
X.Actor[N.id[which(data.Model[N.id,]$Gender=='F')]] = X.F[chunk[[i]]]
X.Actor[N.id[which(data.Model[N.id,]$Gender=='M')]] = X.M[chunk[[i]]]
X.Partner[N.id[which(data.Model[N.id,]$Gender=='F')]] = X.M[chunk[[i]]]
X.Partner[N.id[which(data.Model[N.id,]$Gender=='M')]] = X.F[chunk[[i]]]
}  

# Create a data frame
data.Model = data.frame(cbind(data.Model,Y,X.Actor,X.Partner)) 

return(data=data.Model)
}