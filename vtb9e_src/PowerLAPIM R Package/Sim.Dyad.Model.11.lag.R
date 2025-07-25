Sim.Dyad.Model.11.lag = function(N0.dyad,N1.dyad,T.obs,  
c.F0,c.F1,c.M0,c.M1,rho.YF0,rho.YF1,rho.YM0,rho.YM1,
a.FF0,a.FF1,a.FF02,a.FF12,p.MF0,p.MF1,p.MF02,p.MF12,
a.MM0,a.MM1,a.MM02,a.MM12,p.FM0,p.FM1,p.FM02,p.FM12,
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

# Create number of observations: T.obs + T.burning
T.burning = 1000
T.total = T.burning + T.obs

# Create variables observations and subjects
data.Model = cbind(expand.grid(Obs=1:T.obs,subject.ID=1:N.subject),
Z=c(rep(0,T.obs*N0.subject),rep(1,T.obs*N1.subject)),dyad.ID=c(unlist(lapply(1:N0.dyad, function(i) rep(i,2*T.obs))),unlist(lapply((N0.dyad+1):N.dyad, function(i) rep(i,2*T.obs)))),Gender=c(rep(c(rep('F',T.obs),rep('M',T.obs)),N0.dyad),rep(c(rep('F',T.obs),rep('M',T.obs)),N1.dyad)),Female=c(rep(c(rep(1,T.obs),rep(0,T.obs)),N0.dyad),rep(c(rep(1,T.obs),rep(0,T.obs)),N1.dyad)),Male=c(rep(c(rep(0,T.obs),rep(1,T.obs)),N0.dyad),rep(c(rep(0,T.obs),rep(1,T.obs)),N1.dyad)))

# Simulate error within-person errors

Sigma.eps = diag(2)
Sigma.eps[lower.tri(Sigma.eps, diag=FALSE)] = rho.eps.FM
Sigma.eps = pmax(Sigma.eps, t(Sigma.eps), na.rm=TRUE)
Sigma.eps = diag(c(sigma.eps.F,sigma.eps.M))%*%Sigma.eps%*%diag(c(sigma.eps.F,sigma.eps.M))
E = mvrnorm(T.total*N.dyad, mu=c(0,0), Sigma.eps)
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
V = rbind(V,matrix(unlist(lapply(1:length(var.diag.nu), function(p) rep(V.j[j,p],T.total))), ncol=length(var.diag.nu), byrow=F))
}
colnames(V) = c('V.F','V.M')

# Simulate time varying variable X
data.X = expand.grid(Obs=1:T.total,ID=1:N.dyad)
var.diag.X0 = c(sigma.XF0,sigma.XM0)
Sigma.X0 = diag(length(var.diag.X0))
Sigma.X0[lower.tri(Sigma.X0, diag=FALSE)] = rho.X0
Sigma.X0 = pmax(Sigma.X0, t(Sigma.X0), na.rm=TRUE)
Sigma.X0 = diag(var.diag.X0)%*%Sigma.X0%*%diag(var.diag.X0)
data.X0 = mvrnorm(N0.dyad*T.total,c(mu.XF0,mu.XM0),Sigma.X0)
var.diag.X1 = c(sigma.XF1,sigma.XM1)
Sigma.X1 = diag(length(var.diag.X1))
Sigma.X1[lower.tri(Sigma.X1, diag=FALSE)] = rho.X1
Sigma.X1 = pmax(Sigma.X1, t(Sigma.X1), na.rm=TRUE)
Sigma.X1 = diag(var.diag.X1)%*%Sigma.X1%*%diag(var.diag.X1)
data.X1 = mvrnorm(N1.dyad*T.total,c(mu.XF1,mu.XM1),Sigma.X1)
X.F0 = data.X0[,1]
X.F1 = data.X1[,1]
X.M0 = data.X0[,2]
X.M1 = data.X1[,2]
X.F=c(X.F0,X.F1)
X.M=c(X.M0,X.M1)
Z=c(rep(0,N0.dyad*T.total),rep(1,N0.dyad*T.total))
data.X = cbind(data.X,Z,X.F,X.M)

# If is.center.X is equal to TRUE, person-mean centered the predictors
if(is.center.X==TRUE){
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M))
data.X = data.frame(data.X)
}

# Function to get recursive equation
Y.F = rep(0,nrow(data.X))
Y.M = rep(0,nrow(data.X))
n.ID = unique(data.X$ID)
for (i in n.ID){
T.obs.i = which(data.X$ID==i) 

if (is.center.X == TRUE){
# Initialized values
Y.F[T.obs.i[1]] = c.F0 + c.F1*data.X[T.obs.i[1],'Z'] + a.FF0*data.X[T.obs.i[1],'X.F.c'] + a.FF1*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'Z'] + a.FF02*I(data.X[T.obs.i[1],'X.F.c']^2) + a.FF12*I(data.X[T.obs.i[1],'X.F.c']^2)*data.X[T.obs.i[1],'Z'] + p.MF0*data.X[T.obs.i[1],'X.M.c'] + p.MF1*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'Z'] + p.MF02*I(data.X[T.obs.i[1],'X.M.c']^2) + p.MF12*I(data.X[T.obs.i[1],'X.M.c']^2)*data.X[T.obs.i[1],'Z'] + V[T.obs.i[1],'V.F'] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c.M0 + c.M1*data.X[T.obs.i[1],'Z'] + a.MM0*data.X[T.obs.i[1],'X.M.c'] + a.MM1*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'Z'] + a.MM02*I(data.X[T.obs.i[1],'X.M.c']^2) + a.MM12*I(data.X[T.obs.i[1],'X.M.c']^2)*data.X[T.obs.i[1],'Z'] + p.FM0*data.X[T.obs.i[1],'X.F.c'] + p.FM1*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'Z'] + p.FM02*I(data.X[T.obs.i[1],'X.F.c']^2) + p.FM12*I(data.X[T.obs.i[1],'X.F.c']^2)*data.X[T.obs.i[1],'Z'] + V[T.obs.i[1],'V.M'] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c.F0 + c.F1*data.X[t,'Z'] + rho.YF0*Y.F[t-1] + rho.YF1*Y.F[t-1]*data.X[t,'Z'] + a.FF0*data.X[t,'X.F.c'] + a.FF1*data.X[t,'X.F.c']*data.X[t,'Z'] + a.FF02*I(data.X[t,'X.F.c']^2) + a.FF12*I(data.X[t,'X.F.c']^2)*data.X[t,'Z'] + p.MF0*data.X[t,'X.M.c'] + p.MF1*data.X[t,'X.M.c']*data.X[t,'Z'] + p.MF02*I(data.X[t,'X.M.c']^2) + p.MF12*I(data.X[t,'X.M.c']^2)*data.X[t,'Z'] + V[t,'V.F'] + E[t,'E.F']

Y.M[t] = c.M0 + c.M1*data.X[t,'Z'] + rho.YM0*Y.M[t-1] + rho.YM1*Y.M[t-1]*data.X[t,'Z'] + a.MM0*data.X[t,'X.M.c'] + a.MM1*data.X[t,'X.M.c']*data.X[t,'Z'] + a.MM02*I(data.X[t,'X.M.c']^2) + a.MM12*I(data.X[t,'X.M.c']^2)*data.X[t,'Z'] + p.FM0*data.X[t,'X.F.c'] + p.FM1*data.X[t,'X.F.c']*data.X[t,'Z'] + p.FM02*I(data.X[t,'X.F.c']^2) + p.FM12*I(data.X[t,'X.F.c']^2)*data.X[t,'Z'] + V[t,'V.M'] + E[t,'E.M']
}}

if (is.center.X == FALSE){
# Initialized values
Y.F[T.obs.i[1]] = c.F0 + c.F1*data.X[T.obs.i[1],'Z'] + a.FF0*data.X[T.obs.i[1],'X.F'] + a.FF1*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'Z'] + a.FF02*I(data.X[T.obs.i[1],'X.F']^2) + a.FF12*I(data.X[T.obs.i[1],'X.F']^2)*data.X[T.obs.i[1],'Z'] + p.MF0*data.X[T.obs.i[1],'X.M'] + p.MF1*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'Z'] + p.MF02*I(data.X[T.obs.i[1],'X.M']^2) + p.MF12*I(data.X[T.obs.i[1],'X.M']^2)*data.X[T.obs.i[1],'Z'] + V[T.obs.i[1],'V.F'] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c.M0 + c.M1*data.X[T.obs.i[1],'Z'] + a.MM0*data.X[T.obs.i[1],'X.M'] + a.MM1*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'Z'] + a.MM02*I(data.X[T.obs.i[1],'X.M']^2) + a.MM12*I(data.X[T.obs.i[1],'X.M']^2)*data.X[T.obs.i[1],'Z'] + p.FM0*data.X[T.obs.i[1],'X.F'] + p.FM1*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'Z'] + p.FM02*I(data.X[T.obs.i[1],'X.F']^2) + p.FM12*I(data.X[T.obs.i[1],'X.F']^2)*data.X[T.obs.i[1],'Z'] + V[T.obs.i[1],'V.M'] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c.F0 + c.F1*data.X[t,'Z'] + rho.YF0*Y.F[t-1] + rho.YF1*Y.F[t-1]*data.X[t,'Z'] + a.FF0*data.X[t,'X.F'] + a.FF1*data.X[t,'X.F']*data.X[t,'Z'] + a.FF02*I(data.X[t,'X.F']^2) + a.FF12*I(data.X[t,'X.F']^2)*data.X[t,'Z'] + p.MF0*data.X[t,'X.M'] + p.MF1*data.X[t,'X.M']*data.X[t,'Z'] + p.MF02*I(data.X[t,'X.M']^2) + p.MF12*I(data.X[t,'X.M']^2)*data.X[t,'Z'] + V[t,'V.F'] + E[t,'E.F']

Y.M[t] = c.M0 + c.M1*data.X[t,'Z'] + rho.YM0*Y.M[t-1] + rho.YM1*Y.M[t-1]*data.X[t,'Z'] + a.MM0*data.X[t,'X.M'] + a.MM1*data.X[t,'X.M']*data.X[t,'Z'] + a.MM02*I(data.X[t,'X.M']^2) + a.MM12*I(data.X[t,'X.M']^2)*data.X[t,'Z'] + p.FM0*data.X[t,'X.F'] + p.FM1*data.X[t,'X.F']*data.X[t,'Z'] + p.FM02*I(data.X[t,'X.F']^2) + p.FM12*I(data.X[t,'X.F']^2)*data.X[t,'Z'] + V[t,'V.M'] + E[t,'E.M']
}}}

data.Y = cbind(Y.F,Y.M)
colnames(data.Y) = c('Y.F','Y.M')

T.total.i = NULL
for (i in n.ID){
T.total.i = c(T.total.i,which(data.X$ID==i)[-seq(1:T.burning)])  
}

# Create a data frame for T.obs
data.X = data.X[T.total.i,]
data.Y = data.Y[T.total.i,]

Y.F = data.Y[,'Y.F']
Y.M = data.Y[,'Y.M']
X.F = as.matrix(data.X[,'X.F'])
X.M = as.matrix(data.X[,'X.M'])

# Create a data frame
X.Actor = rep(0,nrow(data.Model))
X.Partner = rep(0,nrow(data.Model))
Y = rep(0,nrow(data.Model))

chunk = split(1:(N.dyad*T.obs), factor(sort(rank(1:(N.dyad*T.obs))%%N.dyad)))

for (i in 1:N.dyad){
N.id = which(data.Model$dyad.ID==i)
Y[N.id[which(data.Model[N.id,]$Gender=='F')]] = Y.F[chunk[[i]]]
Y[N.id[which(data.Model[N.id,]$Gender=='M')]] = Y.M[chunk[[i]]]
X.Actor[N.id[which(data.Model[N.id,]$Gender=='F')]] = X.F[chunk[[i]]]
X.Actor[N.id[which(data.Model[N.id,]$Gender=='M')]] = X.M[chunk[[i]]]
X.Partner[N.id[which(data.Model[N.id,]$Gender=='F')]] = X.M[chunk[[i]]]
X.Partner[N.id[which(data.Model[N.id,]$Gender=='M')]] = X.F[chunk[[i]]]
}

# Create a data frame
data.Model = data.frame(cbind(data.Model,Y,X.Actor,X.Partner)) 

# Create lag variable
Y.lag = rep(0,nrow(data.Model))
n.subject = unique(data.Model$subject.ID)
for (j in n.subject){
Y.lag[which(data.Model$subject.ID==j)] = shift(data.Model$Y[which(data.Model$subject.ID==j)])
}

data.Model = cbind(data.Model,Y.lag=Y.lag)

return(data=data.Model)
}