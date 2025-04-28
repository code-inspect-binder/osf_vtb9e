Sim.Dyad.Model.5 = function(N.dyad,T.obs,  
c.F,c.M,a.FF,p.MF,a.MM,p.FM,
b.F,b.M,b.FF,b.MF,b.MM,b.FM,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W){

# Number of subject
N.subject = 2*N.dyad

# Create variables observations and subjects
data.Model = cbind(expand.grid(Obs=1:T.obs,subject.ID=1:N.subject),
dyad.ID=unlist(lapply(1:N.dyad, function(i) rep(i,2*T.obs))),Gender=rep(c(rep('F',T.obs),rep('M',T.obs)),N.dyad),
Female=rep(c(rep(1,T.obs),rep(0,T.obs)),N.dyad),Male=rep(c(rep(0,T.obs),rep(1,T.obs)),N.dyad))

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
var.diag.X = c(sigma.XF,sigma.XM)
Sigma.X = diag(length(var.diag.X))
Sigma.X[lower.tri(Sigma.X, diag=FALSE)] = rho.X
Sigma.X = pmax(Sigma.X, t(Sigma.X), na.rm=TRUE)
Sigma.X = diag(var.diag.X)%*%Sigma.X%*%diag(var.diag.X)
data.X = cbind(data.X,mvrnorm(N.dyad*T.obs,c(mu.XF,mu.XM),Sigma.X))
colnames(data.X) = c('Obs','ID','X.F','X.M')
X.F = data.X[,'X.F']
X.M = data.X[,'X.M']
W.dyad = rnorm(N.dyad*T.obs,mu.W,sigma.W)
data.X = cbind(data.X,W.dyad)

# If is.center.X is equal to TRUE, person-mean centered the predictors
if(is.center.X==TRUE){
if(is.center.W==TRUE){  
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M), 
       W.dyad.c = W.dyad - mean(W.dyad))
# Simulate Dependent Variables
Y.F = c.F + a.FF*data.X[,'X.F.c'] + p.MF*data.X[,'X.M.c'] + b.F*data.X[,'W.dyad.c'] + b.FF*data.X[,'X.F.c']*data.X[,'W.dyad.c'] + b.MF*data.X[,'X.M.c']*data.X[,'W.dyad.c'] + V[,'V.F'] + E[,'E.F']

Y.M = c.M + a.MM*data.X[,'X.M.c'] + p.FM*data.X[,'X.F.c'] + b.M*data.X[,'W.dyad.c'] + b.MM*data.X[,'X.M.c']*data.X[,'W.dyad.c'] + b.FM*data.X[,'X.F.c']*data.X[,'W.dyad.c'] + V[,'V.M'] + E[,'E.M']
}
if(is.center.W==FALSE){
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M))
# Simulate Dependent Variables
Y.F = c.F + a.FF*data.X[,'X.F.c'] + p.MF*data.X[,'X.M.c'] + b.F*data.X[,'W.dyad'] + b.FF*data.X[,'X.F.c']*data.X[,'W.dyad'] + b.MF*data.X[,'X.M.c']*data.X[,'W.dyad'] + V[,'V.F'] + E[,'E.F']

Y.M = c.M + a.MM*data.X[,'X.M.c'] + p.FM*data.X[,'X.F.c'] + b.M*data.X[,'W.dyad'] + b.MM*data.X[,'X.M.c']*data.X[,'W.dyad'] + b.FM*data.X[,'X.F.c']*data.X[,'W.dyad'] + V[,'V.M'] + E[,'E.M']  
}}  

if(is.center.X==FALSE){
if(is.center.W==TRUE){  
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(W.dyad.c = W.dyad - mean(W.dyad))
# Simulate Dependent Variables
Y.F = c.F + a.FF*data.X[,'X.F'] + p.MF*data.X[,'X.M'] + b.F*data.X[,'W.dyad.c'] + b.FF*data.X[,'X.F']*data.X[,'W.dyad.c'] + b.MF*data.X[,'X.M']*data.X[,'W.dyad.c'] + V[,'V.F'] + E[,'E.F']

Y.M = c.M + a.MM*data.X[,'X.M'] + p.FM*data.X[,'X.F'] + b.M*data.X[,'W.dyad.c'] + b.MM*data.X[,'X.M']*data.X[,'W.dyad.c'] + b.FM*data.X[,'X.F']*data.X[,'W.dyad.c'] + V[,'V.M'] + E[,'E.M']
}
if(is.center.W==FALSE){
# Simulate Dependent Variables
Y.F = c.F + a.FF*data.X[,'X.F'] + p.MF*data.X[,'X.M'] + b.F*data.X[,'W.dyad'] + b.FF*data.X[,'X.F']*data.X[,'W.dyad'] + b.MF*data.X[,'X.M']*data.X[,'W.dyad'] + V[,'V.F'] + E[,'E.F']

Y.M = c.M + a.MM*data.X[,'X.M'] + p.FM*data.X[,'X.F'] + b.M*data.X[,'W.dyad'] + b.MM*data.X[,'X.M']*data.X[,'W.dyad'] + b.FM*data.X[,'X.F']*data.X[,'W.dyad'] + V[,'V.M'] + E[,'E.M']
}}

# Create a data frame
X.Actor = rep(0,nrow(data.Model))
X.Partner = rep(0,nrow(data.Model))
W = rep(0,nrow(data.Model))
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
W[N.id[which(data.Model[N.id,]$Gender=='F')]] = W.dyad[chunk[[i]]]
W[N.id[which(data.Model[N.id,]$Gender=='M')]] = W.dyad[chunk[[i]]]
}

# Create a data frame
data.Model = data.frame(cbind(data.Model,Y,X.Actor,X.Partner,W)) 

return(data=data.Model)
}