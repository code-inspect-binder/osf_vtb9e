Sim.Dyad.Model.6.lag = function(N.dyad,T.obs,  
c,rho.Y,a,p,b,b.a,b.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
mu.W,sigma.W,is.center.X,is.center.W){

# Number of subject
N.subject = 2*N.dyad

# Create number of observations: T.obs + T.burning
T.burning = 1000
T.total = T.burning + T.obs

# Create variables observations and subjects
data.Model = cbind(expand.grid(Obs=1:T.obs,subject.ID=1:N.subject),
dyad.ID=unlist(lapply(1:N.dyad, function(i) rep(i,2*T.obs))),Gender=rep(c(rep('F',T.obs),rep('M',T.obs)),N.dyad),
Female=rep(c(rep(1,T.obs),rep(0,T.obs)),N.dyad),Male=rep(c(rep(0,T.obs),rep(1,T.obs)),N.dyad))

# Simulate error within-person errors

Sigma.eps = diag(2)
Sigma.eps[lower.tri(Sigma.eps, diag=FALSE)] = rho.eps.FM
Sigma.eps = pmax(Sigma.eps, t(Sigma.eps), na.rm=TRUE)
Sigma.eps = diag(c(sigma.eps.F,sigma.eps.M))%*%Sigma.eps%*%diag(c(sigma.eps.F,sigma.eps.M))
E = mvrnorm(T.total*N.dyad, mu=c(0,0), Sigma.eps)
colnames(E) = c('E.F','E.M')

# Simulate error level-2
# Simulate between-subject random effect
V.j = rnorm(N.dyad,0,sigma.nu)

V = NULL
for (j in 1:N.dyad){
V = c(V,rep(V.j[j],T.total))
}

# Simulate time varying variable X
data.X = expand.grid(Obs=1:T.total,ID=1:N.dyad)
var.diag.X = c(sigma.XF,sigma.XM)
Sigma.X = diag(length(var.diag.X))
Sigma.X[lower.tri(Sigma.X, diag=FALSE)] = rho.X
Sigma.X = pmax(Sigma.X, t(Sigma.X), na.rm=TRUE)
Sigma.X = diag(var.diag.X)%*%Sigma.X%*%diag(var.diag.X)
data.X = cbind(data.X,mvrnorm(N.dyad*T.total,c(mu.XF,mu.XM),Sigma.X))
colnames(data.X) = c('Obs','ID','X.F','X.M')
X.F = data.X[,'X.F']
X.M = data.X[,'X.M']
W.dyad = rnorm(N.dyad*T.total,mu.W,sigma.W)
data.X = cbind(data.X,W.dyad)

# If is.center.X is equal to TRUE, person-mean centered the predictors
if(is.center.X==TRUE){
if (is.center.W==TRUE){
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M), 
       W.dyad.c = W.dyad - mean(W.dyad))
data.X = data.frame(data.X)
}
if (is.center.W==FALSE){
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M))
data.X = data.frame(data.X)
}}

# Function to get recursive equation
Y.F = rep(0,nrow(data.X))
Y.M = rep(0,nrow(data.X))
n.ID = unique(data.X$ID)
for (i in n.ID){
T.obs.i = which(data.X$ID==i) 

if (is.center.X == TRUE){
if (is.center.W==TRUE){
# Initialized values
Y.F[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.F.c'] + p*data.X[T.obs.i[1],'X.M.c'] + b*data.X[T.obs.i[1],'W.dyad.c'] + b.a*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'W.dyad.c'] + b.p*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'W.dyad.c'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.M.c'] + p*data.X[T.obs.i[1],'X.F.c'] + b*data.X[T.obs.i[1],'W.dyad.c'] + b.a*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'W.dyad.c'] + b.p*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'W.dyad.c'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c + rho.Y*Y.F[t-1] + a*data.X[t,'X.F.c'] + p*data.X[t,'X.M.c'] + b*data.X[t,'W.dyad.c'] + b.a*data.X[t,'X.F.c']*data.X[t,'W.dyad.c'] + b.p*data.X[t,'X.M.c']*data.X[t,'W.dyad.c'] + V[t] + E[t,'E.F']

Y.M[t] = c + rho.Y*Y.M[t-1] + a*data.X[t,'X.M.c'] + p*data.X[t,'X.F.c'] + b*data.X[t,'W.dyad.c'] + b.a*data.X[t,'X.M.c']*data.X[t,'W.dyad.c'] + b.p*data.X[t,'X.F.c']*data.X[t,'W.dyad.c'] + V[t] + E[t,'E.M']
}}
if (is.center.W==FALSE){
# Initialized values
Y.F[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.F.c'] + p*data.X[T.obs.i[1],'X.M.c'] + b*data.X[T.obs.i[1],'W.dyad'] + b.a*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'W.dyad'] + b.p*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'W.dyad'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.M.c'] + p*data.X[T.obs.i[1],'X.F.c'] + b*data.X[T.obs.i[1],'W.dyad'] + b.a*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'W.dyad'] + b.p*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'W.dyad'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c + rho.Y*Y.F[t-1] + a*data.X[t,'X.F.c'] + p*data.X[t,'X.M.c'] +  b*data.X[t,'W.dyad'] + b.a*data.X[t,'X.F.c']*data.X[t,'W.dyad'] + b.p*data.X[t,'X.M.c']*data.X[t,'W.dyad'] + V[t] + E[t,'E.F']

Y.M[t] = c + rho.Y*Y.M[t-1] + a*data.X[t,'X.M.c'] + p*data.X[t,'X.F.c'] + b*data.X[t,'W.dyad'] + b.a*data.X[t,'X.M.c']*data.X[t,'W.dyad'] + b.p*data.X[t,'X.F.c']*data.X[t,'W.dyad'] + V[t] + E[t,'E.M']
}}}


if (is.center.X == FALSE){
if (is.center.W==TRUE){
# Initialized values
Y.F[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.F'] + p*data.X[T.obs.i[1],'X.M'] + b*data.X[T.obs.i[1],'W.dyad.c'] + b.a*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'W.dyad.c'] + b.p*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'W.dyad.c'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.M'] + p*data.X[T.obs.i[1],'X.F'] + b*data.X[T.obs.i[1],'W.dyad.c'] + b.a*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'W.dyad.c'] + b.p*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'W.dyad.c'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c + rho.Y*Y.F[t-1] + a*data.X[t,'X.F'] + p*data.X[t,'X.M'] + b*data.X[t,'W.dyad.c'] + b.a*data.X[t,'X.F']*data.X[t,'W.dyad.c'] + b.p*data.X[t,'X.M']*data.X[t,'W.dyad.c'] + V[t] + E[t,'E.F']

Y.M[t] = c + rho.Y*Y.M[t-1] + a*data.X[t,'X.M'] + p*data.X[t,'X.F'] + b*data.X[t,'W.dyad.c'] + b.a*data.X[t,'X.M']*data.X[t,'W.dyad.c'] + b.p*data.X[t,'X.F']*data.X[t,'W.dyad.c'] + V[t] + E[t,'E.M']
}}
if (is.center.W==FALSE){
# Initialized values
Y.F[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.F'] + p*data.X[T.obs.i[1],'X.M'] + b*data.X[T.obs.i[1],'W.dyad'] + b.a*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'W.dyad'] + b.p*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'W.dyad'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c + a*data.X[T.obs.i[1],'X.M'] + p*data.X[T.obs.i[1],'X.F'] + b*data.X[T.obs.i[1],'W.dyad'] + b.a*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'W.dyad'] + b.p*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'W.dyad'] + V[T.obs.i[1]] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c + rho.YF*Y.[t-1] + a*data.X[t,'X.F'] + p*data.X[t,'X.M'] + b*data.X[t,'W.dyad'] + b.a*data.X[t,'X.F']*data.X[t,'W.dyad'] + b.p*data.X[t,'X.M']*data.X[t,'W.dyad'] + V[t] + E[t,'E.F']

Y.M[t] = c + rho.Y*Y.M[t-1] + a*data.X[t,'X.M'] + p*data.X[t,'X.F'] + b*data.X[t,'W.dyad'] + b.a*data.X[t,'X.M']*data.X[t,'W.dyad'] + b.p*data.X[t,'X.F']*data.X[t,'W.dyad'] + V[t] + E[t,'E.M']
}}}}

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
W.dyad = as.matrix(data.X[,'W.dyad']) 

# Create a data frame
X.Actor = rep(0,nrow(data.Model))
X.Partner = rep(0,nrow(data.Model))
W = rep(0,nrow(data.Model))
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
W[N.id[which(data.Model[N.id,]$Gender=='F')]] = W.dyad[chunk[[i]]]
W[N.id[which(data.Model[N.id,]$Gender=='M')]] = W.dyad[chunk[[i]]]
}

# Create a data frame
data.Model = data.frame(cbind(data.Model,Y,X.Actor,X.Partner,W)) 

# Create lag variable
Y.lag = rep(0,nrow(data.Model))
n.subject = unique(data.Model$subject.ID)
for (j in n.subject){
Y.lag[which(data.Model$subject.ID==j)] = shift(data.Model$Y[which(data.Model$subject.ID==j)])
}

data.Model = cbind(data.Model,Y.lag=Y.lag)

return(data=data.Model)
}