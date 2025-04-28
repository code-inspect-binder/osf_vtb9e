Sim.Dyad.Model.15.lag = function(N.dyad,T.obs,  
c.F,c.M,rho.YF,rho.YM,a.FF,a.FF2,p.MF,p.MF2,a.MM,a.MM2,p.FM,p.FM2,
d.F,d.M,d.FF,d.FF2,d.MF,d.MF2,d.MM,d.MM2,d.FM,d.FM2,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu.F,sigma.nu.M,rho.nu.F.M,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X){

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
var.diag.X = c(sigma.XF,sigma.XM)
Sigma.X = diag(length(var.diag.X))
Sigma.X[lower.tri(Sigma.X, diag=FALSE)] = rho.X
Sigma.X = pmax(Sigma.X, t(Sigma.X), na.rm=TRUE)
Sigma.X = diag(var.diag.X)%*%Sigma.X%*%diag(var.diag.X)
data.X = cbind(data.X,mvrnorm(N.dyad*T.total,c(mu.XF,mu.XM),Sigma.X))
colnames(data.X) = c('Obs','ID','X.F','X.M')
X.F = data.X[,'X.F']
X.M = data.X[,'X.M']
D.dyad = rbinom(N.dyad*T.total,1,prob.D)
data.X = cbind(data.X,D.dyad)

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
Y.F[T.obs.i[1]] = c.F + a.FF*data.X[T.obs.i[1],'X.F.c'] + a.FF2*I(data.X[T.obs.i[1],'X.F.c']^2) + p.MF*data.X[T.obs.i[1],'X.M.c'] + p.MF2*I(data.X[T.obs.i[1],'X.M.c']^2) + d.F*data.X[T.obs.i[1],'D.dyad'] + d.FF*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'D.dyad'] + d.FF2*I(data.X[T.obs.i[1],'X.F.c']^2)*data.X[T.obs.i[1],'D.dyad'] + d.MF*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'D.dyad'] + d.MF2*I(data.X[T.obs.i[1],'X.M.c']^2)*data.X[T.obs.i[1],'D.dyad'] + V[T.obs.i[1],'V.F'] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c.M + a.MM*data.X[T.obs.i[1],'X.M.c'] + a.MM2*I(data.X[T.obs.i[1],'X.M.c']^2) + p.FM*data.X[T.obs.i[1],'X.F.c'] + p.FM2*I(data.X[T.obs.i[1],'X.F.c']^2) + d.M*data.X[T.obs.i[1],'D.dyad'] + d.MM*data.X[T.obs.i[1],'X.M.c']*data.X[T.obs.i[1],'D.dyad'] + d.MM2*I(data.X[T.obs.i[1],'X.M.c']^2)*data.X[T.obs.i[1],'D.dyad'] + d.FM*data.X[T.obs.i[1],'X.F.c']*data.X[T.obs.i[1],'D.dyad'] + d.FM2*I(data.X[T.obs.i[1],'X.F.c']^2)*data.X[T.obs.i[1],'D.dyad'] + V[T.obs.i[1],'V.M'] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c.F + rho.YF*Y.F[t-1] + a.FF*data.X[t,'X.F.c'] + a.FF2*I(data.X[t,'X.F.c']^2) + p.MF*data.X[t,'X.M.c'] + p.MF2*I(data.X[t,'X.M.c']^2) + d.F*data.X[t,'D.dyad'] + d.FF*data.X[t,'X.F.c']*data.X[t,'D.dyad'] + d.FF2*I(data.X[t,'X.F.c']^2)*data.X[t,'D.dyad'] + d.MF*data.X[t,'X.M.c']*data.X[t,'D.dyad'] + d.MF2*I(data.X[t,'X.M.c']^2)*data.X[t,'D.dyad'] + V[t,'V.F'] + E[t,'E.F']

Y.M[t] = c.M + rho.YM*Y.M[t-1] + a.MM*data.X[t,'X.M.c'] + a.MM2*I(data.X[t,'X.M.c']^2) + p.FM*data.X[t,'X.F.c'] + p.FM2*I(data.X[t,'X.F.c']^2) + d.M*data.X[t,'D.dyad'] + d.MM*data.X[t,'X.M.c']*data.X[t,'D.dyad'] + d.MM2*I(data.X[t,'X.M.c']^2)*data.X[t,'D.dyad'] + d.FM*data.X[t,'X.F.c']*data.X[t,'D.dyad'] + d.FM2*I(data.X[t,'X.F.c']^2)*data.X[t,'D.dyad'] + V[t,'V.M'] + E[t,'E.M']
}}

if (is.center.X == FALSE){
# Initialized values
Y.F[T.obs.i[1]] = c.F + a.FF*data.X[T.obs.i[1],'X.F'] + a.FF2*I(data.X[T.obs.i[1],'X.F']^2) + p.MF*data.X[T.obs.i[1],'X.M'] + p.MF2*I(data.X[T.obs.i[1],'X.M']^2) + d.F*data.X[T.obs.i[1],'D.dyad'] + d.FF*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'D.dyad'] + d.FF2*I(data.X[T.obs.i[1],'X.F']^2)*data.X[T.obs.i[1],'D.dyad'] + d.MF*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'D.dyad'] + d.MF2*I(data.X[T.obs.i[1],'X.M']^2)*data.X[T.obs.i[1],'D.dyad'] + V[T.obs.i[1],'V.F'] + E[T.obs.i[1],'E.F']

Y.M[T.obs.i[1]] = c.M + a.MM*data.X[T.obs.i[1],'X.M'] + a.MM2*I(data.X[T.obs.i[1],'X.M']^2) + p.FM*data.X[T.obs.i[1],'X.F'] + p.FM2*I(data.X[T.obs.i[1],'X.F']^2) + d.M*data.X[T.obs.i[1],'D.dyad'] + d.MM*data.X[T.obs.i[1],'X.M']*data.X[T.obs.i[1],'D.dyad'] + d.MM2*I(data.X[T.obs.i[1],'X.M']^2)*data.X[T.obs.i[1],'D.dyad'] + d.FM*data.X[T.obs.i[1],'X.F']*data.X[T.obs.i[1],'D.dyad'] + d.FM2*I(data.X[T.obs.i[1],'X.F']^2)*data.X[T.obs.i[1],'D.dyad'] + V[T.obs.i[1],'V.M'] + E[T.obs.i[1],'E.M']

for (t in T.obs.i[-1]){
# Simulate Dependent Variables
Y.F[t] = c.F + rho.YF*Y.F[t-1] + a.FF*data.X[t,'X.F'] + a.FF2*I(data.X[t,'X.F']^2) + p.MF*data.X[t,'X.M'] + p.MF2*I(data.X[t,'X.M']^2) + d.F*data.X[t,'D.dyad'] + d.FF*data.X[t,'X.F']*data.X[t,'D.dyad'] + d.FF2*I(data.X[t,'X.F']^2)*data.X[t,'D.dyad'] + d.MF*data.X[t,'X.M']*data.X[t,'D.dyad'] + d.MF2*I(data.X[t,'X.M']^2)*data.X[t,'D.dyad'] + V[t,'V.F'] + E[t,'E.F']

Y.M[t] = c.M + rho.YM*Y.M[t-1] + a.MM*data.X[t,'X.M'] + a.MM2*I(data.X[t,'X.M']^2) + p.FM*data.X[t,'X.F'] + p.FM2*I(data.X[t,'X.F']^2) + d.M*data.X[t,'D.dyad'] + d.MM*data.X[t,'X.M']*data.X[t,'D.dyad'] + d.MM2*I(data.X[t,'X.M']^2)*data.X[t,'D.dyad'] + d.FM*data.X[t,'X.F']*data.X[t,'D.dyad'] + d.FM2*I(data.X[t,'X.F']^2)*data.X[t,'D.dyad'] + V[t,'V.M'] + E[t,'E.M']
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
D.dyad = as.matrix(data.X[,'D.dyad']) 

# Create a data frame
X.Actor = rep(0,nrow(data.Model))
X.Partner = rep(0,nrow(data.Model))
D = rep(0,nrow(data.Model))
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
D[N.id[which(data.Model[N.id,]$Gender=='F')]] = D.dyad[chunk[[i]]]
D[N.id[which(data.Model[N.id,]$Gender=='M')]] = D.dyad[chunk[[i]]]
}

# Create a data frame
data.Model = data.frame(cbind(data.Model,Y,X.Actor,X.Partner,D)) 

# Create lag variable
Y.lag = rep(0,nrow(data.Model))
n.subject = unique(data.Model$subject.ID)
for (j in n.subject){
Y.lag[which(data.Model$subject.ID==j)] = shift(data.Model$Y[which(data.Model$subject.ID==j)])
}

data.Model = cbind(data.Model,Y.lag=Y.lag)

return(data=data.Model)
}