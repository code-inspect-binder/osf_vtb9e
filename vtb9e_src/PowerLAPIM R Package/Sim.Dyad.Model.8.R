Sim.Dyad.Model.8 = function(N.dyad,T.obs,  
c,a,p,d,d.a,d.p,
sigma.eps.F,sigma.eps.M,rho.eps.FM,
sigma.nu,
mu.XF,sigma.XF,mu.XM,sigma.XM,rho.X,
prob.D,is.center.X){

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
V.j = rnorm(N.dyad,0,sigma.nu)

V = NULL
for (j in 1:N.dyad){
V = c(V,rep(V.j[j],T.obs))
}

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
D.dyad = rbinom(N.dyad*T.obs,1,prob.D)
data.X = cbind(data.X,D.dyad)

# If is.center.X is equal to TRUE, person-mean centered the predictors
if(is.center.X==TRUE){
data.X <- data.X %>% 
group_by(ID) %>% 
mutate(X.F.c = X.F - mean(X.F),
       X.M.c = X.M - mean(X.M))
# Simulate Dependent Variables
Y.F = c + a*data.X[,'X.F.c'] + p*data.X[,'X.M.c'] + d*data.X[,'D.dyad'] + d.a*data.X[,'X.F.c']*data.X[,'D.dyad'] + d.p*data.X[,'X.M.c']*data.X[,'D.dyad'] + V + E[,'E.F']

Y.M = c + a*data.X[,'X.M.c'] + p*data.X[,'X.F.c'] + d*data.X[,'D.dyad'] + d.a*data.X[,'X.M.c']*data.X[,'D.dyad'] + d.p*data.X[,'X.F.c']*data.X[,'D.dyad'] + V + E[,'E.M']
}

if(is.center.X==FALSE){
# Simulate Dependent Variables
Y.F = c + a*data.X[,'X.F'] + p*data.X[,'X.M'] + d*data.X[,'D.dyad'] + d.a*data.X[,'X.F']*data.X[,'D.dyad'] + d.p*data.X[,'X.M']*data.X[,'D.dyad'] + V + E[,'E.F']

Y.M = c + a*data.X[,'X.M'] + p*data.X[,'X.F'] + d*data.X[,'D.dyad'] + d.a*data.X[,'X.M']*data.X[,'D.dyad'] + d.p*data.X[,'X.F']*data.X[,'D.dyad'] + V + E[,'E.M']
}

# Create a data frame
X.Actor = rep(0,nrow(data.Model))
X.Partner = rep(0,nrow(data.Model))
D = rep(0,nrow(data.Model))
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
D[N.id[which(data.Model[N.id,]$Gender=='F')]] = D.dyad[chunk[[i]]]
D[N.id[which(data.Model[N.id,]$Gender=='M')]] = D.dyad[chunk[[i]]]
}

# Create a data frame
data.Model = data.frame(cbind(data.Model,Y,X.Actor,X.Partner,D)) 

return(data=data.Model)
}