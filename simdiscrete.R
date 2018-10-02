##### Simulation using CLIME/gLASSO/neighbourhood selection #########
##### HBIC for parameter selection ##################################
##### chisq test for inference ######################################
#### all discrete case ############################################

library(igraph)
library(mvtnorm)

n <-200
p <-500
N<-p*(p-1)/2
nei <- 2*p/50


####generate graph #################################
set.seed(1000)
sw <- sample_smallworld(dim=1, size=p, nei=2, p=0.05)
adj <- as_adjacency_matrix(sw,sparse = FALSE)
plot(sw)


## generate covariance and precision matrices###
edge <- sum(adj)/2
dg <- adj[upper.tri(adj)]


id <- runif(edge,min = 0.2,max = 0.7)
ome <- adj

ome[upper.tri(ome)& ome==1] <- id
ome[lower.tri(ome)] <- t(ome)[lower.tri(t(ome))]

diag(ome) <-1

eigvals <- eigen(ome, only.values=T)$values
eigvals

perturb <- abs(min(eigvals))+0.01
ome <- ome + diag(p)*perturb

cov <- solve(ome)
cor <- cov2cor(cov)

### generate data########################
mu <- rep(0,p)

z <- rmvnorm(n,mu,cor)
###### L=2 ##########
######All discrete###
## x is the discretized variables###
x <- z
C1 <- runif(p,min = 0.30,max =0.85 )
C1r <- round(C1,digits = 2)
C2 <- runif(p,min = 1.5, max =3)
C2r <- round(C2,digits = 2)
for(i in 1:n){
  for(j in 1:p){
    x[i,j] <- sum(z[i,j]>C1r[j])+sum(z[i,j]>C2r[j])
  }
}


###estimate statistics#####
###R_ij######
R <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    a=0
    for (k in 1:(n-1)){
      for (l in (k+1):n){
        a = 2*(x[k,i]-x[l,i])*(x[k,j]-x[l,j])/(n*(n-1));
        R[i,j] = a + R[i,j];
        R[j,i] = R[i,j]
      }
    }
  }
}

##estimate delta_j###
### the importance of this step is to add an pertubation term#####
delta_1 <- rep(0,p)
delta_2 <- rep(0,p)
for(j in 1:p){
  b = sum(x[,j]==2)/n + 0.001
  delta_2[j] <- qnorm(1-b)
  d = sum(x[,j])/n
  delta_1[j] <- qnorm(1-(d-b))
}

## find the relationship between R and pho##

trans_dis <- function(x,delta1,delta2){
  corr = matrix(c(1,x,x,1),2,2,byrow = T)
  pb <- pmvnorm(lower = c(delta1,delta2), upper = c(Inf,Inf), corr = corr) 
  pb[1]
}

Pho <- matrix(0,p,p)
diag(Pho) <-1
for(i in 1:(p-1)){
  for (j in (i+1):p){
    
    u <- R[i,j]/2 + (pnorm(-delta_1[i])+pnorm(-delta_2[i]))*(pnorm(-delta_1[j])+pnorm(-delta_2[j]))
    fun <- function(x) {
      a <- trans_dis(1,delta_1[i],delta_1[j]) + trans_dis(1,delta_1[i],delta_2[j]) + trans_dis(1,delta_2[i],delta_1[j])+ trans_dis(1,delta_2[i],delta_2[j]) 
      b <- trans_dis(-1,delta_1[i],delta_1[j]) + trans_dis(-1,delta_1[i],delta_2[j]) + trans_dis(-1,delta_2[i],delta_1[j])+ trans_dis(-1,delta_2[i],delta_2[j]) 
      if (u>a) { u <- a}
      else if (u < b) { u <- b}
      
      trans_dis(x,delta_1[i],delta_1[j]) + trans_dis(x,delta_1[i],delta_2[j]) + trans_dis(x,delta_2[i],delta_1[j])+ trans_dis(x,delta_2[i],delta_2[j]) - u
    }
    a<-uniroot(fun,c(-1,1))
    
    Pho[i,j] <- a$root
    Pho[j,i] <- Pho[i,j]
  }
}


sig1 <- cor(x)
## latent transformation####
## add perbutation###
## this is not a good way, it is better to use the projection method to find a semidefinite matrix ###
eigvals <- eigen(Pho, only.values=T)$values
perturb <- max(max(eigvals) - p*min(eigvals), 0)/(p-1)
sig2 <- Pho + diag(p)*perturb

#### CLIME method ##########################
library(flare)
####Pearson-CLIME####
#nlambda <- 50
lambda1 <- seq(0.6,0,-0.02)
pclime <- sugm(sig1, lambda = lambda1,  nlambda = NULL , lambda.min.ratio = NULL,
               method = "clime", sym = "or", shrink=NULL,
               prec = 1e-4, max.ite = 1e4, standardize = FALSE,
               perturb = TRUE, verbose = TRUE)

#### HBIC FOR THE OPTIMAL PARAMETER #############
HBIC_p <- mat.or.vec(length(lambda1),1)
for(h in 1:length(lambda1)){
  icov <- pclime$icov[[h]]
  HBIC_p[h] <- sum(diag(sig1 %*% icov)) - log(det(icov)) + log(log(n)) * log(p) * sum(icov[upper.tri(icov)]!=0)/n
}

# HBIC_p
M <- pclime$icov[[which.min(HBIC_p)]]

####Latent-CLIME####
#nlambda <- 50
lambda2 <- seq(1.1,0,-0.02)

lclime <- sugm(sig2, lambda2, lambda.min.ratio = NULL,
               method = "clime", sym = "or", shrink=NULL,
               prec = 1e-4, max.ite = 1e4, standardize = FALSE,
               perturb = TRUE, verbose = TRUE)

#### HBIC FOR THE OPTIMAL PARAMETER #############
HBIC_l <- mat.or.vec(length(lambda2),1)
for(h in 1:length(lambda2)){
  icov <- lclime$icov[[h]]
  HBIC_l[h] <- sum(diag(Pho %*% icov)) - log(det(icov)) + log(log(n)) * log(p) * sum(icov[upper.tri(icov)]!=0)/n
}

M <- lclime$icov[[which.min(HBIC_l)]]



#####gLasso ###################
library(huge)
pgl <- huge(sig1,lambda = lambda1,method = "glasso")

lgl <- huge(sig2,lambda = lambda2, method = "glasso")

#### HBIC FOR THE OPTIMAL PARAMETER #############
HBIC_lg <- mat.or.vec(length(lambda2),1)
for(h in 1:length(lambda2)){
  icov <- as.matrix(lgl$icov[[h]])
  HBIC_lg[h] <- sum(diag(Pho %*% icov)) - log(det(icov)) + log(log(n)) * log(p) * sum(icov[upper.tri(icov)]!=0)/n
}

M <- as.matrix(lgl$icov[[which.min(HBIC_lg)]])

##### neighbourhood selection ######
lnl <- huge(sig2,lambda = lambda2)

HBIC_ln <- mat.or.vec(length(lambda2),1)
for(h in 1:length(lambda2)){
  icov <- as.matrix(lnl$beta[[h]])
  HBIC_ln[h] <- sum(diag(Pho %*% icov)) - log(det(icov)) + log(log(n)) * log(p) * sum(icov[upper.tri(icov)]!=0)/n
}

M <- as.matrix(lnl$beta[[which.min(HBIC_ln)]])



