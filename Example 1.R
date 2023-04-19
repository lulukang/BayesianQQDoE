# Example for fixed eta

#.libPaths('~/data/Rlib')
library(AlgDesign)
library(lhs)
setwd('~/Dropbox/BayesianQQ/Code')
#setwd('~/data/Working/Bayesian QQ')
source('DOE_QQopt.R')
source('Eval_QQopt.R')

# One-dim Example 1. 
n <- 14
p2 <- 0
p3_c <- 0
p3_q <- 1
zeta <- 1/2
rho <- 0
order <- 2
eta_mean<-c(1,1,0)
maxit <- 300
eps <- 1e-2

# find the optimal design for fixed eta value
Opt <- -Inf
eta<-eta_mean
for (i in 1:10) {
     DOE<-DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, maxit, eps, return.weight=T)
     if (DOE$opt > Opt) DOE0 <- DOE
}

#----------Make some comparisons------------------------------------
# This design is actually not as good as the optimal design. 
DOE2<-c(-1,-1,-1,0,0,0,1,1,1,1,1,1,1,1)
Eval_QQopt(DOE2, p2=0, p3_c=0, p3_q=1, zeta=1/2,rho=0, order=2, eta=eta, return.weight=T)

# D-optimal design with 14 run from AlgDesign package. 
DOE3<-optFederov(~quad(.),as.matrix(seq(from=-1,to=1,length=20),ncol=1),nTrials=14, crit='D', maxIteration=100) 

# Find the optimal continous design
spaceD <- maximinLHS(100, 3)
eta_range <- c(3, 1.5, 1.5)
eta <- matrix(0, nrow=100, ncol=3)
for (i in 1:100) {
eta[i,] <- spaceD[i,]*eta_range*2+(eta_mean-eta_range)
}

Freq <- numeric(3)
for (l in 1:100){
    Opt <- -Inf
    for (i  in 1:10) {
        DOE <- DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta[l,], maxit, eps, return.weight=F)
        if (DOE$opt > Opt) DOE0 <- DOE
    }
    Freq <- Freq + c(sum(DOE0$index==1), sum(DOE0$index==2), sum(DOE0$index==3))
}
# 467 468 465

#-----------------------

# 
n <- 24
p2 <- 3
p3_c <- 1
p3_q <- 1
zeta <- 1/2

order <- 2
eta <- runif(22)*2
maxit <- 1000
eps <- 1e-3

rho <- 0
Opt <- -Inf
for (i in 1:10) {
     DOE<-DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, maxit, eps, return.weight=T)
     if (DOE$opt > Opt) DOE0 <- DOE
}

rho <-0.2
Opt <- -Inf
for (i in 1:10) {
     DOE<-DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, maxit, eps, return.weight=T)
     if (DOE$opt > Opt) DOE1 <- DOE
}


