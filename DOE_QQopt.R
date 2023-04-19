# Exchange Algorithm to find the exact optimal design (with replications)
# This algorithm is for factorial design, meaning that the candidate points are full factorial design points, and the model contains the effects for the factorial design models, and are also generated in the function. 
# using Q(X|eta) in the BayesainQQ D-opt Design Criterion, for mixed 2 and 3-level quanlitative and quantitative factors. 
# First: random sampling with probability inversely proportional with the deletion function is used when select potential point for exchange. 
# Then: cyclically search is performaed to further improve the design
# INPUT
# n: number of experimental run size, including replications. 
# p2: number of 2-level factors
# p3_c: number of 3-level categorical variables (qualitative factors)
# p3_q: number of 3-level quantitative factors
# zeta: the ratio for prior correlation matrix
# rho: the noise-to-signal ratio
# order: the index of the effects that are assumed in the model
# eta: the fixed eta value (a vector)
# intD: initial design. It is not provided by default. intD must be the indices of the design points from the candidate set. 
# maxit: maximum number of iterations allowed
# eps: convergence criteiron on the standard devation of the last 10 opt value
# return.weight: if the three weight w0, w1, and w2 will be returned for all candidate points and the selected design points

# Output
# index: index of the points selected in the design
# design: returned design
# opt: the optimal design criterion value

DOE_QQopt<-function(n, p2, p3_c, p3_q, zeta, rho, order, eta, intD=NULL, ka=0.85, maxit=500, eps=1e-4,return.misc=F) {

p <- p2+p3_c+p3_q
N <- 2^p2*3^(p3_c+p3_q)

# F is the full model matrix for all the candidate points;
F2 <- matrix(c(1,-1,1,1),nrow=2,ncol=2,byrow=T)
F3 <- matrix(c(1, -(1.5)^0.5, 0.5^0.5, 1, 0, -2^0.5, 1, 1.5^0.5, 0.5^0.5), nrow=3, ncol=3, byrow=T)
a1 <- (1+zeta)/2;
a2 <- 1/3*(1+2*zeta);
a3 <- 1/9*(3+4*zeta+2*zeta^4);
R2 <- diag(c(a1,(1-zeta)/2))
R3_c <- diag(c(a2, (1-zeta)/3, (1-zeta)/3))
R3_q <- diag(c(a3, 1/3*(1-zeta^4), 1/9*(3-4*zeta+zeta^4)))
R3_q[1,3] <- 1/9*sqrt(2)*zeta*(zeta^3-1)
R3_q[3,1] <- R3_q[1,3]
one2 <- matrix(c(1,1),nrow=2,ncol=1,byrow=T)
one3 <- matrix(c(1,1,1),nrow=3,ncol=1,byrow=T)
D2 <- matrix(c(-1,1), nrow=2, ncol=1, byrow=T)
D3 <- matrix(c(-1,0,1), nrow=3, ncol=1, byrow=T)
eff_order <- 0
tempF <- 1
tempR <- 1
candidate <- matrix(1,nrow=1,ncol=1)
if (p2 > 0) 
for (i in 1:p2) {
    tempF <- kronecker(F2,tempF)
    tempR <- kronecker(R2,tempR)
    candidate <- cbind(kronecker(one2, candidate), kronecker(D2,matrix(1,nrow=nrow(candidate), ncol=1)))
    eff_order <- kronecker(matrix(c(0, 1), nrow=1), matrix(1,nrow=1, ncol= 2^(i-1)))+kronecker(t(one2),eff_order)
}
if (p3_c >0)
for (i in 1:p3_c) {
    tempF <- kronecker(F3,tempF)
    tempR <- kronecker(R3_c,tempR)
    candidate <- cbind(kronecker(one3,candidate), kronecker(D3,matrix(1,nrow=nrow(candidate), ncol=1)))
    eff_order <- kronecker(matrix(c(0,1,1), nrow=1),matrix(1,nrow=1, ncol=2^p2*3^(i-1)))+kronecker(t(one3),eff_order)
}
if (p3_q > 0)
for (i in 1:p3_q) {
    tempF <- kronecker(F3,tempF)
    tempR <- kronecker(R3_q,tempR)
    candidate <- cbind(kronecker(one3,candidate), kronecker(D3,matrix(1,nrow=nrow(candidate), ncol=1)))
    eff_order <- kronecker(matrix(c(0,1,2), nrow=1),matrix(1,nrow=1, ncol=2^p2*3^(p3_c+i-1)))+kronecker(t(one3),eff_order)
}
candidate <- candidate[,-1]
colnames(candidate) <- colnames(candidate, do.NULL=F, prefix='X')
eff_ix <- eff_order<=order
no_eff<-sum(eff_ix)
if (no_eff>n) {
    stop('Not enough run size for the model!')
    return(NULL)
}
Fm<- tempF[,eff_ix]
R <- 1/(a1^p2*a2^p3_c*a3^p3_q)*tempR[eff_ix,eff_ix]
invR <- solve(R)
link <- Fm%*%eta

if (N<no_eff) {
    stop('Value of eta not appropriate. Too few useful candidate points.')
    return(NULL)
}
w1 <- exp(link)/(1+exp(link))
w2 <- 1-w1
w0 <- w1*w2

# construct the initial design by deleting the "worst" design point one-by-one to the n, then augment the design to size n
if (is.null(intD)) {
ix<-1:N
size_ix<-N
Fc <- Fm[ix,]
M0 <- try(solve(t(Fc)%*%diag(w0[ix])%*%Fc), silent=TRUE)
M1 <- try(solve(t(Fc)%*%diag(w1[ix])%*%Fc+rho*invR), silent=TRUE)
M2 <- try(solve(t(Fc)%*%diag(w2[ix])%*%Fc+rho*invR), silent=TRUE)
if (any(c(class(M0),class(M1),class(M2))=='try-error')) {
    stop('Eta value not appropriate. Singular information matrix!')
    return(NULL)
} 
while(size_ix>no_eff) {
	v0 <- numeric(size_ix)
	v1 <- numeric(size_ix)
	v2 <- numeric(size_ix)
	for (i in 1:size_ix) {
		v0[i] <- Fc[i,] %*% M0 %*% Fc[i,]
		v1[i] <- Fc[i,] %*% M1 %*% Fc[i,]
		v2[i] <- Fc[i,] %*% M2 %*% Fc[i,]
	}
	suppressWarnings(del<- -(log(1-w0[ix]*v0)+1/2*log(1-w1[ix]*v1)+1/2*log(1-w2[ix]*v2)))
	if (any(is.nan(del))) del[is.nan(del)] <- Inf
	kick_ix<-which(del==min(del))
	if (length(kick_ix)>1) kick_ix<-sample(kick_ix,size=1)
	kick<-ix[kick_ix]
    ix <- ix[-kick_ix]
    size_ix<-size_ix-1
    Fc <- Fm[ix,]
    # update the inverse
    temp1<-Fm[kick,]%*%M0
    M0 <- M0+(w0[kick]/(1-w0[kick]*v0[kick_ix]))*t(temp1)%*%temp1
    temp1<-Fm[kick,]%*%M1
    M1 <- M1+(w1[kick]/(1-w1[kick]*v1[kick_ix]))*t(temp1)%*%temp1
    temp1<-Fm[kick,]%*%M2
    M2 <- M2+(w2[kick]/(1-w2[kick]*v2[kick_ix]))*t(temp1)%*%temp1
}
prop<-numeric(size_ix)
ka<-log(1-ka)
for (i in 1:size_ix) {
	prop[i]<-1+ceiling(ka/(log(max(w1[ix[i]],w2[ix[i]]))))
}
prop<-prop/sum(prop)
addix<-sample(ix,size=n-no_eff,replace=T, prob=prop)

ix<-sort(c(ix,addix))
}
else ix <- intD

Fc <- Fm[ix,]
M0 <- solve(t(Fc)%*%diag(w0[ix])%*%Fc)
M1 <- solve(t(Fc)%*%diag(w1[ix])%*%Fc+rho*invR)
M2 <- solve(t(Fc)%*%diag(w2[ix])%*%Fc+rho*invR)
opt <- -log(det(M0))-1/2*(log(det(M1))+log(det(M2)))

# start the iteration
convergence <- F
noit <- 1
updated <- TRUE
opt.track<-numeric(maxit)
opt.track[noit]<-opt
while (noit <= maxit) {
   # find out which one should be kicked out from current design by random sampling. 
   if (updated) {
      v0 <- numeric(N)
      v1 <- numeric(N)
      v2 <- numeric(N)
      for (i in 1:N) {
         v0[i] <- Fm[i,] %*% M0 %*% Fm[i,]
         v1[i] <- Fm[i,] %*% M1 %*% Fm[i,]
         v2[i] <- Fm[i,] %*% M2 %*% Fm[i,]
      }
      suppressWarnings(del <- -(log(1-w0[ix]*v0[ix])+1/2*log(1-w1[ix]*v1[ix])+1/2*log(1-w2[ix]*v2[ix])))
      if (any(is.nan(del))) del[is.nan(del)] <- Inf
   }
   kick <- sample(ix,size=1,prob=1/del)
   kick_ix <- which(ix==kick)
   if (length(kick_ix) > 1) kick_ix <- sample(kick_ix, size=1)

   # find out the best point to pick for exchange. 
   delta <- numeric(N)
   delta0 <- numeric(N)
   delta1 <- numeric(N)
   delta2 <- numeric(N)
   others <- setdiff(1:N, kick)
   for (j in others) {
   	   delta0[j] <- (1+w0[j]*v0[j])*(1-w0[kick]*v0[kick])+w0[kick]*w0[j]*(Fm[kick,]%*% M0 %*% Fm[j,])^2
   	   delta1[j] <- (1+w1[j]*v1[j])*(1-w1[kick]*v1[kick])+w1[kick]*w1[j]*(Fm[kick,]%*% M1 %*% Fm[j,])^2
   	   delta2[j] <- (1+w2[j]*v2[j])*(1-w2[kick]*v2[kick])+w2[kick]*w2[j]*(Fm[kick,]%*% M2 %*% Fm[j,])^2
   }
   delta[others] <- log(delta0[others]) +1/2*(log(delta1[others])+log(delta2[others]))
   delta_opt <- max(delta[others])
   pick <- others[delta[others]==delta_opt]
   if (length(pick)>1) pick <- sample(pick,size=1)
   # update the design and formula.
   if (delta_opt>0) {
   	   opt <- opt+delta_opt
       ix[kick_ix] <- pick
   	   temp0 <- M0%*%Fm[pick,]
   	   temp01 <- M0%*%Fm[kick,]
   	   TempS0<-temp0%*%t(temp01)
   	   M0 <- M0 - (1/delta0[pick])*(w0[pick]*(1-w0[kick]*v0[kick])*temp0%*%t(temp0)-w0[kick]*(1+w0[pick]*v0[pick])*temp01%*%t(temp01)+w0[pick]*w0[kick]*drop(Fm[kick,]%*%M0%*%Fm[pick,])*(TempS0+t(TempS0)))
   	   temp1 <- M1%*%Fm[pick,]
   	   temp11 <- M1%*%Fm[kick,]
   	   TempS1<-temp1%*%t(temp11)
   	   M1 <- M1 - (1/delta1[pick])*(w1[pick]*(1-w1[kick]*v1[kick])*temp1%*%t(temp1)-w1[kick]*(1+w1[pick]*v1[pick])*temp11%*%t(temp11)+w1[pick]*w1[kick]*drop(Fm[kick,]%*%M1%*%Fm[pick,])*(TempS1+t(TempS1)))
   	   temp2 <- M2%*%Fm[pick,]
   	   temp21 <- M2%*%Fm[kick,]
   	   TempS2<-temp2%*%t(temp21)
   	   M2 <- M2 - (1/delta2[pick])*(w2[pick]*(1-w2[kick]*v2[kick])*temp2%*%t(temp2)-w2[kick]*(1+w2[pick]*v2[pick])*temp21%*%t(temp21)+w2[pick]*w2[kick]*drop(Fm[kick,]%*%M2%*%Fm[pick,])*(TempS2+t(TempS2)))
   	   updated <- TRUE
   }
   else updated <- FALSE
   noit <- noit+1
   opt.track[noit]<-opt
   if (noit>10) {
   	  if (sd(opt.track[(noit-10):noit])<eps) {
         convergence <- T
         break
      }
   }
}
##----cyclically check every point to see if further improvement is possible-------------
ix <- sort(ix)
any.change <- T
updated <- T
circle.no <-0
while (any.change && (noit+circle.no < maxit)) {
	any.change <- F
	for (kick_ix in 1:n) {
	   # cyclically update the design points
       if (updated) {
          v0 <- numeric(N)
          v1 <- numeric(N)
          v2 <- numeric(N)
          for (i in 1:N) {
            v0[i] <- Fm[i,] %*% M0 %*% Fm[i,]
            v1[i] <- Fm[i,] %*% M1 %*% Fm[i,]
            v2[i] <- Fm[i,] %*% M2 %*% Fm[i,]
          }
       }
       kick <- ix[kick_ix]
   
       # find out the best point to pick for exchange. 
       delta <- numeric(N)
       delta0 <- numeric(N)
       delta1 <- numeric(N)
       delta2 <- numeric(N)
       others <- setdiff(1:N, kick)
       for (j in others) {
   	       delta0[j] <- (1+w0[j]*v0[j])*(1-w0[kick]*v0[kick])+w0[kick]*w0[j]*(Fm[kick,]%*% M0 %*% Fm[j,])^2
   	       delta1[j] <- (1+w1[j]*v1[j])*(1-w1[kick]*v1[kick])+w1[kick]*w1[j]*(Fm[kick,]%*% M1 %*% Fm[j,])^2
   	       delta2[j] <- (1+w2[j]*v2[j])*(1-w2[kick]*v2[kick])+w2[kick]*w2[j]*(Fm[kick,]%*% M2 %*% Fm[j,])^2
       }
       delta[others] <- log(delta0[others]) +1/2*(log(delta1[others])+log(delta2[others]))
       delta_opt <- max(delta[others])
       pick <- others[delta[others]==delta_opt]
       if (length(pick)>1) pick <- sample(pick,size=1)
       # update the design and formula.
       if (delta_opt>0) {
   	      opt <- opt+delta_opt
          ix[kick_ix] <- pick
   	      temp0 <- M0%*%Fm[pick,]
   	      temp01 <- M0%*%Fm[kick,]
   	      TempS0<-temp0%*%t(temp01)
   	      M0 <- M0 - (1/delta0[pick])*(w0[pick]*(1-w0[kick]*v0[kick])*temp0%*%t(temp0)-w0[kick]*(1+w0[pick]*v0[pick])*temp01%*%t(temp01)+w0[pick]*w0[kick]*drop(Fm[kick,]%*%M0%*%Fm[pick,])*(TempS0+t(TempS0)))
   	      temp1 <- M1%*%Fm[pick,]
   	      temp11 <- M1%*%Fm[kick,]
   	      TempS1<-temp1%*%t(temp11)
   	      M1 <- M1 - (1/delta1[pick])*(w1[pick]*(1-w1[kick]*v1[kick])*temp1%*%t(temp1)-w1[kick]*(1+w1[pick]*v1[pick])*temp11%*%t(temp11)+w1[pick]*w1[kick]*drop(Fm[kick,]%*%M1%*%Fm[pick,])*(TempS1+t(TempS1)))
   	      temp2 <- M2%*%Fm[pick,]
   	      temp21 <- M2%*%Fm[kick,]
   	      TempS2<-temp2%*%t(temp21)
   	      M2 <- M2 - (1/delta2[pick])*(w2[pick]*(1-w2[kick]*v2[kick])*temp2%*%t(temp2)-w2[kick]*(1+w2[pick]*v2[pick])*temp21%*%t(temp21)+w2[pick]*w2[kick]*drop(Fm[kick,]%*%M2%*%Fm[pick,])*(TempS2+t(TempS2)))
   	      if (!any.change) any.change <- T
   	      updated <- TRUE
        }
        else updated <- FALSE
   }
   circle.no <-circle.no+1
   opt.track[noit+circle.no]<-opt
}

#-----------construct the returning object--------------------------------------------------------
ix <- sort(ix)
unique.ix <- unique(ix)
rep.no <- table(ix)
design <- cbind(rep.no, candidate[unique.ix,])
if (return.misc) {
    W.all <- cbind(w0, w1, w2)
    colnames(W.all) <- c('w0','w1','w2')
    W.doe <- W.all[unique.ix,]
    optDOE<-list(design=design, index=unique.ix, opt=opt, candidate=candidate, model_all=Fm, model_doe=Fm[ix,], W.all=W.all, W.doe=W.doe, convergence=convergence, opt.track=opt.track[1:(noit+circle.no)], noit=noit, circle.no=circle.no)
}
else {
    optDOE<-list(design=design, index=unique.ix, opt=opt, convergence=convergence, circle.no=circle.no)
}
    return(optDOE)
}