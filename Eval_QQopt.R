# Evaluate the optimal design criterion for given design
# using Q(X|eta) in the BayesainQQ D-opt Design Criterion, for mixed 2 and 3-level quanlitative and quantitative factors. 
# INPUT
# D: design to be evaluated
# ix: the indices of the design from candidate set, optional  
# p2: number of 2-level factors
# p3_c: number of 3-level categorical variables (qualitative factors)
# p3_q: number of 3-level quantitative factors
# zeta: the ratio for prior correlation matrix
# rho: the noise-to-signal ratio
# order: the index of the effects that are assumed in the model
# eta: the fixed eta value (a vector)
# return.weight: if the three weight w0, w1, and w2 will be returned for all candidate points and the selected design points

# Output: 
# crit: design criterion

Eval_QQopt<-function(D, n, ix=NULL, rep.no=NULL, approximate=FALSE, p2, p3_c, p3_q, zeta, rho, order, eta, return.misc=F) {
p <- p2+p3_c+p3_q
if (p==1) D <- as.matrix(D,ncol=1)
if (!approximate & missing(n)) n <- nrow(D)
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
if (p==1) candidate<-as.matrix(candidate,ncol=1)
colnames(candidate) <- colnames(candidate, do.NULL=F, prefix='X')
eff_ix <- eff_order<=order
Fm<- tempF[,eff_ix]
R <- 1/(a1^p2*a2^p3_c*a3^p3_q)*tempR[eff_ix,eff_ix]
invR <- solve(R)
link <- Fm%*%eta
w1 <- exp(link)/(1+exp(link))
w2 <- 1-w1
w0 <- w1*w2

if (!approximate) {
   if (sum(eff_ix)>n) {
      stop('Not enough run size for the model!')
      return(NULL)	
   }
   if (is.null(D) && is.null(ix)) {
   	  stop('Not sufficient design information!')
      return(NULL)	
   }   
   if (is.null(ix)) {
      ix <- numeric(n)
      for (i in 1:n) {
          for (j in 1:N) {
              if (sum(D[i,]==candidate[j,])==p) {
                 ix[i] <- j
                 break
              }
          }
       }
       ix <- sort(ix) 
       ix.unique<-unique(ix)
       if ((length(ix[ix>0]) < n)) {
       	  if (is.null(rep.no)) {
       	  	 stop('Missing the replications for each design point!')
             return(NULL)	
       	  }
       	  else n_i <- rep.no
       }
       else n_i<-table(ix)                      
   }
   else {
   	  ix <- sort(ix)
   	  if (length(ix) < n ) {
   	      if (is.null(rep.no)) {
       	  	 stop('Missing the replications for each design point!')
             return(NULL)	
       	  }
       	  else n_i <- rep.no
      }
      else n_i <- table(ix)
      ix.unique <- unique(ix)
   }
   
   Fc <- Fm[ix.unique,]   
   K0 <- t(Fc)%*%diag(w0[ix.unique]*n_i)%*%Fc
   K1 <- t(Fc)%*%diag(w1[ix.unique]*n_i)%*%Fc+rho*invR
   K2 <- t(Fc)%*%diag(w2[ix.unique]*n_i)%*%Fc+rho*invR
   }
else {
   K0<-n*t(Fm)%*%diag(drop(w0)*drop(D))%*%Fm
   K1 <- n*t(Fm)%*%diag(drop(w1)*drop(D))%*%Fm+rho*invR
   K2 <- n*t(Fm)%*%diag(drop(w2)*drop(D))%*%Fm+rho*invR
}
crit<-log(abs(det(K0)))+1/2*(log(det(K1))+log(det(K2)))
if (return.misc) {
    W.all <- cbind(w0, w1, w2)
    colnames(W.all) <- c('w0','w1','w2')
    if (!approximate) {
       design <- cbind(n_i, candidate[ix.unique,])
       return(list(design=design, index=ix.unique, crit=crit, candidate=candidate, model_all=Fm, model_doe=Fm[ix,], W.all=W.all, W.sel=W.all[ix.unique,]))
    }
    else
       return(list(design=D, crit=crit,candidate=candidate, model_all=Fm, W.all=W.all))
}
else
    return(crit)
}