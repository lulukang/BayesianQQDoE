.libPaths('~/data/Rlib')
library(AlgDesign)
library(lhs)
#setwd('~/Documents/Dropbox/BayesianQQ/Code')
setwd('~/data/Working/BayesianQQ')
source('DOE_QQopt.R')
source('Eval_QQopt.R')

#-----------------------------Fixed Eta QQ DOE---------------------------------------------------------
n <- 66
p2 <- 3
p3_c <- 1
p3_q <- 1
zeta <- 1/2

order <- 2
#eta <- numeric(22)
#order_1st<-c(1,2,3,5,8,16)
#order_2nd<-setdiff(1:22, order_1st)
#eta[order_1st] <- runif(6,min=-1,max=1)
#eta[order_2nd] <- runif(22-6,min=-0.3, max=0.3)
maxit <- 500
eps <-1e-4

#write.table(eta,'arti_eta_fix.txt',sep=',',col.names=F, row.names=F,quote=F)
eta<-read.table('arti_eta_fix.txt',sep=',',header=F)
eta<-eta$V1

rho <- 0
B<-10
opt <- -Inf
for (j in 1:B) {
     DOE<-DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, ka=0.8, return.misc=T)
     if (DOE$opt > opt) {
     	DOE0 <- DOE
     	opt <- DOE$opt
     }
}
# the opt value is 133.9972

rho <- 0.3
B<-30
opt <- -Inf
for (j in 1:B) {
     DOE<-DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, ka=0.8)
     if (DOE$opt > opt) {
     	DOE1 <- DOE
     	opt <- DOE$opt
     }
}
# the opt value is 135.553

# compare with naive design 
# 1. GLM design for n=66
candi<-read.table('arti_candi.txt',sep=',',header=F)
model_candi<-read.table('arti_model.txt',sep=',',header=F)
model_candi<-as.matrix(model_candi)
link<-model_candi%*%eta
w1<-exp(link)/(1+exp(link))
w0<-w1*(1-w1)
model_eta<-diag(drop(sqrt(w0)))%*%model_candi
DOE.glm <- optFederov(~.-1, data=model_eta, nTrials=n,approximate=TRUE)
opt.glm1 <- Eval_QQopt(D=NULL, n=66, ix=DOE.glm$rows, rep.no=DOE.glm$design[,1], approximate=F, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0, order=2, eta=eta)
# opt=132.4668
opt.glm2 <- Eval_QQopt(D=NULL, n=66, ix=DOE.glm$rows, rep.no=DOE.glm$design[,1], approximate=F, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0.3, order=2, eta=eta)
# opt=134.0737
# 2. LM design for n=66
DOE.lm <- optFederov(~.-1, data=model_candi, nTrials=n, approximate=TRUE)
opt.lm1 <- Eval_QQopt(D=NULL, n=66, ix=DOE.lm$rows, rep.no=DOE.lm$design[,1], approximate=F, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0, order=2, eta=eta)
# opt=131.6846
opt.lm2 <- Eval_QQopt(D=NULL, n=66, ix=DOE.lm$rows, rep.no=DOE.lm$design[,1], approximate=F, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0.3, order=2, eta=eta)
# opt=133.2188

# 3. GLM (44) + LM (22) design
DOE.naive1<- optFederov(~.-1, data=model_eta, nTrials=44,approximate=TRUE)
DOE.naive2<- optFederov(~.-1, data=model_candi, nTrials=22, approximate=TRUE)
naive.ix <-sort(c(rep(DOE.naive1$rows, times=DOE.naive1$design[,1]), rep(DOE.naive2$rows, times=DOE.naive2$design[,1])))
DOE.naive<-list(ix=unique(naive.ix),rep.no=table(naive.ix))
opt.naive1<-Eval_QQopt(D=NULL, n=66, ix=DOE.naive$ix, rep.no=DOE.naive$rep.no, approximate=F, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0, order=2, eta=eta)
# opt=132.9859
opt.naive2<-Eval_QQopt(D=NULL, n=66, ix=DOE.naive$ix, rep.no=DOE.naive$rep.no, approximate=F, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0.3, order=2, eta=eta)
# opt=134.6432

rep.rho0 <- numeric(72)
rep.rho0[DOE0$index] <- DOE0$design[,1]
rep.rho03<- numeric(72)
rep.rho03[DOE1$index] <- DOE1$design[,1]
rep.glm <- numeric(72)
rep.glm[DOE.glm$rows] <- DOE.glm$design[,1]
rep.lm <- numeric(72)
rep.lm[DOE.lm$rows] <- DOE.lm$design[,1]
rep.naive <- numeric(72)
rep.naive[DOE.naive$ix] <-DOE.naive$rep.no

write.table(cbind(candi,rep.rho0,rep.rho03,rep.glm,rep.lm,rep.naive),'arti_fix_does.txt', col.names=F,row.names=F,sep='\t',quote=F)

designs<-read.table('arti_fix_does.txt', header=F)

write.table(read.table('arti_fix_does.txt', header=F),'arti_fix_does2.txt',sep='&',eol='\\\n', quote=F,col.names=F,row.names=T)
