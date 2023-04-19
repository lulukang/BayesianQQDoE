.libPaths('~/data/Rlib')
library(AlgDesign)
library(lhs)
#setwd('~/Documents/Dropbox/BayesianQQ/Code')
setwd('~/data/Working/BayesianQQ')
source('DOE_QQopt.R')
source('Eval_QQopt.R')

#-------------------------Generate Global Baysian QQ Design---------------------------------------------------------------- 
B <- 500

order_1st<-c(1,2,3,5,8,16)
order_2nd<-setdiff(1:22, order_1st)
upper <- numeric(22)
lower <- numeric(22)
upper[order_1st] <- rep(2,6)
upper[order_2nd] <- rep(0.5,16)
lower[order_1st] <- rep(-2,6)
lower[order_2nd] <- rep(-0.5,16)

ETA_D<-maximinLHS(B,22)
minxx<-apply(ETA_D,2,min)
maxxx<-apply(ETA_D,2,max)
ETA_D<-scale(ETA_D,minxx,maxxx-minxx)
ETA_V<-scale(ETA_D,F,1/(upper-lower))
ETA_V<-scale(ETA_V,-lower,F)  
write.table(ETA_V, 'arti_eta.txt',sep=',',col.names=F,row.names=F,quote=F)

BB <- 10
n <- 66
p2 <- 3
p3_c <- 1
p3_q <- 1
zeta <- 1/2
order <- 2

rho <- 0
Design<-matrix(0,nrow=72, ncol=B)
opt.values<-numeric(B)
eta_notwork<-numeric(0)

for (j in 1:B) {
    eta<-ETA_V[j,]
    opt<--Inf
    for (i in 1:BB) {
       DOE<-try(DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, ka=0.8), silent=TRUE)
       if (class(DOE)=='try-error') {
          eta_notwork <- c(eta_notwork,j)
          DOE0 <- NULL
          break
       }
       if (DOE$opt>opt) {
       	  DOE0 <- DOE
       	  opt <- DOE$opt
       }
    }
    Design[DOE0$index,j]<-DOE0$design[,1]
    opt.values[j] <- DOE0$opt
}

write.table(Design,'arti_global_rep_rho0.txt',sep=',',col.names=F,row.names=F,quote=F)
# write.table(eta_notwork, 'arti_bad_eta_rho0.txt',sep=',',col.names=F,row.names=F,quote=F)
# no bad eta values for this set of ETA values.
write.table(opt.values,'arti_global_opt_rho0.txt',sep=',',col.names=F,row.names=F,quote=F)


rho <- 0.3
Design<-matrix(0,nrow=72, ncol=B)
opt.values<-numeric(B)
eta_notwork<-numeric(0)

for (j in 1:B) {
    eta<-ETA_V[j,]
    opt<--Inf
    for (i in 1:BB) {
       DOE<-try(DOE_QQopt(n, p2, p3_c, p3_q, zeta, rho, order, eta, ka=0.8), silent=TRUE)
       if (class(DOE)=='try-error') {
          eta_notwork <- c(eta_notwork,j)
          DOE0 <- NULL
          break
       }
       if (DOE$opt>opt) {
       	  DOE0 <- DOE
       	  opt <- DOE$opt
       }
    }
    Design[DOE0$index,j]<-DOE0$design[,1]
    opt.values[j] <- DOE0$opt
}

write.table(Design,'arti_global_rep_rho03.txt',sep=',',col.names=F,row.names=F,quote=F)
# write.table(eta_notwork, 'arti_bad_eta_rho03.txt',sep=',',col.names=F,row.names=F,quote=F)
# no bad eta values for this set of ETA values.
write.table(opt.values,'arti_global_opt_rho03.txt',sep=',',col.names=F,row.names=F,quote=F)

candi<-read.table('arti_candi.txt',sep=',',header=F)
model_candi<-read.table('arti_model.txt',sep=',',header=F)
model_candi<-as.matrix(model_candi)
Design<-matrix(0,nrow=72, ncol=B)
opt.values<-matrix(0, nrow=B,ncol=2)
for (j in 1:B) {
    eta<-ETA_V[j,]
    link<-model_candi%*%eta
    w1<-exp(link)/(1+exp(link))
    w0<-w1*(1-w1)
    model_eta<-diag(drop(sqrt(w0)))%*%model_candi
    DOE.naive1<- optFederov(~.-1, data=model_eta, nTrials=44,approximate=TRUE)
    DOE.naive2<- optFederov(~.-1, data=model_candi, nTrials=22, approximate=TRUE)
    naive.ix <-sort(c(rep(DOE.naive1$rows, times=DOE.naive1$design[,1]), rep(DOE.naive2$rows, times=DOE.naive2$design[,1])))
    DOE.naive<-list(ix=unique(naive.ix),rep.no=table(naive.ix))
    Design[DOE.naive$ix, j]<-DOE.naive$rep.no
    opt.values[j,1]<-Eval_QQopt(D=NULL, n=66, ix=DOE.naive$ix, rep.no=DOE.naive$rep.no, approximate=F, p2=p2, p3_c=p3_c, p3_q=p3_q, zeta=zeta, rho=0, order=2, eta=eta)
    opt.values[j,2]<-Eval_QQopt(D=NULL, n=66, ix=DOE.naive$ix, rep.no=DOE.naive$rep.no, approximate=F, p2=p2, p3_c=p3_c, p3_q=p3_q, zeta=zeta, rho=0.3, order=2, eta=eta)
}
write.table(Design,'arti_global_rep_naive.txt',sep=',',col.names=F,row.names=F,quote=F)
write.table(opt.values,'arti_global_opt_naive.txt',sep=',',col.names=c('rho0','rho03'),row.names=F,quote=F)

#---------------------Compare the designs--------------------------------------------------------

Design.rho0<-read.table('arti_global_rep_rho0.txt',sep=',',header=F)
Design.rho03<-read.table('arti_global_rep_rho03.txt',sep=',',header=F)
Design.naive<-read.table('arti_global_rep_naive.txt',sep=',',header=F)

opt.rho0<-read.table('arti_global_opt_rho0.txt',sep=',',header=F)$V1
opt.rho03<-read.table('arti_global_opt_rho03.txt',sep=',',header=F)$V1
opt.naive<-read.table('arti_global_opt_naive.txt',sep=',',header=T)
local.eff<-matrix(0,nrow=length(opt.rho0),ncol=2)
local.eff[,1] <- exp(1/22*(opt.rho0-opt.naive$rho0))
local.eff[,2] <- exp(1/22*(opt.rho03-opt.naive$rho03))

d.rho0 <- rowSums(Design.rho0)
d.rho0 <- d.rho0/sum(d.rho0)
d.rho03 <- rowSums(Design.rho03)
d.rho03 <- d.rho03/sum(d.rho03)
d.naive <- rowSums(Design.naive)
d.naive <- d.naive/sum(d.naive)

par(mfrow=c(1,3),cex=1.1)
bar_color<-rep('grey',72)
bar_color[sort(d.rho0,decreasing=T,index.return=T)$ix[1:22]]<-'blue'
barplot(d.rho0, main=expression(paste('(a) ',rho==0)), col=bar_color, xlab="Candidate Design Points", ylab='Frequency')

bar_color<-rep('grey',72)
bar_color[sort(d.rho03,decreasing=T,index.return=T)$ix[1:22]]<-'blue'
barplot(d.rho03, main=expression(paste('(c) ',rho==0.3)), col=bar_color, xlab="Candidate Design Points", ylab='Frequency')

bar_color<-rep('grey',72)
bar_color[sort(d.naive,decreasing=T,index.return=T)$ix[1:22]]<-'blue'
barplot(d.naive, main='(c) Global Combined Design', col=bar_color, xlab="Candidate Design Points", ylab='Frequency')

B <- 100
order_1st<-c(1,2,3,5,8,16)
order_2nd<-setdiff(1:22, order_1st)
upper <- numeric(22)
lower <- numeric(22)
upper[order_1st] <- rep(2,6)
upper[order_2nd] <- rep(0.5,16)
lower[order_1st] <- rep(-2,6)
lower[order_2nd] <- rep(-0.5,16)

ETA_D<-maximinLHS(B,22)
minxx<-apply(ETA_D,2,min)
maxxx<-apply(ETA_D,2,max)
ETA_D<-scale(ETA_D,minxx,maxxx-minxx)
ETA_V<-scale(ETA_D,F,1/(upper-lower))
ETA_V<-scale(ETA_V,-lower,F)  
write.table(ETA_V, 'arti_eta2.txt',sep=',',col.names=F,row.names=F,quote=F)

ETA_V<-read.table('arti_eta2.txt',sep=',',header=F)
ETA_V<-as.matrix(ETA_V)
crit_value<-matrix(0,nrow=B,ncol=4)
for (j in 1:B) {
    crit_value[j,1] <- Eval_QQopt(D=d.rho0, n=66, approximate=TRUE, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0, order=2, eta=drop(ETA_V[j,]))
    crit_value[j,2] <- Eval_QQopt(D=d.rho03, n=66, approximate=TRUE, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0.3, order=2, eta=drop(ETA_V[j,]))
    crit_value[j,3] <- Eval_QQopt(D=d.naive, n=66, approximate=TRUE, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0, order=2, eta=drop(ETA_V[j,]))
    crit_value[j,4] <- Eval_QQopt(D=d.naive, n=66, approximate=TRUE, p2=3, p3_c=1, p3_q=1, zeta=1/2, rho=0.3, order=2, eta=drop(ETA_V[j,]))
}

global.eff<-matrix(0,nrow=nrow(crit_value),ncol=2)
global.eff[,1] <- exp(1/22*(crit_value[,1]-crit_value[,3]))
global.eff[,2] <- exp(1/22*(crit_value[,2]-crit_value[,4]))

par(mfrow=c(1,2))
hist(local.eff[,1], breaks=10, xlim=c(1, 1.25), xlab='efficiencies',main='(a)')
hist(local.eff[,2], breaks=10, xlim=c(1, 1.25), xlab='efficiencies',main='(b)')
hist(global.eff[,1], breaks=10, xlim=c(0.95, 1.2), xlab='efficiencies',main='(a)')
hist(global.eff[,2], breaks=10, xlim=c(0.95, 1.2), xlab='efficiencies',main='(b)')

