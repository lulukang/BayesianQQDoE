.libPaths('~/data/Rlib')
library(AlgDesign)
library(lhs)
#setwd('~/Documents/Dropbox/BayesianQQ/Code')
setwd('~/data/Working/BayesianQQ')
source('DOE_QQopt.R')
source('Eval_QQopt.R')

#-------------------------Generate Global Baysian QQ Design---------------------------------------------------------------- 
B <- 500
q <- 21

order_0<-1
order_1st_s<-c(2,4)
order_1st_n<-c(7,11,16)
order_2nd<-setdiff(1:21, c(order_0,order_1st_s,order_1st_n))
upper <- numeric(q)
lower <- numeric(q)
upper[order_0] <- 6
lower[order_0] <- 0
upper[order_1st_s] <- c(5,5)
lower[order_1st_s] <- c(1,1)
upper[order_1st_n] <- rep(1,3)
lower[order_1st_n] <- rep(-1,3)
upper[order_2nd] <- rep(0.3,21-6)
lower[order_2nd] <- rep(-0.3,21-6)

ETA_D<-maximinLHS(B,q)
minxx<-apply(ETA_D,2,min)
maxxx<-apply(ETA_D,2,max)
ETA_D<-scale(ETA_D,minxx,maxxx-minxx)
ETA_V<-scale(ETA_D,F,1/(upper-lower))
ETA_V<-scale(ETA_V,-lower,F)  

write.table(ETA_V, 'etch_eta.txt',sep=',',col.names=F,row.names=F,quote=F)

ETA_V<-read.table('etch_eta.txt',sep=',',header=F)
ETA_V<-as.matrix(ETA_V)

#BB <- 30
n <- 21*6
p2 <- 0
p3_c <- 0
p3_q <- 5
zeta <- 1/2
order <- 2

rho <- 0.5
Design<-matrix(0,nrow=3^5, ncol=B)
opt.values<-numeric(B)
eta_notwork<-numeric(0)

for (j in 144:B) {
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
    Design[DOE$index,j]<-DOE$design[,1]
    opt.values[j] <- DOE$opt
}
write.table(Design,'etch_global_rep_rho05.txt',sep=',',col.names=F,row.names=F,quote=F)
# write.table(eta_notwork, 'etch_bad_eta_rho05.txt',sep=',',col.names=F,row.names=F,quote=F)
# not bad eta values
write.table(opt.values,'etch_global_opt_rho05.txt',sep=',',col.names=F,row.names=F,quote=F)

candi<-read.table('etch_candi.txt',sep=',',header=F)
model_candi<-read.table('etch_model.txt',sep=',',header=F)
model_candi<-as.matrix(model_candi)
Design<-matrix(0,nrow=3^5, ncol=B)
opt.values<-numeric(B)
for (j in 1:B) {
    eta<-ETA_V[j,]
    link<-model_candi%*%eta
    w1<-exp(link)/(1+exp(link))
    w0<-w1*(1-w1)
    model_eta<-diag(drop(sqrt(w0)))%*%model_candi
    DOE.naive1<- optFederov(~.-1, data=model_eta, nTrials=21*4,approximate=TRUE)
    DOE.naive2<- optFederov(~.-1, data=model_candi, nTrials=21*2, approximate=TRUE)
    naive.ix <-sort(c(rep(DOE.naive1$rows, times=DOE.naive1$design[,1]), rep(DOE.naive2$rows, times=DOE.naive2$design[,1])))
    DOE.naive<-list(ix=unique(naive.ix),rep.no=table(naive.ix))
    Design[DOE.naive$ix, j]<-DOE.naive$rep.no
    opt.values[j]<-Eval_QQopt(D=NULL, n=21*6, ix=DOE.naive$ix, rep.no=DOE.naive$rep.no, approximate=F, p2=p2, p3_c=p3_c, p3_q=p3_q, zeta=zeta, rho=0.5, order=2, eta=eta)
}
write.table(Design,'etch_global_rep_naive.txt',sep=',',col.names=F,row.names=F,quote=F)
write.table(opt.values,'etch_global_opt_naive.txt',sep=',',col.names=F,row.names=F,quote=F)

#---------------------Compare the designs--------------------------------------------------------

Design.rho05<-read.table('etch_global_rep_rho05.txt',sep=',',header=F)
Design.naive<-read.table('etch_global_rep_naive.txt',sep=',',header=F)

opt.rho05<-read.table('etch_global_opt_rho05.txt',sep=',',header=F)$V1
opt.naive<-read.table('etch_global_opt_naive.txt',sep=',',header=F)$V1
local.eff <- exp(1/21*(opt.rho05-opt.naive))

d.rho05 <- rowSums(Design.rho05)
d.rho05 <- d.rho05/sum(d.rho05)
d.naive <- rowSums(Design.naive)
d.naive <- d.naive/sum(d.naive)

par(mfrow=c(1,2),cex=1.1)
bar_color<-rep('grey',72)
bar_color[sort(d.rho05,decreasing=T,index.return=T)$ix[1:21]]<-'blue'
barplot(d.rho05, main=expression(paste('(a) ',rho==0.5)), col=bar_color, xlab="Candidate Design Points", ylab='Frequency')

bar_color<-rep('grey',72)
bar_color[sort(d.naive,decreasing=T,index.return=T)$ix[1:21]]<-'blue'
barplot(d.naive, main='(b) Global Combiend Design', col=bar_color, xlab="Candidate Design Points", ylab='Frequency')

B <- 100
q <- 21
order_0<-1
order_1st_s<-c(2,4)
order_1st_n<-c(7,11,16)
order_2nd<-setdiff(1:21, c(order_0,order_1st_s,order_1st_n))
upper <- numeric(q)
lower <- numeric(q)
upper[order_0] <- 6
lower[order_0] <- 0
upper[order_1st_s] <- c(5,5)
lower[order_1st_s] <- c(1,1)
upper[order_1st_n] <- rep(1,3)
lower[order_1st_n] <- rep(-1,3)
upper[order_2nd] <- rep(0.3,21-6)
lower[order_2nd] <- rep(-0.3,21-6)

ETA_D<-maximinLHS(B,q)
minxx<-apply(ETA_D,2,min)
maxxx<-apply(ETA_D,2,max)
ETA_D<-scale(ETA_D,minxx,maxxx-minxx)
ETA_V<-scale(ETA_D,F,1/(upper-lower))
ETA_V<-scale(ETA_V,-lower,F)  

write.table(ETA_V, 'etch_eta2.txt',sep=',',col.names=F,row.names=F,quote=F)

ETA_V<-read.table('etch_eta2.txt',sep=',',header=F)
ETA_V<-as.matrix(ETA_V)

ffd<-read.table('etch_FFD.txt',sep='\t',header=F)
ffd<-ffd-1
ffd<-as.matrix(ffd)
ffd<-ffd[rep(1:27,each=5),]

crit_value<-matrix(0,nrow=B,ncol=3)
for (j in 1:B) {
    crit_value[j,1] <- Eval_QQopt(D=d.rho05, n=21*6, approximate=TRUE, p2=0, p3_c=0, p3_q=5, zeta=1/2, rho=0.5, order=2, eta=drop(ETA_V[j,]))
    crit_value[j,2] <- Eval_QQopt(D=d.naive, n=21*6, approximate=TRUE, p2=0, p3_c=0, p3_q=5, zeta=1/2, rho=0.5, order=2, eta=drop(ETA_V[j,]))
    crit_value[j,3] <- Eval_QQopt(D=ffd, n=135, approximate=FALSE, p2=0, p3_c=0, p3_q=5, zeta=1/2, rho=0.5, order=2, eta=drop(ETA_V[j,]))
}

global.eff1 <- exp(1/21*(crit_value[,1]-crit_value[,2]))
global.eff2 <- exp(1/21*(crit_value[,1]-crit_value[,3]))

par(mfrow=c(1,3),cex=1.5)
hist(local.eff, breaks=10, xlim=c(1, 1.20), xlab='efficiencies',main='(a)')
hist(global.eff1, breaks=10, xlim=c(0.95, 1.20), xlab='efficiencies',main='(b)')
hist(global.eff1, breaks=10, xlim=c(0.95, 1.20), xlab='efficiencies',main='(c)')


