q<-6
pi_1<-0.2
pi_N<-0.8

n<-12:100
bound<-numeric(length(n))

for (i in n) {
	bound[i-11]<-(pi_1/0.5)^i*(pnorm((i-q-0.5*i)/(0.5*sqrt(i)))-pnorm((q-1-0.5*i)/(0.5*sqrt(i))))
}

