library(statcomp)
library(plot3D)

#Example:

x = arima.sim(model=list(ar = 0.3), n = 10^4)
opd = ordinal_pattern_distribution(x = x, ndemb = 3)
?ordinal_pattern_distribution
opd
permutation_entropy(opd)
permutation_entropy
shannon_entropy

#Functions:

permutation_extropy = function(opd,type="SH") {
  prob = opd / sum(opd)
  opd.prob = 1 - prob
  if(type=="GS") H = - sum(sapply(prob, FUN=function(prob) if (prob >= 1.e-30) return(prob^2 - prob) else return(0)))
  if(type=="SH") H = - sum(sapply(opd.prob, FUN=function(prob) if (prob >= 1.e-30) return(prob * log(prob)) else return(0)))
  return(H)
}

permutation_extropy(opd)

permutation_MEDdivergence = function(opd1,opd2,omega=0.5,type="GS") {
  Pbi <- function(i) 1 - sum(opd1[1:(i-1)])
  Qbi <- function(i) 1 - sum(opd2[1:(i-1)])
  n = length(opd1)

  if(type=="GS"){
  	phi <- function(x) x^2 #-x
	A = B = 0
	for(i in 1:n){
		A = A + phi((1 - Pbi(i))/(1 - Qbi(i)))*(1 - Qbi(i))
		B = B * phi(Pbi(i)/Qbi(i))*Qbi(i)
	}
  }

  if(type=="SH"){
  	phi <- function(x) x*log(x)
	A = B = 0
	for(i in 1:n){
		C1 = (1 - Pbi(i))/(1 - Qbi(i))
		if(C1>0) A = A + phi(C1)*(1 - Qbi(i))
		C2 = Pbi(i)/Qbi(i)
		if(C2>0) B = B * phi(C2)*Qbi(i)
	}

  }

  MED = omega*A + (1-omega)*B
  return(MED)
}


### Chaotic maps

# Logistic map

logistic <- function(c,n=200){
  x <- rep(NA,n)
  x[1]=0.1
  for (i in 2:n){
    x[i] = c * x[i-1] * (1 - x[i-1])
  }
  return(x)
}

n = 1000
c <- seq(0,4,0.05)
E <- ExGS <- ExSH <- rep(NA,length(c))

for(i in 1:length(c)){
	Z <- logistic(c[i],n) 
	opd = ordinal_pattern_distribution(x = Z, ndemb = 3)
	E[i] = permutation_entropy(opd)
	ExGS[i] = permutation_extropy(opd, type="GS")
	ExSH[i] = permutation_extropy(opd, type="SH")
}

E
ExGS
ExSH

plot(c, ExSH, type="l", lwd=2, col="red", ylab="Measures", main="a)")
lines(c, ExGS, col="blue", lwd=2)
lines(c, E, lwd=2)
legend("topleft", c("PE","GS-PEx", "PEx"), col=c("black","blue","red"), lwd=2, bty="n")

n = 1000
c <- seq(1.2,3.345,0.05)
D <- MEDGS <- MEDSH <- matrix(NA,length(c),length(c))

for(i in 1:length(c)) for(j in 1:length(c)){
	Z1 <- logistic(c[i],n) 
	opd1 = ordinal_pattern_distribution(x = Z1, ndemb = 3)
	Z2 <- logistic(c[j],n)
	opd2 = ordinal_pattern_distribution(x = Z2, ndemb = 3)
	MEDGS[i,j] = permutation_MEDdivergence(opd1,opd2,omega=0.5,type="GS") #/10
	MEDSH[i,j] = permutation_MEDdivergence(opd1,opd2,omega=0.5,type="SH") #/10
}

image2D(MEDGS,c,c,xlab=expression(c[1]),ylab=expression(c[2]),main="b)")
image2D(MEDSH,c,c,xlab=expression(c[1]),ylab=expression(c[2]),main="c)")

# Chevyshev

Chevy <- function(a,n=200){
  x <- rep(NA,n)
  x[1]=0.3
  for (i in 2:n){
    x[i]=cos((a^2)*acos(x[i-1]))
  }
  return(x)
}

n = 1000
a <- seq(0.1,5,0.05)
E <- ExGS <- ExSH <- rep(NA,length(c))

for(i in 1:length(a)){
	Z <- Chevy(a[i],n) 
	opd = ordinal_pattern_distribution(x = Z, ndemb = 3)
	E[i] = permutation_entropy(opd)
	ExGS[i] = permutation_extropy(opd, type="GS")
	ExSH[i] = permutation_extropy(opd, type="SH")
}

E
ExGS
ExSH

plot(a, E, type="l", lwd=2, col="black", ylab="Measures", main="a)")
lines(a, ExGS, col="blue", lwd=2)
lines(a, ExSH, lwd=2, col="red")
legend("bottomright", c("PE","GS-PEx", "PEx"), 
col=c("black","blue","red"), lwd=2, bty="n")

n = 1000
a <- seq(1.7,3.8,0.05)
MEDGS <- MEDSH <- matrix(NA,length(a),length(a))

for(i in 1:length(c)) for(j in 1:length(c)){
	Z1 <- Chevy(a[i],n) 
	opd1 = ordinal_pattern_distribution(x = Z1, ndemb = 3)
	Z2 <- Chevy(a[j],n) 
	opd2 = ordinal_pattern_distribution(x = Z2, ndemb = 3)
	MEDGS[i,j] = permutation_MEDdivergence(opd1,opd2,omega=0.5,type="GS") #/10
	MEDSH[i,j] = permutation_MEDdivergence(opd1,opd2,omega=0.5,type="SH") #/10
}

image2D(MEDGS,a,a,xlab=expression(c[1]),ylab=expression(c[2]),main="b)")
image2D(MEDSH,a,a,xlab=expression(c[1]),ylab=expression(c[2]),main="c)")



### Covid-19 data

data <- read.csv("WHO-COVID-19-global-daily-data.csv", sep=",", header=TRUE)
attach(data)
names(data)

table(data$Country)
A = data[data$Country == "India", c(1,7)]
fix(A)
names(A)
year <- substring(A$Date_reported, 1, last = 4)
A2 <- data.frame(year=as.numeric(year),nd=as.numeric(A[,2]))
A2

x = A2[A2$year==2022,2]
opd = ordinal_pattern_distribution(x = x, ndemb = 3)
permutation_entropy(opd)
permutation_extropy(opd)

### Divergences x year

Ind = data[data$Country == "India", c(1,7)]
year <- substring(Ind$Date_reported, 1, last = 4)
A.Ind <- data.frame(year=as.numeric(year),nd=as.numeric(Ind[,2]))

USA = data[data$Country == "United States of America", c(1,7)]
year <- substring(USA$Date_reported, 1, last = 4)
A.USA <- data.frame(year=as.numeric(year),nd=as.numeric(USA[,2]))

Chin = data[data$Country == "China", c(1,7)]
year <- substring(Chin$Date_reported, 1, last = 4)
A.Chin <- data.frame(year=as.numeric(year),nd=as.numeric(Chin[,2]))

# 2020:

x.Ind = A.Ind[A.Ind$year==2020,2]
opd.Ind = ordinal_pattern_distribution(x = x.Ind, ndemb = 3)
x.USA = A.USA[A.USA$year==2020,2]
opd.USA = ordinal_pattern_distribution(x = x.USA, ndemb = 3)
x.China = A.Chin[A.Chin$year==2020,2]
opd.China = ordinal_pattern_distribution(x = x.China, ndemb = 3)

GS20_1 = permutation_extropy(opd.Ind, type="GS")
GS20_2 = permutation_extropy(opd.USA, type="GS")
GS20_3 = permutation_extropy(opd.China, type="GS")

GS20_12 = permutation_MEDdivergence(opd.Ind,opd.USA,omega=0.5,type="GS")
GS20_13 = permutation_MEDdivergence(opd.Ind,opd.China,omega=0.5,type="GS")
GS20_23 = permutation_MEDdivergence(opd.USA,opd.China,omega=0.5,type="GS")

# 2021:

x.Ind = A.Ind[A.Ind$year==2021,2]
opd.Ind = ordinal_pattern_distribution(x = x.Ind, ndemb = 3)
x.USA = A.USA[A.USA$year==2021,2]
opd.USA = ordinal_pattern_distribution(x = x.USA, ndemb = 3)
x.China = A.Chin[A.Chin$year==2021,2]
opd.China = ordinal_pattern_distribution(x = x.China, ndemb = 3)

GS21_1 = permutation_extropy(opd.Ind, type="GS")
GS21_2 = permutation_extropy(opd.USA, type="GS")
GS21_3 = permutation_extropy(opd.China, type="GS")

GS21_12 = permutation_MEDdivergence(opd.Ind,opd.USA,omega=0.5,type="GS")
GS21_13 = permutation_MEDdivergence(opd.Ind,opd.China,omega=0.5,type="GS")
GS21_23 = permutation_MEDdivergence(opd.USA,opd.China,omega=0.5,type="GS")


# 2022:

x.Ind = A.Ind[A.Ind$year==2022,2]
opd.Ind = ordinal_pattern_distribution(x = x.Ind, ndemb = 3)
x.USA = A.USA[A.USA$year==2022,2]
opd.USA = ordinal_pattern_distribution(x = x.USA, ndemb = 3)
x.China = A.Chin[A.Chin$year==2022,2]
opd.China = ordinal_pattern_distribution(x = x.China, ndemb = 3)

GS22_1 = permutation_extropy(opd.Ind, type="GS")
GS22_2 = permutation_extropy(opd.USA, type="GS")
GS22_3 = permutation_extropy(opd.China, type="GS")

GS22_12 = permutation_MEDdivergence(opd.Ind,opd.USA,omega=0.5,type="GS")
GS22_13 = permutation_MEDdivergence(opd.Ind,opd.China,omega=0.5,type="GS")
GS22_23 = permutation_MEDdivergence(opd.USA,opd.China,omega=0.5,type="GS")

### Plots

years <- c(2020,2021,2022)

A <- rbind(c(GS20_1,GS21_1,GS22_1),  
c(GS20_2,GS21_2,GS22_2),
c(GS20_3,GS21_3,GS22_3)
)

plot(years, A[1,],type="o", ylim=c(0.6,1), xaxt='n', ylab="PExt-GS", lwd=2, 
main="a)")
axis(1,years)
lines(years, A[2,],type="o",col="blue",lwd=2)
lines(years, A[3,],type="o",col="red",lwd=2)
legend("top", c("India","USA","China"), col=c("black","blue","red"), 
lwd=2, horiz = TRUE, bty='n')


A <- rbind(c(GS20_12,GS21_12,GS22_12),  
c(GS20_13,GS21_13,GS22_13),
c(GS20_23,GS21_23,GS22_23)
)

plot(years, A[1,],type="o", ylim=c(370,1310), xaxt='n', 
ylab=expression(MED[phi]), lwd=2, main="b)")
axis(1,years)
lines(years, A[2,],type="o",col="blue",lwd=2)
lines(years, A[3,],type="o",col="red",lwd=2)
legend("topright", c("India-USA","India-China","USA-China"), 
col=c("black","blue","red"), lwd=2, bty='n')

