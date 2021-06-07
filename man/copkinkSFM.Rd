\name{copkinkSFM}
\alias{copkinkSFM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copula based nonlinear Stochastic frontier Model
}
\description{
In the standard stochastic frontier model, the two-sided error term V and the one-sided technical inefficiency error term W are assumed to be independent. In this paper, we relax this assumption by modeling the dependence between V and W using copulas. Nine copula families are considered and their parameters are estimated using maximum simulated likelihood.Also the nonlinear or kink effect of the model is provided.
}
\usage{
copkinkSFM(Y,X,z=NULL,family,RHO,LB,UB,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y}{ vector of dependent variable}
 \item{X}{ matrix of independent variable with kink effect}
\item{z}{ matrix of independent variable without kink effect}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2,... (see, Vinecopula package)}
\item{RHO}{ The initail value of the copula parameter}
\item{LB}{  The lower bound of the copula parameter}
\item{UB}{ The upper bound of the copula parameter}
\item{type}{ nonlinear structure, type=c("kink", "threshold")}

}
\details{
This model is used to explain both nonlinear and linear effects of independent variables on dependent variable)

}
\value{
\item{result }{The result contain the estimated parameters, standard errors, t-stat, and p-value}
\item{AIC }{Akaiki Information Criteria}
\item{BIC }{Bayesian Information Criteria}
\item{Loglikelihood }{Maximum Log-likelihood function}
}
\references{
Maneejuk, P., Yamaka, W., & Sriboonchitta, S. (2017). Analysis of global competitiveness using copula-based stochastic frontier kink model. In Robustness in Econometrics (pp. 543-559). Springer, Cham.
}
\author{
Woraphon Yamaka and Paravee Maneejuk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Required packages
library(truncnorm)
library(mvtnorm)
library("VineCopula")
library("frontier")
library(CDVine)
library(copula)
library(matrixStats)
library(matrixcalc)
library(sn)
library(e1071)
library(randtoolbox)
library("signal")
library(numDeriv)
library(moments)
library("rngWELL")
library(Matrix)
library("signal")
library("mnormt")
library("modeest")



#example simulation data

set.seed(2125)
nob=100
alpha=c(0.5,-0.8,0.7)
sim=sfa.simu(nob=nob,alpha=alpha,sigV=0.5,sigU=1,kink=0.7,family=1,rho=0.8,type="kink")
data=data.frame(sim$Y,sim$X)

# Select familty  copula upper and lower bouubd ( look at CDVine package)
# family=1   # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....
#Gaussian (-.99, .99)
#Student t (-.99, .99)
#Clayton (0.1, Inf)
#sn  (-.99, .99)


family=1  # normal=1 , Ind=0
type="kink"  # threshold, "kink"", "linear""
Y=data[,1]
z=NULL
X=as.matrix(data[,2])
plot(X,Y)

model=copkinkSFM(Y=Y,X=X,family=family,z=NULL,RHO=0.5,LB=-0.99,UB=0.99, type=type)
model
model0 <- sfa( Y~ X)
summary( model0 )

##plot kink -------------
B=model$result[1:3]
bt1=c(B[1],B[2],B[3])
thres=model$result[,1][6]
x=X[,1]
dx = seq(min(x),max(x),by=0.01)		           # Grid on regression function for display
G1 = cbind(1,neg.part(dx-thres),pos.part(dx-thres))  # for kink
yf1 = G1%*%bt1

#The parameter estimates from this fitted model
par(bg = 'gray99')
plot(x,Y, xlim = c(min(x),max(x)),pch =19,col="orange",
 xlab = 'X', ylab = 'Y', main = 'Simulation n=300')
lines(dx,ts(yf1), lty=20,type="l",lwd=6,col="2")
grid()

minE=min(yf1)
pointx=thres
points(pointx,minE,col="blue",bg="blue",pch=22,lwd=5)



##plot Thresold -------------
B=model$result[1:3]
bt1=c(B[1],B[2],B[3])
thres=model$result[,1][6]
x=X[,1]
dx = seq(min(x),max(x),by=0.01)		           # Grid on regression function for display
N1=length(dx)
G1 = cbind(1,Tlow(dx,thres))        	       # for kink
G2 = cbind(1,Tup(dx,thres))          		 # for kink
G = cbind(1,Tlow(dx,thres),Tup(dx,thres))  # for kink
yf = c(G%*%bt1)

N=which(round(dx,2)==round(thres,2))-1
dx1 = dx[1:N]
dx2 = dx[(N+1):N1]
#The parameter estimates from this fitted model
par(bg = 'gray99')
plot(x,Y, xlim = c(min(x),max(x)),pch =19,col="orange",
 xlab = 'X', ylab = 'Y', main = 'Simulation n=100')
lines(dx1,ts(yf)[1:N], lty=20,type="l",lwd=6,col="2")
lines(dx2,ts(yf)[(N+1):N1]	, lty=20,type="l",lwd=6,col="2")
grid()


##plot linear -------------
B=model$result[1:2]
bt1=c(B[1],B[2])
x=z[,1]
dx = seq(min(x),max(x),by=0.01)		           # Grid on regression function for display
G1 = cbind(1,dx)  # for kink
yf1 = G1%*%bt1

#The parameter estimates from this fitted model
par(bg = 'gray99')
plot(x,Y, xlim = c(min(x),max(x)),pch =19,col="orange",
 xlab = 'X', ylab = 'Y', main = 'Simulation n=100')
lines(dx,ts(yf1), lty=20,type="l",lwd=6,col="2")
grid()


### TE
theta=model$result[,1]
tekink=TE(theta,Y,X,z=NULL,family=family,type=type)
tekink=rev(sort(tekink))
plot(tekink,lty=1,col="white",xlab = 'observation', ylab = 'TE value', main = "Technical Efficiency", ylim=c(0,1))
lines(tekink, lty=1,type="l",col="blue",lwd=2)
legend("topright",merge=FALSE,fill=c("blue"), col=c("blue"),
text.col=c("red","blue","black"),lwd=2,legend=c("Kink"),cex=1, bg='gray90')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
