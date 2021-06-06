\name{sfa.simu}
\alias{sfa.simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simuation data from Copula based nonlinear Stochastic frontier Model
}
\description{
This function is used to simulate the dataset
}
\usage{
sfa.simu(nob,alpha,sigV,sigU,kink,family,rho,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{nob}{ Number of simulated data}
 \item{alpha}{ vector of true parameters}
\item{sigV}{ Variance of V}
\item{sigU}{ Variance of U}
\item{kink}{ kink or threshold parameter}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2,... (see, Vinecopula package)}
\item{rho}{  copula parameter}
\item{type}{ nonlinear structure, type=c("kink", "threshold")}
}
\details{
Simulation dataset
}
\value{
\item{out }{ Y=Dependent variable, x=matrix of independent variables with nonlinear effects}
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
##===========================================
library(truncnorm)
library(CDVine)
library(copula)
library(frontier)
library(matrixStats)
library(matrixcalc)
library(sn)
library(e1071)
library(randtoolbox)
set.seed(2125)
nob=300
n=nob
alpha=c(0.5,-0.8,0.7)
sim=sfa.simu(nob=nob,alpha=alpha,sigV=0.5,sigU=1,kink=0.7,family="sn",rho=0.8, type="kink")

sim=sfa.simu(nob=nob,alpha=alpha,sigV=0.5,sigU=1,kink=0.7,family="sn",rho=0.8, type="threshold")

data=data.frame(sim$Y,sim$X)
y=data[,1]
x=data[,2]
windows()
plot(x,y)

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
