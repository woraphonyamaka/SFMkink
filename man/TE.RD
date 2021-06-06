\name{TE}
\alias{TE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Technical efficiency measure.
}
\description{
Computing and plotting the technical efficiency.
}
\usage{
TE(theta,Y,X,z=NULL,family,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{theta}{ The estimated parameters form the model}
\item{Y}{ vector of dependent variable}
\item{X}{ matrix of independent variable with kink effect}
\item{z}{ matrix of independent variable without kink effect}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
\item{type}{ nonlinear structure, type=c("kink", "threshold")"}
}
\details{
Computing and plotting the technical efficiency.
}
\value{
item{Output }{Technical efficiency series}
item{plot }{Plot of Technical efficiency}
}
\references{
Maneejuk, P., Yamaka, W., & Sriboonchitta, S. (2017). Analysis of global competitiveness using copula-based stochastic frontier kink model. In Robustness in Econometrics (pp. 543-559). Springer, Cham.
}
\author{
Worphon Yamaka and Paravee Maneejuk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
### TE
theta=model$result[,1]
tekink=TE1(theta,Y,X,z=NULL,family=family,type=type)
tekink=rev(sort(tekink))
te=rev(sort(te))
plot(tekink,lty=1,col="white",xlab = 'observation', ylab = 'TE value', main = "Technical Efficiency", ylim=c(0,1))
lines(tekink, lty=1,type="l",col="blue",lwd=2)
lines(te, lty=1,type="l",col="red", lwd=2)
lines(te, lty=1,type="l",col="black", lwd=2)

legend("topright",merge=FALSE,fill=c("red","blue","black"), col=c("red","blue","black"),
text.col=c("red","blue","black"),lwd=2,legend=c("Kink","Threshold","Linear"),cex=1, bg='gray90')
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
