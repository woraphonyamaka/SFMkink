
#sn copula density


### from "sn" ver.1.2-0 ###
T.Owen <- function(h, a, jmax=50, cut.point=8){
  T.int <-function(h, a, jmax, cut.point){
    fui <- function(h,i) (h^(2*i))/((2^i)*gamma(i+1))
    seriesL <- seriesH <- NULL
    i <- 0:jmax
    low<- (h <= cut.point)
    hL <- h[low]
    hH <- h[!low]
    L <- length(hL)
    if (L > 0) {
      b <- outer(hL, i, fui)
      cumb <- apply(b, 1, cumsum)
      b1 <- exp(-0.5*hL^2) * t(cumb)
      matr <- matrix(1, jmax+1, L) - t(b1)
      jk <- rep(c(1,-1), jmax)[1:(jmax+1)]/(2*i+1)
      matr <- t(matr*jk) %*% a^(2*i+1)
      seriesL <- (atan(a) - as.vector(matr))/(2*pi)
    }
    if (length(hH) > 0) seriesH <-
      atan(a)*exp(-0.5*(hH^2)*a/atan(a)) * (1+0.00868*(hH*a)^4)/(2*pi)
    series <- c(seriesL, seriesH)
    id <- c((1:length(h))[low],(1:length(h))[!low])
    series[id] <- series # re-sets in original order
    series
  }
  if(!is.vector(a) | length(a)>1) stop("'a' must be a vector of length 1")
  if(!is.vector(h)) stop("'h' must be a vector")
  aa <- abs(a)
  ah <- abs(h)
  if(is.na(aa)) stop("parameter 'a' is NA")
  if(aa==Inf) return(sign(a)*0.5*pnorm(-ah)) # sign(a): 16.07.2007
  if(aa==0) return(rep(0,length(h)))
  na <- is.na(h)
  inf <- (ah == Inf)
  ah <- replace(ah,(na|inf),0)
  if(aa <= 1)
    owen <- T.int(ah,aa,jmax,cut.point)
  else
    owen<- (0.5*pnorm(ah) + pnorm(aa*ah)*(0.5-pnorm(ah))
            - T.int(aa*ah,(1/aa),jmax,cut.point))
  owen <- replace(owen,na,NA)
  owen <- replace(owen,inf,0)
  return(owen*sign(a))
}

########
## quantile function of standard skew-Normal distribution
## using "NR" solver
qssn <- function(p, zeta = 0, tol = 1e-08, ...){
  max.q <- sqrt(qchisq(p,1));
  min.q <- -sqrt(qchisq(1-p,1));
  if(zeta == Inf) return(as.numeric(max.q))
  if(zeta == -Inf) return(as.numeric(min.q))
  na <- is.na(p) | (p < 0) | (p > 1)
  zero <- (p == 0)
  one <- (p == 1)
  p <- replace(p, (na | zero | one), 0.5)
  delta <- zeta/sqrt(1+zeta^2);
  ex <- sqrt(2/pi)*delta;
  vx <- 1-ex^2;
  sx <- 0.5*(4-pi)*(ex^3);
  kx <- 2*(pi-3)*(ex^4);
  g1 <- sx/(vx^(3/2));
  g2 <- kx/(vx^2);
  x <- qnorm(p)
  x <- (x + (x^2 - 1) * g1/6 + x * (x^2 - 3) * g2/24 -
          x * (2 * x^2 - 5) * g1^2/36)
  x <- ex + sqrt(vx) * x
  px <- pssn(x, zeta=zeta, ...);
  max.err <- 1
  while (max.err > tol) {
    logPDF <- ldssn(x,zeta);
    x1 <- x - (px - p)/exp(logPDF);
    x <- x1
    px <- pssn(x, zeta=zeta, ...)
    max.err <- max(abs(px-p))
    if(is.na(max.err)) stop('failed convergence, try with solver="RFB"')
  }
  x <- replace(x, na, NA)
  x <- replace(x, zero, -Inf)
  x <- replace(x, one, Inf)
  q <- as.numeric(x)
  names(q) <- names(p)
  return(q)
}
## cumulative density function of standard skew-Normal distribution
pssn <- function(x, zeta=0, engine, ...){
  nx <- length(x)
  na <- length(zeta)
  if(missing(engine)) engine <-
    if(na == 1 & nx > 3 & all(zeta*x > -5))
      "T.Owen" else "biv.nt.prob"
  if(engine == "T.Owen") {
    if(na > 1) stop("engine='T.Owen' not compatible with other arguments")
    p <- pnorm(x) - 2 * T.Owen(x, zeta, ...)
  }
  else{ # engine="biv.nt.prob"
    p <- numeric(nx)
    zeta <- cbind(x, zeta)[,2]
    delta <- zeta/sqrt(1+zeta*zeta)
    for(k in seq_len(nx)) {
      if(abs(zeta[k]) == Inf){
        p[k] <- if(zeta[k] > 0)
          2*pnorm(pmax(x[k],0)) - 1
        else
          2*pnorm(pmin(x[k],0))
      }
      else { # SNbook: formula (2.48), p.40
        R <- matrix(c(1, -delta[k], -delta[k], 1), 2, 2)
        p[k] <- 2 * biv.nt.prob(0, rep(-Inf,2), c(x[k], 0), rep(0, 2), R)
      }
    }
  }
  pmin(1, pmax(0, as.numeric(p)))
}
## log-density of standard univariate skew-Normal distribution
ldssn <- function(x,zeta){
  logN <- (-log(sqrt(2*pi)) -x^2/2);
  if(abs(zeta) < Inf)
    logS <- pnorm(zeta*x, log.p=TRUE)
  else
    logS <- log(as.numeric(sign(zeta)*x > 0))
  as.numeric(logN + logS - log(0.5));
}
## interpolating quantiles for skew-Normal copula
ipqssn <- function(udat,zeta,mpoints=150){
  dim <- ncol(udat);
  ix <- matrix(0,nrow=nrow(udat),ncol=dim);
  for(j in 1:dim){
    #j=2
    minmaxu=c(0,1)

    minmaxu <- c(min(udat[,j]),max(udat[,j]));


    minmaxx <- qssn(minmaxu, zeta=zeta[j]);

    #minmaxx[minmaxx==Inf]=maxv
    if(minmaxx[1]==minmaxx[2]){ix[,j]=minmaxx[1]}else{
      xx <- seq(minmaxx[1],minmaxx[2],length.out=mpoints);
      px <- sort(pssn(xx, zeta[j]));
      ix[,j] <- pchip(px, xx, udat[,j]);}
  }
  ix
}

aqssn <- function(udat,zeta){
  dim <- ncol(udat);
  ax <- matrix(0,nrow=nrow(udat),ncol=dim);
  for(j in 1:dim){
    ax[,j] <- qssn(udat[,j], zeta=zeta[j]);
  }
  ax
}
#Fisher transformation for each spherical parameter
trans_Shp<-function(y){pi*exp(y)/(1+exp(y)) } #-inf,inf to 0,pi
backtrans_Shp<-function(y){log(y/(pi-y))} #0,pi to -inf,inf

### delta in Omega using Fisher transformation
trans<-function(y){(exp(y)-1)/(exp(y)+1)} #-inf,inf to -1,1
backtrans<-function(y){log((1+y)/(1-y))} #-1,1 to -inf,inf

#transform R into unconstrained vector

Rtouncvec<-function(R){
  initial_L<-chol(R)
  initialL_shp<-sphcoord(initial_L)
  iniL_uncon<-rep(NA,num)
  for(j in 2:k){for(i in 2:j){iniL_uncon[(j-1)*(j-2)/2+i-1]<-initialL_shp[i,j]}}
  iniL_uncon<-backtrans_Shp(iniL_uncon)
  iniL_uncon
}


#Exchangable R chol decomposition
#give a cholesky decomposition upper triangular matrix for exchangeable R matrix
rho2Rmat<-function(rho){
  L<-matrix(0,k,k)
  L[1,1]=1
  for(t in 1:k){L[t,1]=rho}
  for(j in 2:k){
    for(i in 2:(j)){
      if(i<j){
        cons<-0
        for(s in 1:(i-1)){cons<-cons+L[j,s]*L[i,s]}
        L[j,i]<-(rho-cons)/L[i,i]}

      if(i==j){
        cons1<-0
        for(s in 1:(i-1)){cons1<-cons1+L[j,s]^2}
        L[i,i]=sqrt(1-cons1)
      }
    }
  }
  L[1,1]=1
  L=t(L)
}
#t(rho2Rmat(-0.8))%*%(rho2Rmat(0.2))
#det(rho2Rmat(-1))
#TL<-rho2Rmat(0.9)

#unconstrained vector into R(inverse of the function Rtouncvec)

unconvtoR<-function(L_uncon){
  #L_uncon<-symtovec(R)
  #convert back to the correlation matrix
  ini_sph<-trans_Shp(L_uncon)#initiall_shp
  LM<-unconv2LM(ini_sph)
  Rest<-t(LM)%*%LM
  Rest
}

#transform spherical vector(Angels with range in [0,pi]) into L matrix which is need in the unconvtoR function
unconv2LM<-function(ini_sph1)
{
  LM<-matrix(0,k,k)
  n<-1
  for(j in 2:k)
  {
    ind<-c(n:(n+j-2))#index set of each column
    LM[1:j,j]<-t(polar2rect(1,ini_sph1[ind]))
    n=n+(j-1)#change the pointer to the next column index
    #cat(ind,"\n")
  }
  LM[1,1]<-1
  LM
}


#transform lamda to delta
deltafun<-function(lamdavec){lamdavec/sqrt(1+lamdavec^2)}
lamdafun<-function(deltavec){deltavec/sqrt(1+deltavec^2)}
########

#rho is unconstraint rho

skewnormCopden<-function(theta1,uv){
  #lamda=c(lambda1,lambda2)
  lamda=theta1[2:3]
  rho=theta1[1]
  k <- length(lamda)
  if ( max(deltafun(lamda)) > 0.999 | min(deltafun(lamda)) < -0.999){
    lamda<-lamdafun(runif(k,-0.8,0.8))
  }

  lamdaLkd<-lamda
  #RLkd R matrix for numerator in likelihood function
  #RLkd<-matrix(c(1,rho,rho,1),2,2)
  num=k*(k-1)/2
  Rmat<-matrix(c(1,trans(rho),trans(rho),1),2,2)
  RLkd<-Rmat
  #DeltaLkd Delta matrix for numerator in likelihood function
  DeltaLkd<-diag(1/sqrt(1+lamda^2))
  skewLkd<-t(lamdaLkd)%*%solve(RLkd)%*%solve(DeltaLkd)/sqrt(1+(t(lamdaLkd)%*%solve(RLkd)%*%lamdaLkd)[1,1]) #row vector
  alphaLkd<-t(skewLkd)
  OmegaLkd<-DeltaLkd%*%(RLkd+lamdaLkd%*%t(lamdaLkd))%*%DeltaLkd
  OmegaLkd<-as.matrix(OmegaLkd)
  if(!(is.symmetric.matrix(OmegaLkd))){OmegaLkd<-as.matrix(forceSymmetric(OmegaLkd))}#due to machine precision sometimes it is not a symmetric matrix

  u1<-uv
  x<-ipqssn(u1,zeta=lamdaLkd)
  #x<-aqssn(u1,zeta=lamdaLkd)
  #x1=qsn(u1[,1], xi=0, omega=1, alpha=lamdaLkd[1], tol=1e-8)
  #x2=qsn(u1[,2], xi=0, omega=1, alpha=lamdaLkd[2], tol=1e-8)
  #x<-cbind(x1,x2)
  #x<-qsn(zres, xi = xi[y_index], omega=Omega[y_index,y_index], alpha=lambda[y_index])
  val <- dmsn(x, xi = c(rep(0, k)),  Omega = OmegaLkd, alpha = alphaLkd)
  for (h in 1:k) val <- val / dsn(x[,h], xi = 0, omega=1, alpha=lamdaLkd[h])
  val[is.na(val)]<-0
  return(val)
}
# Copula likelihood
#skewnormCopden(theta = c(rho,lambda1,lambda2),uv=cbind(uvec,vvec))
