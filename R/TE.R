
# This is an  function named 'Copula based kink Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Compute technical efficiency
#
#EX: te1=TE1(coef,Y,X,family=family)


TE<-function(theta,Y,X,z=NULL,family,type){

  if (is.null(z)) {
    exo=NULL
    K2=0
  }else {
    exo=z
    K2=ncol(exo)}

  if (is.null(X)) {
    X=NULL
    K1=1
  }else {
    X=X
    K1=(ncol(X)*2)+1}

  K=K1+K2
  n=length(Y)
  sigmau=abs(theta[K+1])
  sigmav=abs(theta[K+2])
  T=ncol(X)
  if (type =="linear"){
    kink=NULL
  }else{
    kink=theta[c((K+3):(K+2+T))]}
  m=n


  if(type=="kink"){
    XX=matrix(0,n,T*2)
    for (j in 1:T){
      XX[,c(((2*j)-1):(2*j))] = cbind(neg.part(X[,j]-kink[j]),pos.part(X[,j]-kink[j]))}
  } else if (type=="threshold"){
    XX=matrix(0,n,T*2)
    for (j in 1:T){
      XX[,c(((2*j)-1):(2*j))] = cbind(Tlow(X[,j],kink[j]),Tup(X[,j],kink[j]))}
  }else{
    XX=NULL
    K=K1+K2
  }

  xmat=cbind(1,XX,exo)

  w=c(Y-t(theta[1:K])%*%t(xmat))

  set.seed(1988)
  u=replicate(n,rtruncnorm(n, a=0.001, b=Inf, mean = 0, sd = sigmau))
  W=t(replicate(n,w))
  gv=dnorm(u+W,mean=0,sd=sigmav)+0.000001
  gv=matrix(t(gv),nrow=n,ncol=n)
  Gv=pnorm(u+W,mean=0,sd=sigmav)
  Fu=ptruncnorm(u, a=0.0001, b=Inf, mean = 0, sd = sigmau)

  Fu=c(abs(Fu))
  Gv=c(abs(Gv))
  mm=length(Fu)
  for ( i in 1:mm){
    if (is.infinite(Fu[i]))  # control for optimization
      Fu=0.0000000000000001
    if (is.infinite(Gv[i]))  # control for optimization
      Gv=0.000000000000001
    if (is.nan(Fu[i]))  # control for optimization
      Fu=0.00000000000001
    if (is.nan(Gv[i]))  # control for optimization
      Gv=0.00000000000001
    if ((Fu[i]>0.9999))  # control for optimization
      Fu[i]=0.9999
    if ((Gv[i]>0.9999))  # control for optimization
      Gv[i]=0.9999
  }

  if (family==2){
    rho=theta[length(theta)-1]
    df=theta[length(theta)]
    gaucopula=BiCopPDF(Fu, Gv, family=family, par=rho, par2=df)+0.00000001
  } else if (family==0){
    gaucopula=BiCopPDF(Fu, Gv, family=0,par=0)+0.00000001
  } else if (family=="sn"){
    rho=theta[length(theta)-2]
    lam1=theta[length(theta)-1]
    lam2=theta[length(theta)]
    gaucopula=skewnormCopden(theta1 = c(rho,lam1,lam2),uv=cbind(Fu,Gv))
  }else{
    rho=theta[length(theta)]
    gaucopula=BiCopPDF(Fu, Gv, family=family, par=rho, par2=0)+0.00000001}

  gaucopula=matrix(gaucopula,nrow=m,ncol=n)

  hw=1/n*diag(gv%*%gaucopula) # A
  tu=sapply(u,mean)
  tcltcopula=sapply(gaucopula,mean)
  tgv=sapply(t(gv),mean)
  hw2=exp(-tu)*tcltcopula*tgv
  hw2=matrix(hw2,nrow=n,ncol=n)
  hw2=colMeans(hw2) # B
  # technical efficiency TE=A/B
  te=hw2/hw
  n=length(te)
  te1=rev(sort(te))
  plot(te1,lty=1,col="white",xlab = 'observation', ylab = 'TE value', main = "Technical Efficiency")
  lines(te1, lty=1,type="l",col="blue")
  return(te1)
}



