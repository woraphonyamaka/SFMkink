#rm(list=ls())
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
############## Simulation
inilamda_3<-function(Datak=Z12)
{
  n=dim(Datak)[1]
  k=2
  lamdaini<-rep(0,k)
  ubskw<-rep(0,k)
  for(i in 1:k){ubskw[i]<-skewness(Datak[,i])*(n^2)/((n-1)*(n-2))}#
  deltaini<-sign(ubskw)*sqrt((pi/2)*(abs(ubskw))^(2/3)/((abs(ubskw))^(2/3)+((4-pi)/2)^(2/3)))
  for(i in 1:k){
    if(deltaini[i]>=1){deltaini[i]<-0.99}
    if(deltaini[i]<=-1){deltaini[i]<--0.99}
  }
  lamdaini<-deltaini/sqrt(1-deltaini^2)
  lamdaini
}

simCOP_sn=function(rho1=par_cop,lambda1=1,lambda2=2,nsamp,copfam=1){
  ## skew normal copoula data gen
  Rmat<-matrix(c(1,rho1,rho1,1),2,2)
  RLkd<-Rmat
  lamda=c(lambda1,lambda2)
  ###
  lamdaLkd<-lamda
  DeltaLkd<-diag(1/sqrt(1+lamda^2))
  skewLkd<-t(lamdaLkd)%*%solve(RLkd)%*%solve(DeltaLkd)/sqrt(1+(t(lamdaLkd)%*%solve(RLkd)%*%lamdaLkd)[1,1]) #row vector
  alphaLkd<-t(skewLkd)
  OmegaLkd<-DeltaLkd%*%(RLkd+lamdaLkd%*%t(lamdaLkd))%*%DeltaLkd
  OmegaLkd<-as.matrix(OmegaLkd)
  Z12=cbind(rep(0,nsamp),rep(0,nsamp))
  for ( j in 1:nsamp){
    Z12[j,]<-rmsn(n=1, xi=c(0,0), OmegaLkd, alphaLkd,  tau=0, dp=NULL)
  }
  lambda_est<-inilamda_3(Z12)

  Fu=psn(Z12[,1], xi = 0, omega=1, alpha=lamdaLkd[1])
  Fv=psn(Z12[,2], xi = 0, omega=1, alpha=lamdaLkd[2])
  #Fuv=BiCopSim(N=nsamp, family = copfam, par = rho1,par2 = 0)
  ##############################
  res=cbind(Fu,Fv)
  return(res)
}

neg.part <- function(x) x*(x<0)  # Regime 1
pos.part <- function(x) x*(x>0)  # Regime 2

Tlow <- function(x,r) x*(x-r<0)  # Regime 1
Tup  <- function(x,r) x*(x-r>0)   # Regime 2


sfa.simu<-function(nob,alpha,sigV,sigU,kink,family,rho,type)
{
  nmax<-nob

  if (family==2){
    sim = BiCopSim(nob,family=family,rho,4)
  }else if (family=="sn"){
    sim = simCOP_sn(rho1=rho,lambda1=1,lambda2=2,nsamp=nob,copfam=1)
  }else{
    sim =BiCopSim(nob,family=family,rho)
  }

  u=qnorm(sim[,1],0, sd=sigU)
  v=qtruncnorm(sim[,2],a=0, b=Inf, mean = 0, sd = sigV)
  k=1

  x = matrix(rnorm(n*k,kink,2),ncol=k)
  if( type=="kink"){
    x1 = cbind(1,neg.part(x-kink),pos.part(x-kink))
  }else {
    x1 = cbind(1,Tlow (x,kink),Tup (x,kink))}

  y =c(t(alpha)%*%t(x1)+u-v)

  out=list(Y=y,X=x)
  return(out)
}
##===========================================

