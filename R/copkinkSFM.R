# This is an  function named 'Copula based kink Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Main Function
copkinkSFM=function(Y,X,z=NULL,family,RHO,LB,UB,type){

  like<-function(theta,Y,X,z=NULL,family, type){
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
    hw=sum(log(1/m*diag(gv%*%gaucopula)))
    if (is.infinite(hw))  # control for optimization
      hw<--n*100
    if (is.nan(hw))  # control for optimization
      hw<--n*100

    cat("Sum of log Likelihood for kink-SFA ->",sprintf("%4.4f",hw),"\n")


    return(hw) # log likelihood

  }

  ### End Function #############3
  #=================================================


  # start here
  # select familty  copula upper and lower bouubd ( look at CDVine package)
  family=family  # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....
  LB=LB      #lower bound
  UB=UB     # upper bound
  RHO=RHO   # any value in bound


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

  if(type=="kink"|type=="threshold"){
    kink=colMeans(X)
    cc=rep(0.1,K)
    if (family==2){
      lower =c(rep(-Inf,K),0.01,0.01,colMins(X)+abs(colMeans(X)*0.01),LB,4.1)
      upper =c(rep(Inf,K+2),colMaxs(X)-abs(colMeans(X)*0.01),UB,50)
      start0=c(cc,sigmau=1,sigmav=1,kink,rho=RHO,df=4)
    } else if (family=="sn"){
      lower =c(rep(-Inf,K),0.01,0.01,colMins(X)+abs(colMeans(X)*0.01),LB,0,0)
      upper =c(rep(Inf,K+2),colMaxs(X)-abs(colMeans(X)*0.01),UB,10,10)
      start0=c(cc,sigmav=1,sigmau=1,kink,rho=RHO,lam1=1, lam2=2)
    } else if (family==0){
      lower =c(rep(-Inf,K),0.01,0.01,colMins(X)+abs(colMeans(X)*0.01))
      upper =c(rep(Inf,K+2),colMaxs(X)-abs(colMeans(X)*0.01))
      start0=c(cc,sigmav=1,sigmau=1, kink)
    }else{
      lower =c(rep(-Inf,K),0.01,0.01,colMins(X)+abs(colMeans(X)*0.01),LB)
      upper =c(rep(Inf,K+2),colMaxs(X)-abs(colMeans(X)*0.01),UB)
      start0=c(cc,sigmau=1,sigmav=1,kink,rho=RHO)
    }
  }else{
    cc=rep(0.1,K)
    if (family==2){
      lower =c(rep(-Inf,K),0.01,0.01,LB,4.1)
      upper =c(rep(Inf,K+2),UB,50)
      start0=c(cc,sigmau=1,sigmav=1,rho=RHO,df=4)
    } else if (family=="sn"){
      lower =c(rep(-Inf,K),0.01,0.01,LB,0,0)
      upper =c(rep(Inf,K+2),UB,10,10)
      start0=c(cc,sigmav=1,sigmau=1,rho=RHO,lam1=1, lam2=2)
    }else{
      lower =c(rep(-Inf,K),0.01,0.01,LB)
      upper =c(rep(Inf,K+2),UB)
      start0=c(cc,sigmau=1,sigmav=1,rho=RHO)
    }
  }

  theta=start0
  model <- optim(start0,like,Y=Y,X=X,z=z,family=family,type=type,
                 control = list(maxit=100000,fnscale=-1),method="L-BFGS-B",
                 lower =lower,upper =upper, hessian=TRUE )

  if (model$message=="ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"){

    model$hessian[,(K+3)]=rnorm(length(start0))
    model$hessian[(K+3),]=rnorm(length(start0))
  }
  # table of results
  coef<- model$par
  k=length(coef)
  model$se <- sqrt(-diag(solve(model$hessian)))

  for(i in 1:k){
    if (is.nan(model$se[i]))  # control for optimization
      model$se[i] <- sqrt(-diag(solve(-model$hessian)))[i]
  }

  n=length(Y)
  S.E.= model$se
  (paramsWithTs = cbind (model$par , coef/S.E. ) )
  stat=coef/S.E.
  pvalue <- 2*(1 - pnorm(abs(stat)))
  result <- cbind(coef,S.E.,stat,pvalue)
  result
  BIC= -2*model$value+ (log(n)*length(coef))
  AIC = -2*model$value + 2*length(coef)


  output=list(
    result=result,
    AIC=AIC,
    BIC=BIC,
    Loglikelihood=model$value
  )
  output

}

#########END ALL functions
