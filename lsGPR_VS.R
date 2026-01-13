#--------------------------------------------------------------
# LSGPR WITH VARIABLE SELECTION
#-------------------############-------------------------------
# Component                 Code  
#--------------------------------------------------------------
# Inclusion indicators(γ)   gamma = rbinom(p,1,pro) 
# Cntinuous weights (ρ)     rho[k]=runif(1) if included
# Node importance (β)       beta = -log(rho)
# Spike-and-slab prior      if(gamma[k]==0){rho[k]=1} → β=0
# Between-model moves       (spike ↔ slab)
# Within-model moves        (when γ[k]=1)
# Prior inclusion probs     pro = rbeta(1, ...) updated each iteration
# Storage of samples        GAMMA[,g] = gamma
# Node-weighted distances   U.full.beta = lapply(U.full, ...)
#--------------------------------------------------------------

lsGPR_VS<-function(a,a.te,U,U.te,z=NULL,z.te=NULL,Y,init,burn,nrun,prior){

  n=length(c(a))
  nn=length(c(a.te))
  a.full=c(a,a.te)
  U.full=append(U,U.te)
  Ea.full = Eu.full = matrix(0,(n+nn),(n+nn))
  for(i in 2:(n+nn)){
    for(j in 1:(i-1)){
      Ea.full[i,j]=Ea.full[j,i]=(a.full[i]-a.full[j])^2
      Eu.full[i,j]=Eu.full[j,i]=norm(U.full[[i]]-U.full[[j]],"F")^2
    }
  }
  Ea = Ea.full[1:n,1:n]
  sd.a = sd(Ea)
  Ea.full = Ea.full/sd.a
  Ea = Ea.full[1:n,1:n]
  Eu = Eu.full[1:n,1:n]
  sd.u = sd(Eu)
  Eu.full = Eu.full/sd.u
  Eu = Eu.full[1:n,1:n]
  if(!is.null(z)){
    z.full=rbind(z,z.te)
    Ez.full = matrix(0,(n+nn),(n+nn))
    for(i in 2:(n+nn)){
      for(j in 1:(i-1)){
        Ez.full[i,j]=Ez.full[j,i]=sum((z.full[i,]-z.full[j,])^2)
      }
    }
    Ez = Ez.full[1:n,1:n]
    sd.z = sd(Ez)
    Ez.full = Ez.full/sd.z
    Ez = Ez.full[1:n,1:n]
  }
  
  if(missing(init)){init=c(1,1,1,1)}
  if(missing(burn)){burn=5000}
  if(missing(nrun)){nrun=10000}
  
  # if(missing(prior)){prior=list(a.t=5,b.t=1,a.p1=2,b.p1=5,step=0.01,a.p=1,b.p=40)}
  # if(missing(prior)){prior=list(a.t=5,b.t=1,
  #                               a.p1=2,b.p1=50, #psi
  #                               step=0.01, #MH step size
  #                               a.p=1,b.p=25 #node selection
  # )}
  if(missing(prior)){prior=list(
  a.t = 15,    # tau prior
  b.t = 1,
  a.p1 = 2,   # psi1 prior
  b.p1 = 50,
  step = 0.01,  # MH step size
  a.p = 1,      # pi* prior (Beta(1, 25) for sparsity)
  b.p = 10
  )}
  
  # rinvgamma(1,shape=2,scale=50)
  # rbeta(1, shape1=1, shape2=5)
  
  
  tau=init[1]
  psi1=init[2]
  psia=init[3]
  psiu=init[4]
  phi=rep(0,n)
  PSI1=rep(NA,nrun); PSIA=rep(NA,nrun); PSIU=rep(NA,nrun); TAU=rep(NA,nrun)
  LOG=rep(NA,nrun)
  PHI = matrix(NA,n,nrun)
  PRED = matrix(NA,nn,nrun)
  if(!is.null(z)){
    psiz = 1
    PSIZ=rep(NA,nrun)
  }
  
  # add in node selection
  pro = 0.25 # hyperparameter, fixed / beta(1,40)
  #pro = mean(rbeta(10000, shape1=1, shape2=10))
  # pro = 0.25
  # rbeta(1, shape1=prior$a.p, shape2=prior$b.p)
  print(pro)
  # pro = 0.5
  #rinvgamma(1,shape=a.psi1,scale=b.psi1)

  p = nrow(U[[1]]); d = ncol(U[[1]])
  gamma = rbinom(p,1,pro)
  rho = rep(NA,p)
  for(k in 1:p){
    if(gamma[k]==1){rho[k]=runif(1)}
    if(gamma[k]==0){rho[k]=1}
  }
  beta = -log(rho)
  GAMMA = matrix(NA,p,nrun)
  BETA = matrix(NA,p,nrun)
  
  # update Eu and Eu.full
  U.full.beta = lapply(U.full,function(x) hadamard.prod(x,matrix(sqrt(beta),ncol=1)%*%matrix(rep(1,d),nrow=1)))
  Eu.full = matrix(0,(n+nn),(n+nn))
  for(i in 2:(n+nn)){
    for(j in 1:(i-1)){
      Eu.full[i,j]=Eu.full[j,i]=norm(U.full.beta[[i]]-U.full.beta[[j]],"F")^2
    }
  }
  Eu.full = Eu.full/sd.u
  Eu = Eu.full[1:n,1:n]
  pb <- txtProgressBar(min = 0, max = nrun, style = 3)
  #g=1
  for(g in 1:nrun){
    # update psi1
    if(is.null(z)){EE = exp(-psia*Ea-psiu*Eu) + diag(rep(1e-12,n))}
    if(!is.null(z)){EE = exp(-psia*Ea-psiu*Eu-psiz*Ez) + diag(rep(1e-12,n))}
    a.psi1 = prior$a.p1 + 0.5*n
    b.psi1 = prior$b.p1 + 0.5*(matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1))
    psi1 = rinvgamma(1,shape=a.psi1,scale=b.psi1)
    PSI1[g] = psi1
    
    # update psia
    ev1 = diag(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))^2
    log.det = sum(log(ev1))
    quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    rand = exp(rnorm(1,log(psia),prior$step))
    if(is.null(z)){E2 = exp(-rand*Ea-psiu*Eu)}
    if(!is.null(z)){E2 = exp(-rand*Ea-psiu*Eu-psiz*Ez)}
    ev2 = diag(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    acc.prob = min(-0.5*(log.det2 - log.det) - 0.5*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(acc.prob))
    if(acc==1){psia = rand}
    PSIA[g] = psia
    
    # update psiu
    if(is.null(z)){EE = exp(-psia*Ea-psiu*Eu)}
    if(!is.null(z)){EE = exp(-psia*Ea-psiu*Eu-psiz*Ez)}
    ev1 = diag(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))^2
    log.det = sum(log(ev1))
    quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    rand = exp(rnorm(1,log(psiu),prior$step))
    if(is.null(z)){E2 = exp(-psia*Ea-rand*Eu)}
    if(!is.null(z)){E2 = exp(-psia*Ea-rand*Eu-psiz*Ez)}
    ev2 = diag(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    acc.prob = min(-0.5*(log.det2 - log.det) - 0.5*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(acc.prob))
    if(acc==1){psiu = rand}
    PSIU[g] = psiu
    
    if(!is.null(z)){
      # update psiz
      EE = exp(-psia*Ea-psiu*Eu-psiz*Ez)
      ev1 = diag(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))^2
      log.det = sum(log(ev1))
      quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
      rand = exp(rnorm(1,log(psiz),prior$step))
      E2 = exp(-psia*Ea-psiu*Eu-rand*Ez)
      ev2 = diag(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))^2
      log.det2 = sum(log(ev2))
      quad2 = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
      acc.prob = min(-0.5*(log.det2 - log.det) - 0.5*(quad2 - quad), 0)
      acc=rbinom(1,1,exp(acc.prob))
      if(acc==1){psiz = rand}
      PSIZ[g] = psiz
    }
    
    # update tau
    a.tau = prior$a.t + 0.5*n
    b.tau = prior$b.t + 0.5*sum((c(Y)-c(phi))^2)
    tau = rgamma(1,shape=a.tau,rate=b.tau)
    TAU[g] = tau
    
    # calculate the current loglik
    if(is.null(z)){EE.tmp = exp(-psia*Ea-psiu*Eu)}
    if(!is.null(z)){EE.tmp = exp(-psia*Ea-psiu*Eu-psiz*Ez)}
    ev1 = diag(chol(psi1*EE.tmp + tau^(-1)*diag(rep(1,n))))^2
    log.det = sum(log(ev1))
    quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE.tmp + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    loglik.tmp = -0.5*log.det - 0.5*quad
    
    # update node selection
    for(k in 1:p){
      # Between models
      gamma.tmp = gamma
      rho.tmp = rho
      beta.tmp = beta
      if(gamma[k]==1){
        gamma.tmp[k] = 0
        rho.tmp[k] = 1
        beta.tmp[k] = 0
        pro.ratio = log((1-pro)/pro)
      }else if(gamma[k]==0){
        gamma.tmp[k] = 1
        rho.tmp[k] = runif(1)
        beta.tmp[k] = -log(rho.tmp[k])
        pro.ratio = log(pro/(1-pro))
      }
      U.beta = lapply(U,function(x) hadamard.prod(x,matrix(sqrt(beta.tmp),ncol=1)%*%matrix(rep(1,d),nrow=1)))
      Eu.tmp = matrix(0,n,n)
      for(i in 2:n){
        for(j in 1:(i-1)){
          Eu.tmp[i,j]=Eu.tmp[j,i]=norm(U.beta[[i]]-U.beta[[j]],"F")^2
        }
      }
      Eu.tmp = Eu.tmp/sd.u
      if(is.null(z)){EE.tmp = exp(-psia*Ea-psiu*Eu.tmp)}
      if(!is.null(z)){EE.tmp = exp(-psia*Ea-psiu*Eu.tmp-psiz*Ez)}
      ev1 = diag(chol(psi1*EE.tmp + tau^(-1)*diag(rep(1,n))))^2
      log.det = sum(log(ev1))
      quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE.tmp + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
      loglik.tmp2 = -0.5*log.det - 0.5*quad
      acc.prob = min(pro.ratio + loglik.tmp2 - loglik.tmp, 0)
      acc=rbinom(1,1,exp(acc.prob))
      if(acc==1){
        gamma = gamma.tmp
        rho = rho.tmp
        beta = beta.tmp
        loglik.tmp = loglik.tmp2
      }
      
      # Within model, perform only when accept gamma[k]=1
      if(acc==1 && gamma[k]==1){
        rho.tmp[k] = runif(1)
        beta.tmp[k] = -log(rho.tmp[k])
        U.beta = lapply(U,function(x) hadamard.prod(x,matrix(sqrt(beta.tmp),ncol=1)%*%matrix(rep(1,d),nrow=1)))
        Eu.tmp = matrix(0,n,n)
        for(i in 2:n){
          for(j in 1:(i-1)){
            Eu.tmp[i,j]=Eu.tmp[j,i]=norm(U.beta[[i]]-U.beta[[j]],"F")^2
          }
        }
        Eu.tmp = Eu.tmp/sd.u
        if(is.null(z)){EE.tmp = exp(-psia*Ea-psiu*Eu.tmp)}
        if(!is.null(z)){EE.tmp = exp(-psia*Ea-psiu*Eu.tmp-psiz*Ez)}
        ev1 = diag(chol(psi1*EE.tmp + tau^(-1)*diag(rep(1,n))))^2
        log.det = sum(log(ev1))
        quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE.tmp + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
        loglik.tmp2 = -0.5*log.det - 0.5*quad
        acc.prob = min(loglik.tmp2 - loglik.tmp, 0)
        acc=rbinom(1,1,exp(acc.prob))
        if(acc==1){
          rho = rho.tmp
          beta = beta.tmp
          loglik.tmp = loglik.tmp2
        }
      }
    }
    GAMMA[,g] = gamma
    BETA[,g] = beta
    
    # update pro
    pro = rbeta(1, shape1=(prior$a.p+sum(gamma==1)), shape2=(prior$b.p+sum(gamma==0)))
    
    # update Eu and Eu.full
    U.full.beta = lapply(U.full,function(x) hadamard.prod(x,matrix(sqrt(beta),ncol=1)%*%matrix(rep(1,d),nrow=1)))
    Eu.full = matrix(0,(n+nn),(n+nn))
    for(i in 2:(n+nn)){
      for(j in 1:(i-1)){
        Eu.full[i,j]=Eu.full[j,i]=norm(U.full.beta[[i]]-U.full.beta[[j]],"F")^2
      }
    }
    Eu.full = Eu.full/sd.u
    Eu = Eu.full[1:n,1:n]
    
    # update phi
    if(is.null(z)){KK = psi1*exp(-psia*Ea-psiu*Eu) + diag(rep(1e-12,n))}
    if(!is.null(z)){KK = psi1*exp(-psia*Ea-psiu*Eu-psiz*Ez) + diag(rep(1e-12,n))}
    inv.KK = chol2inv(chol(KK))
    v.phi = chol2inv(chol(inv.KK + tau*diag(rep(1,n))))
    m.phi = v.phi%*%(tau*matrix(Y,ncol=1))
    phi = rmvnorm(1,mean=m.phi,sigma=v.phi)
    PHI[,g] = phi
    
    # log likelihood
    K = KK + tau^(-1)*diag(rep(1,n))
    inv.K = chol2inv(chol(K))
    LOG[g] = -0.5*log(det(K)) - 0.5*t(Y)%*%inv.K%*%Y
    
    # prediction
    if(is.null(z)){cov.K = psi1*exp(-psia*Ea.full[1:n,(n+1):(n+nn)]-psiu*Eu.full[1:n,(n+1):(n+nn)])}
    if(!is.null(z)){cov.K = psi1*exp(-psia*Ea.full[1:n,(n+1):(n+nn)]-psiu*Eu.full[1:n,(n+1):(n+nn)]-psiz*Ez.full[1:n,(n+1):(n+nn)])}
    # pred = t(cov.K)%*%inv.K%*%matrix(Y,ncol=1)
    pred = t(cov.K)%*%inv.KK%*%matrix(phi,ncol=1)
    PRED[,g] = pred
    setTxtProgressBar(pb, g)
    if(g%%100==0){
      if(is.null(z)){
        message(cat("iteration:",g,", selected:",sum(gamma),", tau:",tau,", psi1:",psi1,
                    ", psia:",psia,", psiu:",psiu,", loglik:",LOG[g],"\n"))
      }
      if(!is.null(z)){
        message(cat("iteration:",g,", selected:",sum(gamma),", tau:",tau,", psi1:",psi1,
                    ", psia:",psia,", psiu:",psiu,", psiz:",psiz,", loglik:",LOG[g],"\n"))
      }
    }
  }
  close(pb)
  post.tau=mean(TAU[(burn+1):nrun])
  post.psi1=mean(PSI1[(burn+1):nrun])
  post.psia=mean(PSIA[(burn+1):nrun])
  post.psiu=mean(PSIU[(burn+1):nrun])
  if(!is.null(z)){post.psiz=mean(PSIZ[(burn+1):nrun])}
  if(is.null(z)){
    post=c(post.tau,post.psi1,post.psia,post.psiu,sd.a,sd.u)
    return(list(burn=burn, nrun=nrun, post=post, post.loglik=LOG, post.mean=PHI,
                post.pred=PRED, GAMMA=GAMMA, BETA=BETA, TAU=TAU, PSI1=PSI1,
                PSIA=PSIA, PSIU=PSIU))
  }
  if(!is.null(z)){
    post=c(post.tau,post.psi1,post.psia,post.psiu,post.psiz,sd.a,sd.u,sd.z)
    return(list(burn=burn, nrun=nrun, post=post, post.loglik=LOG, post.mean=PHI,
                post.pred=PRED, GAMMA=GAMMA, BETA=BETA, TAU=TAU, PSI1=PSI1,
                PSIA=PSIA, PSIU=PSIU, PSIZ=PSIZ))
  }
}