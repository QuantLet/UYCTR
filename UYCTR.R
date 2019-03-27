# -------------------------------------------------------------------------------
# Unified yield curve Taylor rule estimation functions
# -------------------------------------------------------------------------------

# Function to Calculate Likelihood
Linn     = function(para){
  alpha1 = para[1]; alpha2 = para[2]; beta0 = para[3]
  sQ1    = para[4];  sQ2 = para[5];  like=0
  xf     = matrix(0, nstate, 1)  
  # x filter
  xp     = matrix(0, nstate, 1)  
  # x pred
  Pf     = diag(.1, nstate)     
  # filter cov
  Pp     = diag(.1, nstate)      
  # pred cov
  pi11   <- .75 -> pi22;  pi12 <- .25 -> pi21; pif1 <- .5 -> pif2
  phi      = matrix(0,nstate,nstate)
  phi[1,1] = alpha1; phi[1,2] = alpha2; phi[2,1]=1; phi[4,4]=1
  Ups      = as.matrix(rbind(0,0,beta0,0))
  Q        = matrix(0,nstate,nstate)
  Q[1,1]   = sQ1^2; Q[3,3] = sQ2^2; R=0  
  # R=0 in final model
  # begin filtering 
  for(i in 1:num){
    xp   = phi%*%xf + Ups; Pp = phi%*%Pf%*%t(phi) + Q
    sig1 = as.numeric(M1%*%Pp%*%t(M1) + R)
    sig2 = as.numeric(M2%*%Pp%*%t(M2) + R)
    k1   = Pp%*%t(M1)/sig1; k2 = Pp%*%t(M2)/sig2
    e1   = y[i]-M1%*%xp; e2 = y[i]-M2%*%xp
    pip1 = pif1*pi11 + pif2*pi21; pip2 = pif1*pi12 + pif2*pi22
    den1 = (1/sqrt(sig1))*exp(-.5*e1^2/sig1)
    den2 = (1/sqrt(sig2))*exp(-.5*e2^2/sig2)
    denm = pip1*den1 + pip2*den2
    pif1 = pip1*den1/denm; pif2 = pip2*den2/denm
    pif1 = as.numeric(pif1); pif2 = as.numeric(pif2)
    e1   = as.numeric(e1); e2=as.numeric(e2)
    xf   = xp + pif1*k1*e1 + pif2*k2*e2
    eye  = diag(1, nstate)
    Pf   = pif1*(eye-k1%*%M1)%*%Pp + pif2*(eye-k2%*%M2)%*%Pp
    like = like - log(pip1*den1 + pip2*den2)
    prob[i]<<-pip2; xfilter[,,i]<<-xf; innov.sig<<-c(sig1,sig2)
    yp[i]<<-ifelse(pip1 > pip2, M1%*%xp, M2%*%xp)  
  }
  return(like)   
}

# Estimation
alpha1       = 1.4; alpha2 = -.5; beta0 = .3; sQ1 = .1; sQ2 = .1
est          = function(alphal, alpha2, beta0, sQ1, sQ2){
  init.par   = c(alpha1, alpha2, beta0, sQ1, sQ2)
  (est       = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE, control=list(trace=1,REPORT=1)))
  SE         = sqrt(diag(solve(est$hessian)))
  u          = cbind(estimate=est$par, SE)
  rownames(u)=c('alpha1','alpha2','beta0','sQ1','sQ2'); 
  u
}


#   y(t) = x(t) + v(t);    v(t) ~ iid N(0,V)                     
#   x(t) = x(t-1) + w(t);  w(t) ~ iid N(0,W)                        
#   x(0) ~ N(m0,C0);  V ~ IG(a,b);  W ~ IG(c,d)
#   x(t|t) ~ N(m,C);  x(t|n) ~ N(mm,CC);  x(t|t+1) ~ N(a,R)  


ffbs = function(y,V,W,m0,C0){
  n  = length(y);  
  a  = rep(0,n);  
  R  = rep(0,n)
  m  = rep(0,n);   
  C  = rep(0,n);  
  B  = rep(0,n-1)     
  H  = rep(0,n-1); 
  mm = rep(0,n);  
  CC = rep(0,n)
  x  = rep(0,n); llike = 0.0
  for (t in 1:n){
    if(t==1){a[1] = m0; R[1] = C0 + W
    }else{ a[t]   = m[t-1]; R[t] = C[t-1] + W }
    f      = a[t]
    Q      = R[t] + V
    A      = R[t]/Q
    m[t]   = a[t]+A*(y[t]-f)
    C[t]   = R[t]-Q*A**2
    B[t-1] = C[t-1]/R[t]
    H[t-1] = C[t-1]-R[t]*B[t-1]**2
    llike  = llike + dnorm(y[t],f,sqrt(Q),log=TRUE) }
    mm[n]  = m[n]; CC[n] = C[n]
    x[n]   = rnorm(1,m[n],sqrt(C[n]))
  for (t in (n-1):1){
    mm[t]  = m[t] + C[t]/R[t+1]*(mm[t+1]-a[t+1])
    CC[t]  = C[t] - (C[t]^2)/(R[t+1]^2)*(R[t+1]-CC[t+1])
    x[t]   = rnorm(1,m[t]+B[t]*(x[t+1]-a[t+1]),sqrt(H[t]))  }
  return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike))   }
# Simulate states and data
  set.seed(1); 
  W    = 0.5; 
  V    = 1.0
  n    = 100; 
  m0   = 0.0; 
  C0   = 10.0; 
  x0   = 0
  w    = rnorm(n,0,sqrt(W))
  v    = rnorm(n,0,sqrt(V))
  x    = y = rep(0,n)
  x[1] = x0   + w[1]
  y[1] = x[1] + v[1]
 for (t in 2:n){
  x[t] = x[t-1] + w[t]
  y[t] = x[t] + v[t]   }
  run  = ffbs(y,V,W,m0,C0)
  m    = run$m 
  C    = run$C
  mm   = run$mm
  CC   = run$CC
  L    = m-2*C 
  U1   = m+2*C
  L2   = mm-2*CC
  U2   = mm+2*CC
  N    = 50
  Vs   = seq(0.1,2,length=N)
  Ws   = seq(0.1,2,length=N)
  likes  = matrix(0,N,N)
 for (i in 1:N){
  for (j in 1:N){
    V    = Vs[i]
    W    = Ws[j]
    run  = ffbs(y,V,W,m0,C0)    
    likes[i,j] = run$llike  } 
  }  

# Hyperparameters
a = 0.01; b = 0.01; c = 0.01; d = 0.01
Hyper = function(a, b, c, d){
  set.seed(90210)
  burn  = 10;  M = 1000
  niter = burn + M
  V1    = V;  W1 = W
  draws = NULL
  all_draws = NULL
  for (iter in 1:niter){
  run   = ffbs(y,V1,W1,m0,C0)
  x     = run$x
  V1    = 1/rgamma(1,a+n/2,b+sum((y-x)^2)/2)
  W1    = 1/rgamma(1,c+(n-1)/2,d+sum(diff(x)^2)/2)
  draws = rbind(draws,c(V1,W1,x))    }
  all_draws = draws[,1:2]
  q025  = function(x){quantile(x,0.025)}
  q975  = function(x){quantile(x,0.975)}
  draws = draws[(burn+1):(niter),]
  xs    = draws[,3:(n+2)]
  lx    = apply(xs,2,q025)
  mx    = apply(xs,2,mean)
  ux    = apply(xs,2,q975)
  set.seed(90210)
  burnin  = 100 
  step    = 10   
  M       = 1000  
  niter   = burnin+step*M
  pars    = matrix(0,niter,4)
  xbs     = array(0,c(niter,n,4))
  pr <- progress_text()           
  # displays progress
  pr$init(niter)                  
  for (iter in 1:niter){
    xb = ffbs(y,F,G,V,W11,W22,a1,R1)
    u  = xb[,1] 
    yu = diff(u); xu = u[-n]    
    # for phihat and se(phihat)
    regu = lm(yu~0+xu)                
    # est of beta = phi-1
    phies  = as.vector(coef(summary(regu)))[1:2] + c(1,0) 
    # phi estimate and SE
    dft    = df.residual(regu)   
    G[1,1] = phies[1] + rt(1,dft)*phies[2] 
    # use a t
    V      = 1/rgamma(1, (n0+n)/2, (n0*s20v/2) + sum((y-xb[,1]-xb[,2])^2)/2)
    W11    = 1/rgamma(1, (n0+n-1)/2, (n0*s20w/2) + sum((xb[-1,1]-phies[1]*xb[-n,1])^2)/2)
    W22    = 1/rgamma(1, (n0+ n-3)/2, (n0*s20w/2) + sum((xb[4:n,2] + xb[3:(n-1),2] + 
                                                        xb[2:(n-2),2] +xb[1:(n-3),2])^2)/2)
    xbs[iter,,] = xb
    pars[iter,] = c(G[1,1], sqrt(V), sqrt(W11), sqrt(W22))
    pr$step()             
  }
  return(pr$step())
}

L.mle  = function(theta, sigma, y, x){
  
  if(length(theta)>1){
    m.s      = length(x[, 1])
    x.2      = cbind(1,x)
    A.1      = (-1)*m.s*log(sigma)
    B.1      = (-1)*(1/(2*sigma^2))
    D.1      = t(y - x.2 %*% theta ) %*% (y - x.2 %*% theta)
    E        = A.1+B.1*D.1
  }
  
  if(length(theta)==1){
    m.s      = length(x)
    x.2      = cbind(1,x)
    A.1      = (-1)*m.s*log(sigma)
    B.1      = (-1)*(1/(2*sigma^2))
    D.1      = t(y - x.2 * theta ) %*% (y - x.2 * theta)
    E        = A.1+B.1*D.1
  }
  E
  
}
