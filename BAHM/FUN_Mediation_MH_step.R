
# Only needed for BM model (model=2)

FUN_Mediation_MH_step= function(Model,Data,Prior,Mcmc)
{
  myunireg = function(y,X,betabar,A,nu,q)  # Generate beta and sigma2
  { n=length(y)
    nvar = ncol(X)
    RA = chol(A)
    W = rbind(X, RA)
    z = c(y, as.vector(RA %*% betabar))
    IR = backsolve(chol(crossprod(W)), diag(nvar))
    btilde = crossprod(t(IR)) %*% crossprod(W, z)
    res = z - W %*% btilde   # Traditional
    s = t(res) %*% res
    done=0
    while (done==0)
    { sigma2 = ((nu * q + s)/rchisq(1, nu + n))
    #print(nu)
    #print(q)
    #print(s)
    #print(n)
      done = c(ifelse(sigma2>0.01,1,0))
    }
    beta = btilde + as.vector(sqrt(sigma2)) * IR %*% rnorm(nvar)
    return(list(beta=beta,sigma2=sigma2))
  }
  
  # Inputs
  ModelFlag = Model  # 1 = aggregate, 2 = BM
  nvarX=ncol(Data$X) # including intercept
  nvarM=ncol(Data$m) # no intercept
  nvarZ=ncol(Data$Z) # includes intercept, it is just intercept in the absence of covariates
  m=Data$m
  X=Data$X
  y=Data$y
  Z=Data$Z
  nobs=length(y)
  dima = ncol(X)   # dim of alpha
  dimg = ncol(X)   # dim of gamma
  dimgb = dimg +1
  dimz = ncol(Z)   # dim of alpha
  R=Mcmc$R
  keep=Mcmc$keep
  slambda1=Mcmc$slambda1
	slambda2=Mcmc$slambda2

	#Generating X's for the cases with covariates for the aggregate model
	#DOES NOT WORK NOW
	# Make user generate the interaction variables for the aggregate model is needed before uploading the file
	
  if (missing(Prior)) {
    ma = c(rep(0, dima))     # mean of alpha
    Aa = 0.01 * diag(dima)
    mgb = c(rep(0, (dimgb)))   # mean of betagamma
    Agb = 0.01 * diag(dimgb)
    ml = c(rep(0, dimz))     # mean of alpha
    Al = 0.01 * diag(dimz)
    nu = 5
    qy= c(var(y),var(y))
    qm= c(var(m),var(m))
    #g=3          # this is nu and q for the prior distribution for rho in (52)
  } else {
      ma = Prior$ma
      Aa = Prior$Aa
      mgb = Prior$mgb
      Agb = Prior$Agb
    nu = Prior$nu
    qy = Prior$qy   # vector (M,G)
    qm = Prior$qm   # vector (M,G)
    ml =Prior$ml
	  Al = Prior$Al
    #g=Prior$g      # this is nu and q for the prior distribution for rho in (52)
  }
	
  
  if (is.null(qy)==T) {  
      qy= c(var(y),var(y))  
  } else {  
      qy = c(Prior$qy, Prior$qy)
  }
	
  if (is.null(qm)==T) {  
      qm= c(var(m),var(m))  
  } else {  
      qm = c(Prior$qm, Prior$qm)
  }
  
  set.seed(Mcmc$seed)
  
  # Set up storage
  alphadraw = array(0, dim = c(floor(R/keep), dima, 2))        # 2 classes
  betaMdraw = array(0, dim = c(floor(R/keep), 2)) 
  gammabetaSdraw = array(0, dim = c(floor(R/keep), dimgb))   
  lambdadraw = array(0,dim=c(floor(R/keep),dimz))
  sigma2mdraw = array(0, dim = c(floor(R/keep), 2))
  sigma2ydraw = array(0, dim = c(floor(R/keep), 2))
  rhodraw = matrix(0,ncol=nobs,nrow=floor(R/keep))
  wdraw=matrix(0,ncol=floor(R/keep),nrow=nobs)
  LL_total   = c(rep(0,floor(R/keep)))
  LL = array(0,dim=c(4,nobs,floor(R/keep)))
  temp=c(rep(0,nobs))
  rej=0
	reject=c(rep(0,floor(R/keep)))
  
  
 	# Set up initial values for BM
  if(ModelFlag==2) { rho=c(rep(0.5,nobs))    
	                   w = c(rbinom(nobs,1,rho))
  }
	lambda = matrix(0,ncol=1,nrow=nvarZ)
  
	#  Mediating
  alphaM = matrix(0,ncol=1,nrow=nvarX)
  betaM = matrix(0,ncol=1,nrow=nvarM+1)  
  sigma2yM = 1
  sigma2mM = 1
  
  #  Standard
  alphaD = matrix(0,ncol=1,nrow=nvarX)
  gammabetaD = matrix(0,ncol=1,nrow=nvarM+nvarX)
  sigma2yD = 1
  sigma2mD = 1
  
  itime = proc.time()[3]
  cat("MCMC Iteration (est time to end - min)", fill = TRUE)
 
  #-- adding MH step auto-adjustment for lambda block
  
  done=0 
  while (done==0)
  {
     for(r in 1:R)
     {
       ###  Segment 1 (called M here)
       indexM = c(which(w==1))
       yM = y[indexM]
       mM = m[indexM]
       XM = matrix(X[indexM,],ncol=nvarX)
       if(1){  out_Xtom = myunireg(y=mM,X=XM,betabar=ma,A=Aa,nu=nu,q=qm[1])
              alphaM = out_Xtom$beta
              sigma2mM = out_Xtom$sigma2
       }
       if(1) { out_mtoyM = myunireg(y=yM,X=cbind(rep(1,length(mM)),mM),betabar=mgb[1:2],A=Agb[1,1]*diag(2),nu=nu,q=qy[1])   # one mediator
              betaM = out_mtoyM$beta
              sigma2yM = out_mtoyM$sigma2
       }
       ###  Segment 2 (called D here, but is S=standard segment)
       indexD = c(which(w==0))
       yD = y[indexD]
       mD = m[indexD]
       XD = matrix(X[indexD,],ncol=nvarX)
       if(1){  out_Xtom = myunireg(y=mD,X=XD,betabar=ma,A=Aa,nu=nu,q=qm[2])
               alphaD = out_Xtom$beta
               sigma2mD = out_Xtom$sigma2
       }
       if(1) { out_mtoyD = myunireg(y=yD,X=cbind(XD,mD),betabar=mgb,A=Agb,nu=nu,q=qy[2])
               gammabetaD=out_mtoyD$beta
               sigma2yD = out_mtoyD$sigma2
       }
   
       ####  Updating segment probability / Lambda
       slambda = slambda1
	     lambdanew = c(lambda + slambda * rnorm(nvarZ))
	     
	     rhonew = exp(Z %*% (lambdanew))/(1+exp(Z %*% (lambdanew)))
     	 lognew_lambda = sum(w * log(c(rhonew)) + (1-w)* log(1-c(rhonew)))
	     logold_lambda = sum(w * log(rho) + (1-w)* log(1-rho))
       logknew = -0.5 * (t(c(lambdanew) - ml) %*% Al %*% (c(lambdanew) - ml))
       logkold = -0.5 * (t(c(lambda) - ml) %*% Al %*% (c(lambda) - ml))
     	   alpha = exp(lognew_lambda + logknew - logold_lambda - logkold )
          alpha=ifelse(is.finite(alpha)==T,alpha,-1)
          u = runif(n = 1, min = 0, max = 1)
          if (u < alpha) {  lambda = lambdanew
	                          rho=c(rhonew)               }
          else           {  rej = rej + 1            }			
       
       ####  Updating segment indicators (working with likelihood, not log-lik here)
       LLikD = exp(-0.5*(log(c(2*pi*sigma2yD))+log(c(2*pi*sigma2mD))) - 0.5*(m-X%*%(alphaD))^2/ c(sigma2mD) - 0.5*(y-cbind(X,m)%*%gammabetaD)^2 / c(sigma2yD))
       if(ModelFlag==2)
       { LLikM = exp(-0.5*(log(c(2*pi*sigma2yM))+log(c(2*pi*sigma2mM))) - 0.5*(m-X%*%(alphaM))^2/ c(sigma2mM) - 0.5*(y-cbind(X[,1],m)%*%betaM)^2 / c(sigma2yM)) 
       }  
       temp = rho*LLikM/ (rho*LLikM+(1-rho)*LLikD)   #  pi in documents
       for(hh in 1:nobs)
          { w[hh] = rbinom(1,1,temp[hh]) 
       }
     }
     
     if( rej/R < 0.66 ) {  
       slambda1 = slambda * runif(1,1,1.25)
       done = 0
       rej=0 }
     if( rej/R > 0.74 ) {  
       slambda1 = slambda * runif(1,0.75,1)
       done = 0
       rej=0 }
     else {  
       slambda1 = slambda
       done = 1 }
   }
  #-- end of MH step auto-adjustment
  
  ctime = proc.time()[3]
  #cat(" Total Time Elapsed: ", round((ctime - itime)/60, 2),fill = TRUE)
  return(list(slambda=slambda, RR = rej/R))
}
