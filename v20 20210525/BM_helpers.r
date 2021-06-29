
#------------ prep model inputs --------------#

get_inputs_binary <-function(df, x_vars, y_var, m_var, 
                             Aa_var, Abg_var, nu_var, qy_var, qm_var,
                             R_var, keep_var, g_var)
  {
    df  = as.matrix(df)
    X   = df[,x_vars]
    nhh = nrow(df)
    y   = df[,y_var]
    m   = df[,m_var]
    nvarX = 1 + length(x_vars) # number of X variables selected on the menu. intercept is added here
    nvarM = 1
    
    Data = list(X = cbind(rep(1, nhh), X), 
                y = y, 
                m = matrix(m, ncol=1)
    )
    Prior = list(ma  = c(rep(0, nvarX)),
                 Aa  = Aa_var * diag(nvarX),
                 mgb = c(rep(0, nvarX+1)), 
                 Agb = Abg_var * diag(nvarX+1),
                 nu  = nu_var,
                 qy  = c(qy_var, qy_var),
                 qm  = c(qm_var, qm_var),
                 g   = g_var
    )
    Mcmc = list(Rep=R_var, keep=keep_var)
    
    return(list(Data=Data, Prior=Prior, Mcmc=Mcmc))
}

update_inputs_binary <- function(inputs_binary, seed_var) {
  inputs_binary$Mcmc$seed <- seed_var
  
  return(inputs_binary)
}


#-------------------- run model: 1 loop ------------------#

FUN_Mediation_LCRM_2class_MS_Gibbs_forShinyApp = function(Data, Prior, Mcmc, p = NULL)
{
  myunireg = function(y, X, betabar, A, nu, q)  # Generate beta and sigma2
  { 
    n    = length(y)
    
    nvar = ncol(X)
    RA   = chol(A)
    W    = rbind(X, RA)
    z    = c(y, as.vector(RA %*% betabar))
    IR   = backsolve(chol(crossprod(W)), diag(nvar))
    btilde = crossprod(t(IR)) %*% crossprod(W, z)
    res    = z - W %*% btilde   # Traditional
    s     = t(res) %*% res
    done = 0
    while (done == 0) {
      sigma2 = ((nu * q + s)/rchisq(1, nu + n))
      done = c(ifelse(sigma2>0.01, 1, 0))
    }
    beta = btilde + as.vector(sqrt(sigma2)) * IR %*% rnorm(nvar)
    return(list(beta=beta, sigma2=sigma2))
  }
  
  # Inputs
  nvarX = ncol(Data$X) # including intercept
  nvarM = ncol(Data$m) # no intercept
  m     = Data$m
  X     = Data$X
  y     = Data$y
  nobs  = length(y)
  dima  = ncol(X)   # dim of alpha
  dimg  = ncol(X)   # dim of gamma
  dimgb = dimg +1
  R     = Mcmc$Rep
  keep  = Mcmc$keep
  
  if (missing(Prior)) {
    ma  = c(rep(0, dima))     # mean of alpha
    Aa  = 0.01 * diag(dima)
    mgb = c(rep(0, (dimgb)))   # mean of betagamma
    Agb = 0.01 * diag(dimgb)
    nu  = 5
    qy  = c(var(y),var(y))
    qm  = c(var(m),var(m))
    g   = 3          # this is nu and q for the prior distribution for rho in (52)
  } else {
    ma  = Prior$ma
    Aa  = Prior$Aa
    mgb = Prior$mgb
    Agb = Prior$Agb
    nu  = Prior$nu
    qy  = Prior$qy   # vector (M,G)
    qm  = Prior$qm   # vector (M,G)
    g   = Prior$g    # this is nu and q for the prior distribution for rho in (52)
  }
  
  if (is.null(qy)==T) {  qy = c(var(y),var(y))  }
  else                {  qy = Prior$qy }
  if (is.null(qm)==T) {  qm = c(var(m),var(m))  }
  else                {  qm = Prior$qm }
  
  set.seed(Mcmc$seed)
  
  # Set up storage
  alphadraw = array(0, dim = c(floor(R/keep), dima, 2))        # 2 classes
  betaMdraw = array(0, dim = c(floor(R/keep), 2)) 
  gammabetaSdraw = array(0, dim = c(floor(R/keep), dimgb))   
  sigma2mdraw    = array(0, dim = c(floor(R/keep), 2))
  sigma2ydraw    = array(0, dim = c(floor(R/keep), 2))
  #LLdraw      = matrix(double(floor(R/keep) * nobs), ncol = nobs)
  #LLtotaldraw = c(0,floor(R/keep))
  rhodraw = matrix(0,ncol=1,nrow=floor(R/keep))
  wdraw   = matrix(0,ncol=floor(R/keep),nrow=nobs)
  LL_total = c(rep(0,floor(R/keep)))
  LL       = array(0,dim=c(4,nobs,floor(R/keep)))
  temp     = c(rep(0,nobs))
  
  # Set up initial values
  rho = 0.5           
  w   = c(rbinom(nobs,1,rho))
  #  Mediating
  alphaM = matrix(0,ncol=1,nrow=nvarX)
  betaM  = matrix(0,ncol=1,nrow=nvarM+1)  
  sigma2yM = 1
  sigma2mM = 1
  #  Standard
  alphaD     = matrix(0,ncol=1,nrow=nvarX)
  gammabetaD = matrix(0,ncol=1,nrow=nvarM+nvarX)
  sigma2yD = 1
  sigma2mD = 1
  
  itime = proc.time()[3]
  cat("MCMC Iteration (est time to end - min)", fill = TRUE)
  
  # begin loop
  for(r in 1:R)
  {
    ###  Segment 1 (called M here)
    indexM = c(which(w==1))
    yM = y[indexM]
    mM = m[indexM]
    XM = matrix(X[indexM,],ncol=nvarX)
    if(1){
      out_Xtom = myunireg(y=mM,X=XM,betabar=ma,A=Aa,nu=nu,q=qm[1])
      alphaM   = out_Xtom$beta
      sigma2mM = out_Xtom$sigma2
    }
    if(1){ 
      out_mtoyM = myunireg(y=yM, X=cbind(rep(1,length(mM)), mM), betabar=mgb[1:2], A=Agb[1,1]*diag(2), nu=nu, q=qy[1])   # one mediator
      betaM = out_mtoyM$beta
      sigma2yM = out_mtoyM$sigma2
    }
    
    ###  Segment 2 (called D here, but is S=standard segment)
    indexD = c(which(w==0))
    yD = y[indexD]
    mD = m[indexD]
    XD = matrix(X[indexD,],ncol=nvarX)
    if(1){
      out_Xtom = myunireg(y=mD,X=XD,betabar=ma,A=Aa,nu=nu,q=qm[2])
      alphaD = out_Xtom$beta
      sigma2mD = out_Xtom$sigma2
    }
    
    if(1){ 
      out_mtoyD = myunireg(y=yD,X=cbind(XD,mD),betabar=mgb,A=Agb,nu=nu,q=qy[2])
      gammabetaD=out_mtoyD$beta
      sigma2yD = out_mtoyD$sigma2
    }
    
    ####  Updating segment indicators
    if(1){ 
      LLikM = exp(-0.5*(log(c(2*pi*sigma2yM))+log(c(2*pi*sigma2mM))) - 0.5*(m-X%*%(alphaM))^2/ c(sigma2mM) - 0.5*(y-cbind(X[,1],m)%*%betaM)^2 / c(sigma2yM)) 
      LLikD = exp(-0.5*(log(c(2*pi*sigma2yD))+log(c(2*pi*sigma2mD))) - 0.5*(m-X%*%(alphaD))^2/ c(sigma2mD) - 0.5*(y-cbind(X,m)%*%gammabetaD)^2 / c(sigma2yD)) 
      temp = rho*LLikM/ (rho*LLikM+(1-rho)*LLikD)   #  pi in documents
      for(hh in 1:nobs){
        temppi = ifelse(is.finite(temp[hh])==T,ifelse(temp[hh]>1,1,temp[hh]),0)    # pi in documents
        w[hh] = rbinom(1,1,temppi)
      }
    }
    
    ####  Updating segment probability
    if(1){ 
      rho = rbeta(1,g+sum(w),g+nobs-sum(w))
    }
    
    ####   Save draws
    if (r%%keep == 0){
      mkeep = r/keep
      alphadraw[mkeep,,1] = alphaM
      alphadraw[mkeep,,2] = alphaD
      betaMdraw[mkeep,] = betaM
      gammabetaSdraw[mkeep,] = gammabetaD
      sigma2mdraw[mkeep,] = c(sigma2mM,sigma2mD)
      sigma2ydraw[mkeep,] = c(sigma2yM,sigma2yD)
      rhodraw[mkeep] = rho        # used to be phi
      wdraw[,mkeep] = w
      
      LL_total[mkeep]= sum(log(LLikM[indexM])) + sum(log(LLikD[indexD]))
      LL[1:2,,mkeep] = cbind(LLikM,LLikD)
    }
    
    mkeep = r/keep
    
    if ((r%%100) == 0){
      ctime = proc.time()[3]
      timetoend = ((ctime - itime)/r) * (R - r)
      cat(" ", r, "(", round(timetoend/60, 1), ")",fill=T)
      
      if (!is.null(p)) p()
    }
  }
  
  ctime = proc.time()[3]
  
  cat(" Total Time Elapsed: ", round((ctime - itime)/60, 2),fill = TRUE)
  
  return(list(alphadraw=alphadraw,betaMdraw=betaMdraw,gammabetaSdraw=gammabetaSdraw,sigma2mdraw=sigma2mdraw,
              sigma2ydraw=sigma2ydraw,rhodraw=rhodraw,wdraw=wdraw,LL=LL,LL_total=LL_total,
              Aa=Aa,ma=ma, mgb=mgb,Agb=Agb,nu=nu,qy=qy,qm=qm,g=g,R=R,keep=keep,seed=Mcmc$seed))
}

# # test input and model run code
# #setwd('/Users/kgedney/Dropbox/Shiny APP dev/R Code/Shiny Code/v4')
# testdf <- read.csv('sample_data_Loyalty.csv')
# x_vars <- 'x1'
# y_var  <- 'y1'
# m_var  <- 'avg_m'
# Aa_var  <- 0.01
# Abg_var <- 0.01
# nu_var  <- 5
# qy_var  <- var(testdf[,y_var])
# qm_var  <- var(testdf[,m_var])
# R_var   <- 1000
# seed_var <- 123
# keep_var <- 10
# g_var <- 3
# 
# input_list <- get_inputs_binary(testdf, x_vars, y_var, m_var,
#                                 Aa_var, Abg_var, nu_var, qy_var, qm_var,
#                                 R_var, seed_var, keep_var, g_var)
# 
# output_list_b <- FUN_Mediation_LCRM_2class_MS_Gibbs_forShinyApp(input_list$Data,
#                                                                 input_list$Prior,
#                                                                 input_list$Mcmc)
# output_list_b = list(output_list_b,output_list_b)



#-------------------- run model: 20 loops unparallelized and get outputBM ------------------#


# --- Generate list of seeds ----------------------------------------------------
FUN_Mediation_SeedList_ForShiny  = function(seednum_var, seed_var)
  { set.seed(seed_var)
    seed.list <- c(sort(round(runif(seednum_var, 0, 10000))))
    #seed.index <- seq(1, seednum_var, 1)
    #num_seeds <- length(seed.index)  # should be the same as input$seednum_varalize list
    return(seed.list)
}
  
  

#---- (Mixture MS)  MCMC plots  MULTIX ----------------------------------------------------------------------------
# note: this will just be available as pdf download (no need to print to screen)
# this can be run without the RHat results
# title: MCMC plots
#  NOT DONE???

FUN_PDF_MCMC_Mediation_forShiny = function(dataset,filenamelist,seed.index,seed.list,burnin)
{  #pdf(paste(dataset,"MCMC for BM model", ".pdf", sep = ""), width=20, height=10)
   #plot.new()
   par(oma=c(0,0,2,0));   par(mfcol=c(8,length(seed.index)),mai = c(0.5, 0.3, 0.4, 0.3))
   filename = filenamelist[[1]]
   ylimAm=ylimAs = c( min(filename$alphadraw[-1:-burnin,,])-0.1*min(filename$alphadraw[-1:-burnin,,]),
                      max(filename$alphadraw[-1:-burnin,,])+0.1*max(filename$alphadraw[-1:-burnin,,]))
   ylimGBs = c( min(filename$gammabetaSdraw[-1:-burnin,])-0.1*min(filename$gammabetaSdraw[-1:-burnin,]),
                max(filename$gammabetaSdraw[-1:-burnin,])+0.1*max(filename$gammabetaSdraw[-1:-burnin,]))
   ylimBm =  c( min(filename$betaMdraw[-1:-burnin,])-0.1*min(filename$betaMdraw[-1:-burnin,]),
                max(filename$betaMdraw[-1:-burnin,])+0.1*max(filename$betaMdraw[-1:-burnin,]))
   ysigma =  c( min(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,]))-
                  0.1*min(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,])),
                max(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,]))+
                  0.1*max(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,])))
   ylimLL =  c( min(filename$LL_total[-1:-2])-0.1*min(filename$LL_total[-1:-2]),  # used burnin before
                max(filename$LL_total[-1:-2])+0.1*max(filename$LL_total[-1:-2]))
   for( i in seed.index)
   { filename = filenamelist[[i]]
     matplot(filename$alphadraw[,,1],type='l',col=1:8,main=expression(alpha[M]),ylab=expression(alpha[M]),ylim=ylimAm);
     matplot(filename$alphadraw[,,2],type='l',col=1:8,main=expression(alpha[S]),ylab=expression(alpha[S]),ylim=ylimAs);
     matplot(filename$betaMdraw,type='l',col=1:3,main=expression(beta[M]),ylab=expression(beta[M]),ylim=ylimBm);
     matplot(filename$gammabetaSdraw,type='l',col=1:8,main=expression(gamma[S],beta[S]),ylab=expression(gamma[S],beta[S]),ylim=ylimGBs);
     matplot(filename$sigma2mdraw,type='l',col=1:4,main=expression(sigma[m]^2),ylab=expression(sigma[m]^2),ylim=ysigma);
     matplot(filename$sigma2ydraw,type='l',col=1:4,main=expression(sigma[y]^2),ylab=expression(sigma[y]^2),ylim=ysigma);
     matplot(filename$rhodraw,type='l',col=1:3,main=expression(rho),ylab=expression(rho),ylim=c(0,1));
     plot(filename$LL_total,type='l',col=1,main="LL",ylab="LL",ylim=ylimLL)
   }
   #title(main=paste(dataset),outer=T)
   #pdf(paste("plot", j, ".pdf", sep = ""),width=20, height=10)
  # dev.off()
}



#---- (Mixture MS) Proportion of joint (alpha, beta) in each quadrant for both M and S segments MULTIX ---------------------------------------------------------------------------------
# this ONLY runs for ONE BEST seed, which will be taken from the results of RHat function

FUN_PDF_Mediation_AlphaBetaProportion_MSmixture_forShiny = function(filenamelist,seed.list,x_vars,burnin)
  # seed.list is the list of selected seeds for analysis, not all seeds
{  nvarX = ncol(filenamelist[[1]]$alphadraw)
   QuadrantsCountsM = matrix(0,nrow=4, ncol=nvarX-1)
   QuadrantsCountsS = matrix(0,nrow=4, ncol=nvarX-1)
   DrawsAnalysis = c(seq(burnin+1,length(filenamelist[[seed.list]]$LL_total),1))
   for(j in 2:nvarX)
   {  i = seed.list
      #for(i in seed.list)
      {  ppM=pnM=npM=nnM=0
         ppS=pnS=npS=nnS=0
         for(r in DrawsAnalysis)
         { ppM = ppM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]>0,ifelse(filenamelist[[i]]$betaMdraw[r,2]>0,1,0),0)
           pnM = pnM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]>0,ifelse(filenamelist[[i]]$betaMdraw[r,2]<0,1,0),0)
           npM = npM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]<0,ifelse(filenamelist[[i]]$betaMdraw[r,2]>0,1,0),0)
           nnM = nnM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]<0,ifelse(filenamelist[[i]]$betaMdraw[r,2]<0,1,0),0)
           ppS = ppS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]>0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
           pnS = pnS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]>0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
           npS = npS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]<0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
           nnS = nnS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]<0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
         }
         QuadrantsCountsM[,(j-1)]=c(ppM,pnM,npM,nnM)
         QuadrantsCountsS[,(j-1)]=c(ppS,pnS,npS,nnS)
       }
    }
  ProportionsM <- round(QuadrantsCountsM/length(DrawsAnalysis),4)
  ProportionsS <- round(QuadrantsCountsS/length(DrawsAnalysis),4)
  Proportions <- cbind(ProportionsM,ProportionsS)
  
  temp <- c(rep(0,(nvarX-1)*2))
  for(n in 2:nvarX)
  {
     temp[n-1] <-paste(x_vars[n - 1], " Segment M (mediating)",sep = "")
     temp[(nvarX-1)+n-1] <- paste(x_vars[n - 1], " Segment G (general)",sep = "")
  }
  colnames(Proportions) <- temp
  rownames(Proportions) <- c("I (++)","II (+-)","III (-+)","IV (--)")
  return(Proportions)
}

# testing
#temp_test = FUN_PDF_Mediation_AlphaBetaProportion_MSmixture_forShiny (output_list_b,seed.list=c(1),x_vars,burnin=10)

#---- (MS mixture model) Summary plots for selected solutions/runs MULTIX ---------------------------------------------------------------------------------
# plots for two seeds (again based on the RHat results)
# these will be printed to Shiny but also available as a pdf download
# will be changed to drop all the MCMC details


#-------------- BM model scatterplots (original version)--------------------

FUN_PDF_Mediation_FinalPlots_MSmixture_forShiny_Plot = function(dataset,filenamelist,seed.list,seed.selected,burnin){
  nvarX = ncol(filenamelist[[1]]$alphadraw[,,1])
  #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
  par(oma=c(0,0,2,0));   par(mfcol=c(3*(nvarX-1)+1,length(seed.selected)))
  filename = filenamelist[[seed.selected[1]]]
  for( i in seed.selected)
    { filename = filenamelist[[i]]
      ylimA = c( min(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2])-
                    0.1*min(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2]),
                 max(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2])+
                 0.1*max(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2])
              )
      ylimG = c( min(filename$gammabetaSdraw[-1:-burnin,-1])-0.1*min(filename$gammabetaSdraw[-1:-burnin,-1]),
                 max(filename$gammabetaSdraw[-1:-burnin,-1])+0.1*max(filename$gammabetaSdraw[-1:-burnin,-1]))
      ylimB =  c( min(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)])-
                    0.1*min(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)]),
                  max(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)])+
                    0.1*max(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)])
                )
      breaksCalc = (ylimG[2]-ylimG[1])/25
      for(p in 2:nvarX) {
        plot(filename$alphadraw[-1:-burnin,p,1], filename$betaMdraw[-1:-burnin,2], 
             main=paste("Segment M. Seed=",seed.list[i]),
           xlab=bquote(alpha[M][.(p-1)]),ylab=expression(beta[M]), xlim=ylimA,ylim=ylimB)
        abline(h=0,v=0,col="gray")
      }
      for(p in 2:nvarX) {   
        plot(filename$alphadraw[-1:-burnin,p,2], filename$gammabetaSdraw[-1:-burnin,nvarX+1], main=paste("Segment G. Seed=",seed.list[i]),
           xlab=bquote(alpha[S][.(p-1)]),ylab=expression(beta[S]), xlim=ylimA,ylim=ylimB)
        abline(h=0,v=0,col="gray")
      }
      for(p in 2:nvarX){
         hist(filename$gammabetaSdraw[-1:-burnin,p], xlab=bquote(gamma[S][.(p-1)]), xlim=ylimG, main=paste("Segment G. Seed=",seed.list[i]),
           breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,p]),(max(filename$gammabetaSdraw[-1:-burnin,p])+breaksCalc),breaksCalc )))
         abline(v=0,col="darkred",lwd=2)
      }
      hist(filename$rhodraw[-1:-burnin], main=paste("Rho. Seed=",seed.list[i]), xlab=bquote(rho), xlim=c(0,1),
           breaks=50)
    }
  #title(main=paste(dataset," Final Results. Binary mixture"),outer=T)
  #dev.off()
}

#-------------- BM model scatterplots (updatet version for the App)--------------------
# for one seed only

FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Effects = function(dataset,filenamelist,seed.list,seed.selected,burnin,x_var){
  nvarX = ncol(filenamelist[[1]]$alphadraw[,,1])
   QuadrantsCountsM = matrix(0,nrow=4, ncol=nvarX-1)
   QuadrantsCountsS = matrix(0,nrow=4, ncol=nvarX-1)
   DrawsAnalysis = c(seq(burnin+1,length(filenamelist[[seed.selected]]$LL_total),1))
   for(j in 2:nvarX)
   {  i = seed.selected
      #for(i in seed.list)
      {  ppM=pnM=npM=nnM=0
         ppS=pnS=npS=nnS=0
         for(r in DrawsAnalysis)
         { ppM = ppM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]>0,ifelse(filenamelist[[i]]$betaMdraw[r,2]>0,1,0),0)
           pnM = pnM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]>0,ifelse(filenamelist[[i]]$betaMdraw[r,2]<0,1,0),0)
           npM = npM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]<0,ifelse(filenamelist[[i]]$betaMdraw[r,2]>0,1,0),0)
           nnM = nnM + ifelse(filenamelist[[i]]$alphadraw[r,j,1]<0,ifelse(filenamelist[[i]]$betaMdraw[r,2]<0,1,0),0)
           ppS = ppS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]>0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
           pnS = pnS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]>0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
           npS = npS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]<0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
           nnS = nnS + ifelse(filenamelist[[i]]$alphadraw[r,j,2]<0,ifelse(filenamelist[[i]]$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
         }
         QuadrantsCountsM[,(j-1)]=c(ppM,pnM,npM,nnM)
         QuadrantsCountsS[,(j-1)]=c(ppS,pnS,npS,nnS)
       }
    }
  ProportionsM <- QuadrantsCountsM/length(DrawsAnalysis)
  ProportionsS <- QuadrantsCountsS/length(DrawsAnalysis)
  #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
  #par(oma=c(0,0,2,0));   par(mfcol=c(3*(nvarX-1)+1,length(seed.selected)))
  par(oma=c(0,0,2,0));   par(mfrow=c(nvarX-1,3))
  filename = filenamelist[[seed.selected]]
  ylimA = c( min(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2])-
                0.1*min(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2]),
             max(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2])+
             0.1*max(filename$alphadraw[-1:-burnin,-1,1],filename$alphadraw[-1:-burnin,-1,2])
          )
  ylimG = c( min(filename$gammabetaSdraw[-1:-burnin,-1])-0.1*min(filename$gammabetaSdraw[-1:-burnin,-1]),
             max(filename$gammabetaSdraw[-1:-burnin,-1])+0.1*max(filename$gammabetaSdraw[-1:-burnin,-1]))
  ylimB =  c( min(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)])-
                0.1*min(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)]),
              max(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)])+
                0.1*max(filename$betaMdraw[-1:-burnin,-1],filename$gammabetaSdraw[-1:-burnin,(nvarX+1)])
            )
  breaksCalc = (ylimG[2]-ylimG[1])/25
  for(p in 2:nvarX) {
    plot(filename$alphadraw[-1:-burnin,p,1], filename$betaMdraw[-1:-burnin,2], 
         main=bquote("Scatterplot of " ~ alpha[.(p-1)] ~ " and " ~ beta ~ " for variable " ~ .(x_var[p-1]) ~ ". Segment M"),
         #main=paste(x_var[p-1],"\n Segment M"),
       xlab=bquote(alpha[M][.(p-1)]),ylab=expression(beta[M]), xlim=ylimA,ylim=ylimB)
    abline(h=0,v=0,col="gray")
    if(ylimA[2]>0 & ylimB[2]>0){
      text(x=c(ylimA[2]),  y=c(ylimB[2]), 
         labels = c(paste( "I (",round(ProportionsM[1,p-1],4),")")), pos=2, font=2,cex=1)
    }
    if(ylimA[2]>0 & ylimB[1]<0){
      text(x=c(ylimA[2]),  y=c(ylimB[1]), 
         labels = c(paste( "II (",round(ProportionsM[2,p-1],4),")")), pos=2, font=2,cex=1)
    }
    if(ylimA[1]<0 & ylimB[1]<0){
      text(x=c(ylimA[1]),  y=c(ylimB[1]), 
         labels = c(paste("III (",round(ProportionsM[4,p-1],4),")")), pos=4, font=2,cex=1)
    }
    if(ylimA[1]<0 & ylimB[2]>0){
      text(x=c(ylimA[1]),  y=c(ylimB[2]), 
         labels = c(paste("IV (",round(ProportionsM[3,p-1],4),")")), pos=4, font=2,cex=1)
    }
    plot(filename$alphadraw[-1:-burnin,p,2], filename$gammabetaSdraw[-1:-burnin,nvarX+1],
         main=bquote("Scatterplot of " ~ alpha[.(p-1)] ~ " and " ~ beta ~ ". Variable " ~ .(x_var[p-1]) ~ ". Segment G"),
         xlab=bquote(alpha[S][.(p-1)]),ylab=expression(beta[S]), xlim=ylimA,ylim=ylimB)
    abline(h=0,v=0,col="gray")
    if(ylimA[2]>0 & ylimB[2]>0){
      text(x=c(ylimA[2]),  y=c(ylimB[2]), 
         labels = c(paste( "I (",round(ProportionsS[1,p-1],4),")")), pos=2, font=2,cex=1)
    }
    if(ylimA[2]>0 & ylimB[1]<0){
      text(x=c(ylimA[2]),  y=c(ylimB[1]), 
         labels = c(paste( "II (",round(ProportionsS[2,p-1],4),")")), pos=2, font=2,cex=1)
    }
    if(ylimA[1]<0 & ylimB[1]<0){
      text(x=c(ylimA[1]),  y=c(ylimB[1]), 
         labels = c(paste("III (",round(ProportionsS[4,p-1],4),")")), pos=4, font=2,cex=1)
    }
    if(ylimA[1]<0 & ylimB[2]>0){
      text(x=c(ylimA[1]),  y=c(ylimB[2]), 
         labels = c(paste("IV (",round(ProportionsS[3,p-1],4),")")), pos=4, font=2,cex=1)
    }
    hist(filename$gammabetaSdraw[-1:-burnin,p],
         main=bquote("Histogram of " ~ gamma[.(p-1)] ~ ". Variable " ~ .(x_var[p-1]) ~ ". Segment G"),
         #main=paste(x_var[p-1],"\n Segment G"),
         xlab=bquote(gamma[S][.(p-1)]), xlim=ylimG,
         breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,p]),(max(filename$gammabetaSdraw[-1:-burnin,p])+breaksCalc),breaksCalc )))
    abline(v=0,col="darkred",lwd=2)
  }
 
  #hist(filename$rhodraw[-1:-burnin], main=paste("Rho. Seed=",seed.list[seed.selected]), xlab=bquote(rho), xlim=c(0,1),
  #     breaks=50)
  #title(main=paste(dataset," Final Results. Binary mixture"),outer=T)
  #dev.off()
}

# testing
#debug(FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Effects)
#FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Effects(dataset="",filenamelist=output_list_b,seed.list,seed.selected=c(1),burnin=10)


FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Rho = function(dataset,filenamelist,seed.list,seed.selected,burnin){
  #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
  filename = filenamelist[[seed.selected]]
  hist(filename$rhodraw[-1:-burnin], main="", xlab=bquote(rho), xlim=c(0,1),col = "#75AADB",
       breaks=50)
}
  
FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_meanW = function(dataset,filenamelist,seed.list,seed.selected,burnin){
  #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
  filename = filenamelist[[seed.selected]]
  hist(rowMeans(filename$wdraw[,-1:-burnin]), main="", xlab="Mean w's", xlim=c(0,1),col = "#75AADB",
       breaks=50)
  #title(main=paste(dataset," Final Results. Binary mixture"),outer=T)
  #dev.off()
}


# # NOT DONE
# FUN_PDF_Mediation_FinalPlots_MSmixture_forShiny_PlotMCMC = function(dataset,filenamelist,seed.list,seed.selected,burnin,x_var){
#   nvarX = ncol(filenamelist[[1]]$alphadraw[,,1])
#   #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
#   par(oma=c(0,0,2,0));   par(mfcol=c(max(((nvarX-1)*3),8),length(seed.selected)*2))
#   filename = filenamelist[[seed.selected[1]]]
#   ylimAm=ylimAs = c( min(filename$alphadraw[-1:-burnin,,])-0.1*min(filename$alphadraw[-1:-burnin,,]),
#                      max(filename$alphadraw[-1:-burnin,,])+0.1*max(filename$alphadraw[-1:-burnin,,]))
#   ylimGBs = c( min(filename$gammabetaSdraw[-1:-burnin,])-0.1*min(filename$gammabetaSdraw[-1:-burnin,]),
#                max(filename$gammabetaSdraw[-1:-burnin,])+0.1*max(filename$gammabetaSdraw[-1:-burnin,]))
#   ylimBm =  c( min(filename$betaMdraw[-1:-burnin,])-0.1*min(filename$betaMdraw[-1:-burnin,]),
#                max(filename$betaMdraw[-1:-burnin,])+0.1*max(filename$betaMdraw[-1:-burnin,]))
#   ysigma =  c( min(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,]))-
#                   0.1*min(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,])),
#                max(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,]))+
#                   0.1*max(c(filename$sigma2mdraw[-1:-burnin,],filename$sigma2ydraw[-1:-burnin,])))
#    ylimLL =  c( min(filename$LL_total[-1:-burnin])-0.1*min(filename$LL_total[-1:-burnin]),
#                 max(filename$LL_total[-1:-burnin])+0.1*max(filename$LL_total[-1:-burnin]))
#   for( i in seed.selected)
#     { filename = filenamelist[[i]]
#       ylimA = c( min(filename$alphadraw[-1:-burnin,,])-0.1*min(filename$alphadraw[-1:-burnin,,]),
#                  max(filename$alphadraw[-1:-burnin,,])+0.1*max(filename$alphadraw[-1:-burnin,,]))
#       ylimG = c( min(filename$gammabetaSdraw[-1:-burnin,])-0.1*min(filename$gammabetaSdraw[-1:-burnin,]),
#                  max(filename$gammabetaSdraw[-1:-burnin,])+0.1*max(filename$gammabetaSdraw[-1:-burnin,]))
#       ylimB =  c( min(filename$betaMdraw[-1:-burnin,])-0.1*min(filename$betaMdraw[-1:-burnin,]),
#                   max(filename$betaMdraw[-1:-burnin,])+0.1*max(filename$betaMdraw[-1:-burnin,]))
#       matplot(filename$alphadraw[,,1],type='l',col=1:8,main=expression(alpha[M]),ylab=expression(alpha[M]),ylim=ylimAm);
#       matplot(filename$betaMdraw,type='l',col=1:3,main=expression(beta[M]),ylab=expression(beta[M]),ylim=ylimBm);
#       matplot(filename$alphadraw[,,2],type='l',col=1:8,main=expression(alpha[S]),ylab=expression(alpha[S]),ylim=ylimAs);
#       matplot(filename$gammabetaSdraw,type='l',col=1:8,main=expression(gamma[S],beta[S]),ylab=expression(gamma[S],beta[S]),ylim=ylimGBs);
#       matplot(filename$sigma2mdraw,type='l',col=1:4,main=expression(sigma^2[m]),ylab=expression(sigma^2[m]),ylim=ysigma);
#       matplot(filename$sigma2ydraw,type='l',col=1:4,main=expression(sigma^2[y]),ylab=expression(sigma^2[y]),ylim=ysigma);
#       matplot(filename$rhodraw,type='l',col=1:3,main=expression(rho),ylab=expression(rho),ylim=c(0.1,1));
#       plot(filename$LL_total,type='l',col=1,main="LL",ylab="LL",ylim=ylimLL)
#       done=8
#       if(nvarX>2) { while(done<((nvarX-1)*3))
#                     { plot.new()
#                       done=done+1
#                     }
#       }
#       for(p in 2:nvarX) {
#         plot(filename$alphadraw[-1:-burnin,p,1], filename$betaMdraw[-1:-burnin,2], main=paste("M segment. Seed=",seed.list[i]),
#            xlab=bquote(alpha[M][.(p-1)]),ylab=expression(beta[M]), xlim=ylimA,ylim=ylimB)
#         abline(h=0,v=0,col="gray")
#       }
#       for(p in 2:nvarX) {   plot(filename$alphadraw[-1:-burnin,p,2], filename$gammabetaSdraw[-1:-burnin,nvarX+1], main=paste("S segment. Seed=",seed.list[i]),
#            xlab=bquote(alpha[S][.(p-1)]),ylab=expression(beta[S]), xlim=ylimA,ylim=ylimB)
#         abline(h=0,v=0,col="gray")
#       }
#       breaksCalc = (ylimG[2]-ylimG[1])/25
#       for(p in 2:nvarX){
#          hist(filename$gammabetaSdraw[-1:-burnin,p], main=paste("S segment. Seed=",seed.list[i]), xlab=bquote(gamma[S][.(p-1)]), xlim=ylimG,
#            breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,p]),(max(filename$gammabetaSdraw[-1:-burnin,p])+breaksCalc),breaksCalc )))
#          abline(v=0,col="gray")
#       }
#       if(nvarX==2){ plot.new()  ;   plot.new() ;  plot.new()  ;   plot.new() ;   plot.new() }
#       if(nvarX==3){ plot.new()  ;   plot.new() ;   }
#     }
#   title(main=paste(dataset," Final Results. Binary mixture"),outer=T)
#   dev.off()
# }
# 
# 
# 


#---- (MS mixture model)  Parameters OUT MULTIX ---------------------------------------------------------------------------------
# this is based on the best seed
# outputs HDP intervals

FUN_PDF_Mediation_Parameters_MSmixture_forShiny  = function(filenamelist, seed.list, burnin)
  { filename = filenamelist[[seed.list]]
    nvarX = ncol(filename$alphadraw[-1:-burnin,,1])
    tempCIs_M = rbind(
       cbind(colMeans(filename$alphadraw[-1:-burnin,,1]), 
         t( apply(filename$alphadraw[-1:-burnin,,1],2,hdi,credMass = 0.95))),  # aM
       cbind(colMeans(filename$betaMdraw[-1:-burnin,]),
         t( apply(filename$betaMdraw[-1:-burnin,],2,hdi,credMass = 0.95))),     # bM
       cbind(colMeans(filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),
         t(apply((filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),2, hdi,credMass = 0.95)))  # abM
    )
    tempCIs_S = rbind(
       cbind(colMeans(filename$alphadraw[-1:-burnin,,2]), 
         t( apply(filename$alphadraw[-1:-burnin,,2],2,hdi,credMass = 0.95))),    # aS
       cbind(colMeans(filename$gammabetaSdraw[-1:-burnin,]),
         t( apply(filename$gammabetaSdraw[-1:-burnin,],2,hdi,credMass = 0.95))), # gbS
       cbind(colMeans(filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),
         t(apply((filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),2, hdi,credMass = 0.95))) # abS
    )
    tempCIs_Rho =  matrix(c(mean(filename$rhodraw[-1:-burnin]),
                     t(hdi(filename$rhodraw[-1:-burnin], credMass = 0.95))),ncol=3,nrow=1)
    colnames(tempCIs_M)   <- c("Mean","HPDI lower limit","HPDI upper limit")
    colnames(tempCIs_S)   <- c("Mean","HPDI lower limit","HPDI upper limit")
    colnames(tempCIs_Rho) <- c("Mean","HPDI lower limit","HPDI upper limit")
    
    rownames_list_M = c(rep(0,nrow(tempCIs_M)))
    rownames_list_M[nvarX+1] = "beta_{M_1}"
    rownames_list_M[nvarX+2] = "beta_{M_2}"
    for(i in 1:nvarX) {
      rownames_list_M[i]=paste0("alpha_{M_", i - 1, "}")
      rownames_list_M[nvarX+2+i]=paste0("alphabeta_{M_", i - 1, "}")
    }
    rownames_list_S = c(rep(0,nrow(tempCIs_S)))
    rownames_list_S[nvarX+nvarX+1] = "beta_{G}"
    for(i in 1:nvarX) {
      rownames_list_S[i]=paste0("alpha_{G_", i - 1, "}")
      rownames_list_S[nvarX+i]=paste0("gamma_{G_", i - 1, "}")
      rownames_list_S[nvarX+nvarX+1+i]=paste0("alphabeta_{G_", i - 1, "}")
    }
    rownames_list_Rho = expression(rho)
    rownames(tempCIs_M) <- rownames_list_M
    print(rownames_list_M)
    rownames(tempCIs_S) <- rownames_list_S
    rownames(tempCIs_Rho) <- rownames_list_Rho
   
    return(list(tempCIs_M,tempCIs_S,tempCIs_Rho))
}

# testing
#temp_test = FUN_PDF_Mediation_Parameters_MSmixture_forShiny (output_list_b,seed.list=c(1),burnin=10)

#  Old version
# FUN_PDF_Mediation_Parameters_MSmixture_forShiny  = function(filename,burnin)
# { nvarX = ncol(filename$alphadraw)
#   tempCIs = rbind(
#        t(apply(filename$alphadraw[-1:-burnin,,1],2,quantile,prob=c(0.025,0.05, 0.5,0.95,0.975))),                        # aM
#        t(apply(filename$betaMdraw[-1:-burnin,],2,quantile,prob=c(0.025,0.05, 0.5,0.95,0.975))),                          # bM
#        t(apply((filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),2, quantile,prob=c(0.025,0.05, 0.5,0.95,0.975))),  # abM
#        t(apply(filename$alphadraw[-1:-burnin,,2],2,quantile,prob=c(0.025,0.05, 0.5,0.95,0.975))),                        # aS
#        t(apply(filename$gammabetaSdraw[-1:-burnin,],2,quantile,prob=c(0.025,0.05, 0.5,0.95,0.975))),                     # gbS
#        t(apply((filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),2, quantile,prob=c(0.025,0.05, 0.5,0.95,0.975))),  # abS
#        t(quantile(filename$rhodraw[-1:-burnin], prob=c(0.025,0.05, 0.5,0.95,0.975)))
#        )
#   rownames(tempCIs) = c(rep("a",nvarX),rep("b",2),rep("ab",nvarX),rep("a",nvarX),rep("g",(nvarX)),"b",rep("ab",nvarX),"rho")
#   tempCIs_topaste = NULL
#   tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX*5+4),3],digits=3),"(",format(tempCIs[(nvarX*5+4),1],digits=3),",",format(tempCIs[(nvarX*5+4),5],digits=3),")") )
#   for(i in 2:nvarX){
#     tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[i,3],digits=3),"(",format(tempCIs[i,1],digits=3),",",format(tempCIs[i,5],digits=3),")") )      #a
#   }
#   tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX+2),3],digits=3),"(",format(tempCIs[(nvarX+2),1],digits=3),",",format(tempCIs[(nvarX+2),5],digits=3),")") )      #b
#   for(i in 2:nvarX){
#     tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX+2+i),3],digits=3),"(",format(tempCIs[(nvarX+2+i),1],digits=3),",",format(tempCIs[(nvarX+2+i),5],digits=3),")"))  #ab
#   }
#   for(i in 2:nvarX){
#     tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX*2+2+i),3],digits=3),"(",format(tempCIs[(nvarX*2+2+i),1],digits=3),",",format(tempCIs[(nvarX*2+2+i),5],digits=3),")") )      #a
#   }
#   tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX*4+3),3],digits=3),"(",format(tempCIs[(nvarX*4+3),1],digits=3),",",format(tempCIs[(nvarX*4+3),5],digits=3),")") )      #b
#   for(i in 2:nvarX){
#     tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX*4+3+i),3],digits=3),"(",format(tempCIs[(nvarX*4+3+i),1],digits=3),",",format(tempCIs[(nvarX*4+3+i),5],digits=3),")"))  #ab
#   }
#   for(i in 2:nvarX){
#     tempCIs_topaste = rbind(tempCIs_topaste,
#        paste(format(tempCIs[(nvarX*3+2+i),3],digits=3),"(",format(tempCIs[(nvarX*3+2+i),1],digits=3),",",format(tempCIs[(nvarX*3+2+i),5],digits=3),")") )  #g
#   }
#   return(list(tempCIs = tempCIs,tempCIs_topaste=tempCIs_topaste))
# }



#---- (MS mixture model)  Correlation of Ind Prob ---------------------------------------------------------------------------------
# TBD
FUN_PDF_Mediation_IndProbCorr_forShiny = function(dataset,filenamelist,seed.list,seed.index,pdfW,pdfH,burnin)
{  all_MeanRho  = matrix(0,ncol=1,nrow=length(seed.index))
   IndProb      = matrix(0,nrow=nrow(filenamelist[[1]]$wdraw),ncol=length(seed.index))
      for(i in 1:length(seed.index))
      {  IndProb[,i] = rowMeans(filenamelist[[seed.index[i]]]$wdraw[,-1:-burnin])
         all_MeanRho[i] = mean(filenamelist[[seed.index[i]]]$rhodraw[-1:-burnin])
      }
      pdf(paste(dataset,"_IndProbCorr", ".pdf", sep = ""), width=pdfW, height=pdfH)
      par(mfcol=c(length(seed.index),length(seed.index)), oma=c(0,0,2,0))
      for( i in 1:length(seed.index))
      { for( j in 1:length(seed.index))
        { plot(IndProb[,i],IndProb[,j],xlab=paste("seed=",seed.list[seed.index[i]]),ylab=paste("seed=",seed.list[seed.index[j]]),
               main=paste("r=",format(cor(IndProb[,i],IndProb[,j]),digit=4)),
          xlim=c(0,1),ylim=c(0,1))
        }
      }
   title(main=paste(dataset," Ind Probability Correlation"),outer=T)
   dev.off()
  return(list(IndProb=IndProb,all_MeanRho=all_MeanRho ))
}








