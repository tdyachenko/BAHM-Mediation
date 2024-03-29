


#---- (Aggregate model Results Page) --------------------------#
# report fit LMD NR for the aggregate model. Why do we need a table here?


# --- DIC (both models defined by ModelFlag (2022-02-01) ----------------------------------------------

FUN_DIC_mediation = function(Data, McmcOutput,burnin,ModelFlag)
{
  #load("inputdata.RData")
  Data <- list(y = Data$y, 
               X = as.matrix(Data$X), 
               m = Data$m)
  
  y = Data$y
  X = Data$X
  m = Data$m
  R = length(McmcOutput[[1]]$LL_total[-1:-burnin]) 
  DIC= c(rep(0,length(McmcOutput)))
  
  # DIC = Dbar + pD = Dhat + 2 pD (page 603, eq. 36-37, Spiegelhalter et al, 2002, JRSS)
  # thus, DIC = 2*Dbar - Dhat
  # pD is 'the effective number of parameters', where pD = Dbar - Dhat, or 
  # pD = 0.5*bar(var(D(theta)))  from Gelman et al (2014, eq. 10, p.1002)
  
  # Dhat is a point estimate of the deviance obtained by substituting in the posterior means thetabar
  # Dhat = - 2 * log(p( y | thetabar ))
  # Dbar is the posterior mean of the deviance
  
  # Calculate thetabar: the posterior means of all parameters.
  for(i in 1:length(McmcOutput))
  { m_alphadraw_M = colMeans(McmcOutput[[i]]$alphadraw[-1:-burnin,,1])
  m_alphadraw_S = colMeans(McmcOutput[[i]]$alphadraw[-1:-burnin,,2])
  m_betaMdraw = colMeans(McmcOutput[[i]]$betaMdraw[-1:-burnin,])
  m_gammabetaSdraw = colMeans(McmcOutput[[i]]$gammabetaSdraw[-1:-burnin,])
  #m_lambdadraw = colMeans(McmcOutput[[i]]$lambdadraw[-1:-burnin,])
  m_sigma2mdraw = colMeans(McmcOutput[[i]]$sigma2mdraw[-1:-burnin,])
  m_sigma2ydraw = colMeans(McmcOutput[[i]]$sigma2ydraw[-1:-burnin,])
  m_wdraw = rowMeans(McmcOutput[[i]]$wdraw[,-1:-burnin])
  
  # Calculate D(thetabar)
  LikD = exp(-0.5*(log(c(2*pi*m_sigma2ydraw[2]))+log(c(2*pi*m_sigma2mdraw[2]))) - 
               0.5*(m-X%*%(m_alphadraw_S))^2/ c(m_sigma2mdraw[2]) - 
               0.5*(y-cbind(X,m)%*%m_gammabetaSdraw)^2 / c(m_sigma2ydraw[2]))
  if(ModelFlag==2)
  { LikM = exp(-0.5*(log(c(2*pi*m_sigma2ydraw[1]))+log(c(2*pi*m_sigma2mdraw[1]))) -
                 0.5*(m-X%*%(m_alphadraw_M))^2/ c(m_sigma2mdraw[1]) -
                 0.5*(y-cbind(X[,1],m)%*%m_betaMdraw)^2 / c(m_sigma2ydraw[1])) 
  indexM = which(m_wdraw>0.5)
  Dhat = -2*(sum(log(LikD[-indexM]))+sum(log(LikM[indexM])))
  }
  else{
    Dhat = -2*sum(log(LikD))
  }
  #pD = 0.5*mean(var(McmcOutput[[i]]$LL_total[-1:-burnin]))  # See Gelman et al 2004
  Dbar = mean(-2*McmcOutput[[i]]$LL_total[-1:-burnin])
  DIC[i] = 2*Dbar - Dhat
  #DIC[i] = Dhat + 2 * pD
  }
  return(DIC)
}


# -------------------------------------------------------------------------------
# Proportion of joint (alpha, beta) in each quadrant MULTIX
# -------------------------------------------------------------------------------
# updated 7-13-2021
# -------------------------------------------------------------------------------

FUN_PDF_Mediation_AlphaBetaProportion_forShiny  = function(model,filename, burnin, x_vars)
{ 
  if(model==1)
  {
    nvarX = ncol(filename$alphadraw[,,2])
    QuadrantsCounts = array(0, dim=c(4, (nvarX-1)))
    DrawsAnalysis   = c(seq(from = burnin+1, to = length(filename$LL_total), by = 1))
    pp=pn=np=nn=0
    for(j in 2:(nvarX)){  
      for(r in DrawsAnalysis){ 
        pp = pp + ifelse(filename$alphadraw[r,j,2]>0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
        pn = pn + ifelse(filename$alphadraw[r,j,2]>0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
        np = np + ifelse(filename$alphadraw[r,j,2]<0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
        nn = nn + ifelse(filename$alphadraw[r,j,2]<0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
      }
      QuadrantsCounts[,(j-1)]=c(pp,pn,np,nn)
      pp=pn=np=nn=0
    }
    Proportions <- QuadrantsCounts/length(DrawsAnalysis)
    colnames(Proportions) <- x_vars
    rownames(Proportions) <- c("I (++)","II (+-)","III (-+)","IV (--)")
    # return(list(Proportions = QuadrantsCounts/length(DrawsAnalysis)))
    return(Proportions)
  }
}


# -------------------------------------------------------------------------------
# HDPI MULTIX       updated 7-13-2021 for aggregate only
# -------------------------------------------------------------------------------
# 
FUN_PDF_Mediation_HDPI_forShiny  = function(model,filename, burnin, x_vars,  m_var, CIband)
{ 
  # Aggregate ONLY
  if(model==1) 
  { nvarX = ncol(filename$alphadraw)
    tempa = t(apply(filename$alphadraw[-1:-burnin,,2],2,hdi,credMass = CIband))
    outputTable = cbind(c(colMeans(filename$alphadraw[-1:-burnin,,2])),tempa)
    tempgb = t(apply(filename$gammabetaSdraw[-1:-burnin,],2,hdi,credMass = CIband))
    outputTable = rbind(outputTable,
                        cbind(c(colMeans(filename$gammabetaSdraw[-1:-burnin,])),tempgb))
    
    rownames_list = c(rep(0,nrow(outputTable)))
    rownames_list[1] = "alpha_{0}"
    rownames_list[nvarX+1] = "gamma_{0}"
    rownames_list[2*nvarX+1] = "beta"
    
    for(i in 2:(nvarX)) {
      rownames_list[i]=paste0("alpha_{", i - 1, "}")
      rownames_list[nvarX+i]=paste0("gamma_{", i - 1, "}")
    }
    #outputTable = cbind(rownames_list,outputTable)
    
    colnames(outputTable) <- c("Mean","HPDI Lower limit","HPDI Upper limit")
    rownames(outputTable) <- rownames_list
    
    outputTable <- cbind("Variable" = c("Intercept", x_vars,"Intercept", x_vars, m_var), outputTable)
    
    # return(list(Proportions = QuadrantsCounts/length(DrawsAnalysis)))
    #return(grid.table(outputTable))
    return(outputTable)
    #return(htmlTable(outputTable))
  }
}

# -------------------------------------------------------------------------------
# HDPI MULTIX       updated 7-22-2020 for BM only
# -------------------------------------------------------------------------------
# this is based on the best seed, # outputs HPD intervals

FUN_PDF_Mediation_Parameters_MSmixture_forShiny  = function(filenamelist, x_vars, m_var, z_vars, seed.selected, burnin, CIband)
  { filename = filenamelist[[seed.selected]]
    nvarX = ncol(filename$alphadraw[-1:-burnin,,1])
    tempCIs_M = rbind(
       cbind(colMeans(filename$alphadraw[-1:-burnin,,1]), 
         t( apply(filename$alphadraw[-1:-burnin,,1],2,hdi,credMass = CIband))),  # aM
       cbind(colMeans(filename$betaMdraw[-1:-burnin,]),
         t( apply(filename$betaMdraw[-1:-burnin,],2,hdi,credMass = CIband)))     # bM
       #cbind(colMeans(filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),
         #t(apply((filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),2, hdi,credMass = 0.95)))  # abM
    )
    tempCIs_S = rbind(
       cbind(colMeans(filename$alphadraw[-1:-burnin,,2]), 
         t( apply(filename$alphadraw[-1:-burnin,,2],2,hdi,credMass = CIband))),    # aS
       cbind(colMeans(filename$gammabetaSdraw[-1:-burnin,]),
         t( apply(filename$gammabetaSdraw[-1:-burnin,],2,hdi,credMass = CIband))) # gbS
       #cbind(colMeans(filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),
         #t(apply((filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),2, hdi,credMass = 0.95))) # abS
    )
    tempCIs_Rho =  matrix(c(mean(filename$rhodraw[-1:-burnin]),
                     t(hdi(filename$rhodraw[-1:-burnin], credMass = CIband))),ncol=3,nrow=1)
    tempCIs_lambda = cbind(colMeans(filename$lambdadraw[-1:-burnin, , drop = FALSE]),
         t(apply(filename$lambdadraw[-1:-burnin, , drop = FALSE],2,hdi,credMass = .95)))
    #templambda = t(apply(filename$lambdadraw[-1:-burnin, , drop = FALSE],2,hdi,credMass = CIband))
    #tempCIs_lambda = cbind(c(colMeans(filename$lambdadraw[-1:-burnin, , drop = FALSE])),templambda)
    
    colnames(tempCIs_M)   <- c("Mean","HPDI lower limit","HPDI upper limit")
    colnames(tempCIs_S)   <- c("Mean","HPDI lower limit","HPDI upper limit")
    colnames(tempCIs_Rho) <- c("Mean","HPDI lower limit","HPDI upper limit")
    colnames(tempCIs_lambda) <- c("Mean","HPDI lower limit","HPDI upper limit")
    
    rownames_list_M = c(rep(0,nrow(tempCIs_M)))
    rownames_list_M[nvarX+1] = "beta_{M_0}"
    rownames_list_M[nvarX+2] = "beta_{M_1}"
    for(i in 1:nvarX) {
      rownames_list_M[i]=paste0("alpha_{M_", i - 1, "}")
      #rownames_list_M[nvarX+2+i]=paste0("alphabeta_{M_", i - 1, "}")
    }
    rownames_list_S = c(rep(0,nrow(tempCIs_S)))
    rownames_list_S[nvarX+nvarX+1] = "beta_{S}"
    for(i in 1:nvarX) {
      rownames_list_S[i]=paste0("alpha_{S_", i - 1, "}")
      rownames_list_S[nvarX+i]=paste0("gamma_{S_", i - 1, "}")
      #rownames_list_S[nvarX+nvarX+1+i]=paste0("alphabeta_{S_", i - 1, "}")
    }
    rownames_list_Rho = expression(rho)
    rownames(tempCIs_M) <- rownames_list_M
    # print(rownames_list_M)
    rownames(tempCIs_S) <- rownames_list_S
    rownames(tempCIs_Rho) <- rownames_list_Rho
    
    nvarZ = nrow(tempCIs_lambda)
    rownames_list_Lambda = c(rep(0,nvarZ))
    for(i in 1:nvarZ) {
      rownames_list_Lambda[i]=paste0("lambda_{", i - 1, "}")
    }
    rownames(tempCIs_lambda) <- rownames_list_Lambda
    
    tempCIs_M <- cbind("Variable" = c("Intercept", x_vars,"Intercept", m_var), tempCIs_M)
    tempCIs_S <- cbind("Variable" = c("Intercept", x_vars,"Intercept", x_vars, m_var), tempCIs_S)
    
    print("=====")
    print(z_vars)
    print(tempCIs_lambda)
    print(filename$lambdadraw[100:110, ])
    print("=====")
    
    tempCIs_lambda <- cbind("Variable" = c("Intercept", z_vars), tempCIs_lambda)
    
    return(list(tempCIs_M,tempCIs_S,tempCIs_Rho,tempCIs_lambda))
}

# testing
#temp_test = FUN_PDF_Mediation_Parameters_MSmixture_forShiny (output_list_b,seed.list=c(1),burnin=10)


#-----(Aggregate model)-----------------------#
# function to generate PDF scatterplots of postreior draws of alphas and beta for aggregate model

FUN_PDF_Mediation_ScatterPlots_forShiny = function(model,dataset, filename, burnin, x_vars) {
    # dataset = the name (string) of the data file that needs to be displayed (i.e. "Loyalty data")
    # filename = output file of MCMC (list of objects/tables)
    # burnin = number of saved draws to exclude from plotting
    # x_vars = vector of names of variables in X, excluding intercept
if(model==1){  
  nvarX = ncol(filename$alphadraw[,,2])
  #pdf(paste(dataset, "Posterior Draws.pdf", sep = " "), width=pdfW, height=pdfH)
  par(oma=c(0,0,2,0));
  QuadrantsCounts = array(0, dim=c(4, (nvarX-1)))
  #print(burnin)
  #print(length(filename$LL_total))
  DrawsAnalysis   = c(seq(from = burnin+1, to = length(filename$LL_total), by = 1))
  pp=pn=np=nn=0
  for(j in 2:(nvarX)){  
    for(r in DrawsAnalysis){ 
      # this is the same calculation as in in FUN_PDF_Mediation_AlphaBetaProportion_forShiny above
      # NOTE: can we remove the duplication here? 7-13-2021 TD
      pp = pp + ifelse(filename$alphadraw[r,j,2]>0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
      pn = pn + ifelse(filename$alphadraw[r,j,2]>0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
      np = np + ifelse(filename$alphadraw[r,j,2]<0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]>0,1,0),0)
      nn = nn + ifelse(filename$alphadraw[r,j,2]<0,ifelse(filename$gammabetaSdraw[r,(nvarX+1)]<0,1,0),0)
    }
    QuadrantsCounts[,(j-1)]=c(pp,pn,np,nn)
    pp=pn=np=nn=0
  }
  Proportions <- QuadrantsCounts/length(DrawsAnalysis)
  
  par(mfrow=c((nvarX-1),2))
  ylimB = c( min(filename$gammabetaSdraw[-1:-burnin,nvarX+1])-0.1*min(filename$gammabetaSdraw[-1:-burnin,nvarX+1]),
             max(filename$gammabetaSdraw[-1:-burnin,nvarX+1])+0.1*max(filename$gammabetaSdraw[-1:-burnin,nvarX+1]))
  for(j in 2:(nvarX)){
    ylimA = c( min(filename$alphadraw[-1:-burnin,j,2])-0.1*min(filename$alphadraw[-1:-burnin,j,2]),
               max(filename$alphadraw[-1:-burnin,j,2])+0.1*max(filename$alphadraw[-1:-burnin,j,2]))
    ylimG = c( min(filename$gammabetaSdraw[-1:-burnin,j])-0.1*min(filename$gammabetaSdraw[-1:-burnin,j]),
               max(filename$gammabetaSdraw[-1:-burnin,j])+0.1*max(filename$gammabetaSdraw[-1:-burnin,j]))
    breaksCalc = (ylimG[2]-ylimG[1])/25
    plot(filename$alphadraw[-1:-burnin,j,2], filename$gammabetaSdraw[-1:-burnin,nvarX+1], 
         main=bquote("Scatterplot of " ~ alpha[.(j-1)] ~ " and " ~ beta ~ " for variable " ~ .(x_vars[j-1])),
         xlab=bquote(alpha[.(j-1)]), ylab=expression(beta), xlim=ylimA,ylim=ylimB)
    abline(h=0,v=0,col="gray")
    if(ylimA[2]>0 & ylimB[2]>0){
      text(x=c(ylimA[2]),  y=c(ylimB[2]), 
         labels = c(paste( "I (",round(Proportions[1,j-1],4),")")), pos=2, font=2,cex=0.9)
    }
    if(ylimA[2]>0 & ylimB[1]<0){
      text(x=c(ylimA[2]),  y=c(ylimB[1]), 
         labels = c(paste( "II (",round(Proportions[2,j-1],4),")")), pos=2, font=2,cex=0.9)
    }
    if(ylimA[1]<0 & ylimB[1]<0){
      text(x=c(ylimA[1]),  y=c(ylimB[1]), 
         labels = c(paste("IV (",round(Proportions[4,j-1],4),")")), pos=4, font=2,cex=0.9)
    }
    if(ylimA[1]<0 & ylimB[2]>0){
      text(x=c(ylimA[1]),  y=c(ylimB[2]), 
         labels = c(paste("IV (",round(Proportions[3,j-1],4),")")), pos=4, font=2,cex=0.9)
    }
    hist(filename$gammabetaSdraw[-1:-burnin,j], main=bquote(paste("Histogram of ",gamma[.(j-1)])),
         xlab=bquote(gamma[.(j-1)]), xlim=ylimG,
         breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,j]),(max(filename$gammabetaSdraw[-1:-burnin,j])+breaksCalc),breaksCalc )))
    abline(v=0,col="darkred", lwd=2)
  }
  #title(main=paste(dataset, "Aggregate Model.  Posterior Draws of parameters."),outer=T)
  # dev.off()
  }
}




# -------------------------------------------------------------------------------
# MCMC trace plots 
# -------------------------------------------------------------------------------
# new function for both aggregate and BM models (7-12-2021)
# there is no output of MCMCs for the aggregate model. We don't need it now, but the code is here for future use
# -------------------------------------------------------------------------------
FUN_PDF_MCMC_Mediation_forShiny = function(model,dataset,filenamelist,seed.index,seed.list,burnin)
{  
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
   ylimLL =  c( min(filename$LL_total[-1:-2])-0.01*min(filename$LL_total[-1:-2]),  # used burnin before
                max(filename$LL_total[-1:-2])+0.01*max(filename$LL_total[-1:-2]))
   
  if(model==1)
  {
   #pdf(paste(dataset,"_MCMC_Aggregate", ".pdf", sep = ""), width=pdfW, height=pdfH)
   par(oma=c(0,0,2,0));   par(mfcol=c(5,length(seed.index)),mai = c(0.5, 0.3, 0.4, 0.3))
   for( i in seed.index)
   { filename = filenamelist[[i]]
     matplot(filename$alphadraw[,,2],type='l',col=1:8,main=expression(alpha[S]),ylab=expression(alpha[S]),ylim=ylimAs);
     #matplot(filename$gammabetaMdraw,type='l',col=1:4,main=expression(gamma[M],beta[M]),ylab=expression(gamma[M],beta[M]),ylim=ylimGBs);
     matplot(filename$gammabetaSdraw,type='l',col=1:8,main=expression(gamma[S],beta[S]),ylab=expression(gamma[S],beta[S]),ylim=ylimGBs);
     matplot(filename$sigma2mdraw,type='l',col=1:4,main=expression(sigma[m]^2),ylab=expression(sigma[m]^2),ylim=ysigma);
     matplot(filename$sigma2ydraw,type='l',col=1:4,main=expression(sigma[y]^2),ylab=expression(sigma[y]^2),ylim=ysigma);
     plot(filename$LL_total,type='l',col=1,main="LL",ylab="LL",ylim=ylimLL)
   }
   title(main=paste(dataset),outer=T)
   #pdf(paste("plot", j, ".pdf", sep = ""),width=20, height=10)
   dev.off()
  }
  if(model==2)
  {
   ylimL =  c( min(filename$lambdadraw[-1:-burnin,])-0.1*min(filename$lambdadraw[-1:-burnin,]),  # used burnin before
                max(filename$lambdadraw[-1:-burnin,])+0.1*max(filename$lambdadraw[-1:-burnin,]))
   #pdf(paste(dataset,"_MCMC_MSmixture", ".pdf", sep = ""), width=pdfW, height=pdfH)
   par(oma=c(0,0,2,0));   par(mfcol=c(8,length(seed.index)),mai = c(0.5, 0.3, 0.4, 0.3))
   for( i in seed.index)
   { filename = filenamelist[[i]]
     matplot(filename$alphadraw[,,1],type='l',col=1:8,main=expression(alpha[M]),ylab=expression(alpha[M]),ylim=ylimAm);
     matplot(filename$alphadraw[,,2],type='l',col=1:8,main=expression(alpha[S]),ylab=expression(alpha[S]),ylim=ylimAs);
     matplot(filename$betaMdraw,type='l',col=1:3,main=expression(beta[M]),ylab=expression(beta[M]),ylim=ylimBm);
     matplot(filename$gammabetaSdraw,type='l',col=1:8,main=expression(gamma[S],beta[S]),ylab=expression(gamma[S],beta[S]),ylim=ylimGBs);
     matplot(filename$sigma2mdraw,type='l',col=1:4,main=expression(sigma[m]^2),ylab=expression(sigma[m]^2),ylim=ysigma);
     matplot(filename$sigma2ydraw,type='l',col=1:4,main=expression(sigma[y]^2),ylab=expression(sigma[y]^2),ylim=ysigma);
     matplot(filename$lambdadraw,type='l',col=1:3,main=expression(lambda),ylab=expression(lambda),ylim=ylimL);
     plot(filename$LL_total,type='l',col=1,main="LL",ylab="LL",ylim=ylimLL)
   }
   #title(main=paste(dataset),outer=T)
   #pdf(paste("plot", j, ".pdf", sep = ""),width=20, height=10)
   #dev.off()
  }
}







