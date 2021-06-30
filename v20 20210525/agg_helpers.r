
# # # # testing
# testdf <- read.csv('sample_data_Loyalty.csv')
# x_vars <- 'x1'
# y_var <- 'y1'
# m_var <- 'avg_m'
# Aa_var <- 0.01
# Abg_var <- 0.01
# nu_var <- 5
# qm_var <- var(testdf[,m_var])
# qy_var <- var(testdf[,y_var])
# R_var <- 1000
# seed_var <- 123
# keep_var <- 10
# 
# input_list <- get_inputs_agg(testdf, x_vars, y_var, m_var, Aa_var, Abg_var, nu_var, qm_var, qy_var, R_var, seed_var, keep_var)
# output_list <- FUN_MediationGibbs_oneMediator_Aggregate_forShinyApp(input_list$Data, input_list$Prior, input_list$Mcmc)
# props  <- FUN_PDF_Mediation_AlphaBetaProportion_Aggregate_forShiny(output_list, burnin=100, x_vars)


#------(Aggregate model)--------------------------------#
# prep inputs for aggregate model: returns a list of inputs for Data, Mcmc, Prior

get_inputs_agg <-function(df, x_vars, y_var, m_var, 
                          Aa_var, Abg_var, nu_var, qm_var, qy_var,
                          R_var, seed_var, keep_var){
  
  df = as.matrix(df)
  
  X   = df[,x_vars]
  nhh = nrow(df)
  
  y = df[,y_var, drop=FALSE]
  m = df[,m_var, drop=FALSE]
  
  nvarX = 1 + length(x_vars) # number of X variables selected on the menu. intercept is added here
  nvarM = 1
  
  Data = list(
    X = cbind(rep(1, nhh), X), # datafile with only variables selected as X (independent variable(s))
    y = y, # datafile with only one variable selected as y (dependent variable)
    m = m,  # datafile with only one variable selected as m (mediator))
    Z = cbind(rep(1, nhh))
  ) 
  
  Prior = list(
    ma  = c(rep(0, nvarX)),
    Aa  = Aa_var * diag(nvarX), # "input Aa parameter, default is 0.01"
    mgb = c(rep(0, nvarX+nvarM)),
    Agb = Abg_var * diag(nvarX+nvarM), #"input Abg parameter, default is 0.01"
    nu  = nu_var,       # "input nu parameter, default is 5" 
    qy  = qy_var,  # "input qy parameter, default is var(y)"
    qm  = qm_var   # "input qy parameter, default is var(m)"
  )   
  
  R_draws = R_var
  Mcmc = list(
    R_draws = R_draws, # input value in thousands only, default is 10000
    keep    = keep_var, # input value, default is R/1000
    seed    = seed_var, #randomly generated number, but needs to be sent to the code
    flagBG  = 1, # these are for debugging, they are always 1
    flagA   = 1  # these are for debugging, they are always 1)
  ) 
  
  return(
    list(
      Data=Data,
      Prior=Prior,
      Mcmc=Mcmc
    )
  )
}



#------(Aggregate model)----------------------#
# function to run aggregate model: returns a list of results (8 variables)

FUN_MediationGibbs_oneMediator_Aggregate_forShinyApp = function (Data, Prior, Mcmc) 
{
  X = Data$X
  m = Data$m
  y = Data$y
  nobs = length(y)
  dima = ncol(X)   # dim of alpha
  dimg = ncol(X)   # dim of gamma
  dimbg = dimg + 1
  flagBG = Mcmc$flagBG
  flagA = Mcmc$flagA
  
  if (missing(Prior)) {
    ma = c(rep(0, dima))    	 # mean of alpha
    Aa = 0.01 * diag(dima)
    mbg = c(rep(0, (dimbg)))  	 # mean of betagamma
    Abg = 0.01 * diag(dimbg)
    nu = 5
    qy = var(y)
    qm = var(m)
  }
  else {
    ma = Prior$ma
    Aa = Prior$Aa
    mbg = Prior$mbg
    Abg = Prior$Abg
    nu = Prior$nu
    qy = Prior$qy
    qm = Prior$qm
  }
  R_draws = Mcmc$R_draws
  keep = Mcmc$keep
  
  alphadraw = matrix(double(floor(R_draws/keep) * dima), ncol = dima)
  betadraw  = rep(0, floor(R_draws/keep))
  gammadraw = matrix(double(floor(R_draws/keep) * dimg), ncol = dimg)
  sigma2mxdraw  = c(double(floor(R_draws/keep)))
  sigma2ymxdraw = c(double(floor(R_draws/keep)))    
  alpha         = c(rep(0.1, dima))
  LLdraw        = matrix(double(floor(R_draws/keep) * nobs), ncol = nobs)
  LLtotaldraw   = c(0,floor(R_draws/keep))
  
  sigma2ymx = 1
  sigma2mx  = 1
  
  #cat(" ", fill = TRUE);   cat("Starting Gibbs Sampler", fill = TRUE)    
  #itime = proc.time()[3]
  #cat("MCMC Iteration (est time to end -min) ", fill = TRUE)
  #fsh()
  
  set.seed(Mcmc$seed)
  for (rep in 1:R_draws) {
    
    # First conditional:  beta,gamma
    if(flagBG){
      RA = chol(Abg)
      W  = rbind(cbind(X, m), RA)
      z  = c(y, as.vector(RA %*% mbg))
      IR = backsolve(chol(crossprod(W)), diag(dimbg))
      btilde = crossprod(t(IR)) %*% crossprod(W, z)
      res = z - W %*% btilde
      sbg = t(res) %*% res
      sigma2ymx = (nu * qy + sbg)/rchisq(1, nu + nobs)
      gammabeta = btilde + as.vector(sqrt(sigma2ymx)) * IR %*% rnorm(dimbg)
    }
    # Second conditional:  alpha
    if(flagA){
      RA = chol(Aa)
      W  = rbind(X, RA)
      z  = c(m, as.vector(RA %*% ma))
      IR = backsolve(chol(crossprod(W)), diag(dima))
      atilde = crossprod(t(IR)) %*% crossprod(W, z)
      res = z - W %*% atilde
      sa  = t(res) %*% res
      sigma2mx = (nu * qm + sa)/rchisq(1, nu + nobs)
      alpha    = atilde + as.vector(sqrt(sigma2mx)) * IR %*% rnorm(dima)
      
    }
    if (rep%%keep == 0) {
      mkeep                = rep/keep
      alphadraw[mkeep, ]   = alpha
      betadraw[mkeep]      = gammabeta[dimg+1]
      gammadraw[mkeep, ]   = gammabeta[1:dimg]
      sigma2mxdraw[mkeep]  = sigma2mx
      sigma2ymxdraw[mkeep] = sigma2ymx
      LLdraw[mkeep,]       = -0.5*(log(c(2*pi*sigma2ymx))+log(c(2*pi*sigma2mx))) - 0.5*(m-X%*%(alpha))^2/ c(sigma2mx) -
        0.5*(y-cbind(X,m)%*%gammabeta)^2 / c(sigma2ymx) 
      LLtotaldraw[mkeep]   = sum(LLdraw[mkeep,])
      # if (rep%%500 == 0) {
      #     ctime = proc.time()[3]
      #     timetoend = ((ctime - itime)/rep) * (R - rep)
      #     cat(" ", rep, " (", round(timetoend/60, 1), ")"," ll=",LLtotaldraw[mkeep],
      #         "b=",gammabeta[dimg+1],"a=",alpha,"g=",gammabeta[1:dimg],"symx=",sigma2ymx, "smx=",sigma2mx, fill = TRUE)
      #     #fsh()
      #    }
    }
  }
  #ctime = proc.time()[3]
  #cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2),"\n")
  return(list(alphadraw = alphadraw, 
              betadraw  = betadraw, 
              gammadraw = gammadraw, 
              sigma2mxdraw  = sigma2mxdraw,
              sigma2ymxdraw = sigma2ymxdraw,
              LL      = LLdraw, 
              LLtotal = LLtotaldraw,
              Mcmc    = Mcmc))
}




#---- (Aggregate model Results Page) --------------------------#
# Proportion of joint (alpha, beta) in each quadrant MULTIX

FUN_PDF_Mediation_AlphaBetaProportion_Aggregate_forShiny  = function(filename, burnin, x_vars)
  { 
  nvarX = ncol(filename$alphadraw) - 1
  QuadrantsCounts = array(0, dim=c(4, nvarX))
  DrawsAnalysis   = c(seq(from = burnin+1, to = length(filename$LL_total), by = 1))
  pp=pn=np=nn=0
  for(j in 2:(nvarX+1)){  
    for(r in DrawsAnalysis){ 
      pp = pp + ifelse(filename$alphadraw[r,j,1]>0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
      pn = pn + ifelse(filename$alphadraw[r,j,1]>0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
      np = np + ifelse(filename$alphadraw[r,j,1]<0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
      nn = nn + ifelse(filename$alphadraw[r,j,1]<0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
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

#---- (Aggregate model Results Page) --------------------------#
# Proportion of joint (alpha, beta) in each quadrant MULTIX

FUN_PDF_Mediation_LMD_NR_Aggregate_forShiny  = function(filename, burnin)
  { outputTable = matrix(round(logMargDenNR(filename$LL_total[-1:-burnin]),2),1,1)
    rownames(outputTable) <- c("LMD NR")
    return(outputTable)
}

#---- (Aggregate model Results Page) --------------------------#
# HDP MULTIX

FUN_PDF_Mediation_HDPI_Aggregate_forShiny  = function(filename, burnin, x_vars)
  { nvarX = ncol(filename$alphadraw)
  outputTable = NULL
  for(j in 1:(nvarX)){  
    tempa = hdi(filename$alphadraw[-1:-burnin,j,1], credMass = 0.95)
    outputTable = rbind(outputTable,c(mean(filename$alphadraw[-1:-burnin,j,1]),tempa))
  }
  tempb = hdi(filename$betaMdraw[-1:-burnin,1], credMass = 0.95)
  outputTable = rbind(outputTable,c(mean(filename$betaMdraw[-1:-burnin,1]),tempb))
  for(j in 1:(nvarX)){  
    tempg = hdi(filename$gammabetaSdraw[-1:-burnin,j], credMass = 0.95)
    outputTable = rbind(outputTable,c(mean(filename$gammabetaSdraw[-1:-burnin,j]),tempg))
  }
  
  #data.frame(outputTable)
  
  rownames_list = c(rep(0,nrow(outputTable)))
  rownames_list[1] = "alpha_{0}"
  rownames_list[nvarX+1] = "beta"
  rownames_list[nvarX+2] = "gamma_{0}"
  for(i in 2:(nvarX)) {
    rownames_list[i]=paste0("alpha_{", i - 1, "}")
    rownames_list[nvarX+1+i]=paste0("gamma_{", i - 1, "}")
  }
  #outputTable = cbind(rownames_list,outputTable)
  
  colnames(outputTable) <- c("Mean","Lower limit","Upper limit")
  
  # rownames(outputTable) <- rownames_list
  # colnames(outputTable) <- c("Mean","Lower limit","Upper limit")
  # rownames_list = c(rep(0,nrow(outputTable)))
  # rownames_list[1] = bquote(alpha[0]) 
  # rownames_list[nvarX+1] = paste(expression(beta))
  # rownames_list[nvarX+2] = paste(expression(gamma[0]))
  # for(i in 2:(nvarX)) {
  #   rownames_list[i]=paste(expression(alpha[i-1]))
  #   rownames_list[nvarX+1+i]=paste(expression(gamma[i-1]))
  # }
  rownames(outputTable) <- rownames_list
  # return(list(Proportions = QuadrantsCounts/length(DrawsAnalysis)))
  #return(grid.table(outputTable))
  return(outputTable)
  #return(htmlTable(outputTable))
}


#-----(Aggregate model)-----------------------#
# function to generate PDF scatterplots of postreior draws of alphas and beta for aggregate model

FUN_PDF_Mediation_ScatterPlots_Aggregate_forShiny_Plot = function(dataset, filename, burnin, x_vars) {
    # dataset = the name (string) of the data file that needs to be displayed (i.e. "Loyalty data")
    # filename = output file of MCMC (list of objects/tables)
    # burnin = number of saved draws to exclude from plotting
    # x_vars = vector of names of variables in X, excluding intercept
  nvarX = ncol(filename$alphadraw) - 1
  #pdf(paste(dataset, "Posterior Draws.pdf", sep = " "), width=pdfW, height=pdfH)
  par(oma=c(0,0,2,0));
  QuadrantsCounts = array(0, dim=c(4, nvarX))
  print(burnin)
  print(length(filename$LL_total))
  DrawsAnalysis   = c(seq(from = burnin+1, to = length(filename$LL_total), by = 1))
  pp=pn=np=nn=0
  for(j in 1:nvarX){  
    for(r in DrawsAnalysis){ 
      pp = pp + ifelse(filename$alphadraw[r,j+1,1]>0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
      pn = pn + ifelse(filename$alphadraw[r,j+1,1]>0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
      np = np + ifelse(filename$alphadraw[r,j+1,1]<0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
      nn = nn + ifelse(filename$alphadraw[r,j+1,1]<0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
    }
    QuadrantsCounts[,j]=c(pp,pn,np,nn)
    pp=pn=np=nn=0
  }
  Proportions <- QuadrantsCounts/length(DrawsAnalysis)
  
  par(mfrow=c(nvarX,2))
  ylimB = c( min(filename$betaMdraw[-1:-burnin])-0.1*min(filename$betaMdraw[-1:-burnin,1]),
           max(filename$betaMdraw[-1:-burnin])+0.1*max(filename$betaMdraw[-1:-burnin,1]))
  for(j in 1:nvarX){
    ylimA = c( min(filename$alphadraw[-1:-burnin,j+1,1])-0.1*min(filename$alphadraw[-1:-burnin,j+1,1]),
             max(filename$alphadraw[-1:-burnin,j+1,1])+0.1*max(filename$alphadraw[-1:-burnin,j+1,1]))
    ylimG = c( min(filename$gammabetaSdraw[-1:-burnin,j+1])-0.1*min(filename$gammabetaSdraw[-1:-burnin,j+1]),
             max(filename$gammabetaSdraw[-1:-burnin,j+1])+0.1*max(filename$gammabetaSdraw[-1:-burnin,j+1]))
    breaksCalc = (ylimG[2]-ylimG[1])/25
    plot(filename$alphadraw[-1:-burnin,j+1,1], filename$betaMdraw[-1:-burnin,1], 
         main=bquote("Scatterplot of " ~ alpha[.(j)] ~ " and " ~ beta ~ " for variable " ~ .(x_vars[j])),
         xlab=bquote(alpha[.(j)]), ylab=expression(beta), xlim=ylimA,ylim=ylimB)
    abline(h=0,v=0,col="gray")
    if(ylimA[2]>0 & ylimB[2]>0){
      text(x=c(ylimA[2]),  y=c(ylimB[2]), 
         labels = c(paste( "I (",round(Proportions[1,j],4),")")), pos=2, font=2,cex=0.9)
    }
    if(ylimA[2]>0 & ylimB[1]<0){
      text(x=c(ylimA[2]),  y=c(ylimB[1]), 
         labels = c(paste( "II (",round(Proportions[2,j],4),")")), pos=2, font=2,cex=0.9)
    }
    if(ylimA[1]<0 & ylimB[1]<0){
      text(x=c(ylimA[1]),  y=c(ylimB[1]), 
         labels = c(paste("IV (",round(Proportions[4,j],4),")")), pos=4, font=2,cex=0.9)
    }
    if(ylimA[1]<0 & ylimB[2]>0){
      text(x=c(ylimA[1]),  y=c(ylimB[2]), 
         labels = c(paste("IV (",round(Proportions[3,j],4),")")), pos=4, font=2,cex=0.9)
    }
    hist(filename$gammabetaSdraw[-1:-burnin,j+1], main=bquote(paste("Histogram of ",gamma[.(j)])),
         xlab=bquote(gamma[.(j)]), xlim=ylimG,
         breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,j+1]),(max(filename$gammabetaSdraw[-1:-burnin,j+1])+breaksCalc),breaksCalc )))
    abline(v=0,col="darkred", lwd=2)
  }
  #title(main=paste(dataset, "Aggregate Model.  Posterior Draws of parameters."),outer=T)
# dev.off()
}
