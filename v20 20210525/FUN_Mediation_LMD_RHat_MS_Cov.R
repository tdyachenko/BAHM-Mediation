
#  NOT DEVELOPED FOR MORE THAN 2 SOLUTIONS !!!

library(bayesm);   library(coda); library(gtools)
FUN_Mediation_LMD_RHat_MS_cov = function(datafile,seed.index,burnin,RhatCutoff)
{ 
  # This is the finction without parallel computing
  
  mcmc.list =   function (...)  # function copied from coda package, otherwise doparallel does not see it
  {
    x <- list(...)
    if (length(x) == 1 && is.list(x[[1]])) 
      x <- x[[1]]
    if (!all(unlist(lapply(x, is.mcmc)))) 
      stop("Arguments must be mcmc objects")
    nargs <- length(x)
    if (nargs >= 2) {
      xmcpar <- lapply(x, mcpar)
      if (!all(unlist(lapply(xmcpar, "==", xmcpar[[1]])))) 
        stop("Different start, end or thin values in each chain")
      xnvar <- lapply(x, nvar)
      if (!all(unlist(lapply(xnvar, "==", xnvar[[1]])))) 
        stop("Different number of variables in each chain")
      xvarnames <- lapply(x, varnames, allow.null = FALSE)
      if (!all(unlist(lapply(xvarnames, "==", xvarnames[[1]])))) 
        stop("Different variable names in each chain")
    }
    if (is.R()) 
      class(x) <- "mcmc.list"
    else oldClass(x) <- "mcmc.list"
    return(x)
  }
  
  is.mcmc =   function (x)  # function copied from coda package, otherwise doparallel does not see it
  {
    if (inherits(x, "mcmc")) 
      if (length(dim(x)) == 3) 
        stop("Obsolete mcmc object\nUpdate with a command like\nx <- mcmcUpgrade(x)")
    else TRUE
    else FALSE
  }
  
  mcmc=function (data = NA, start = 1, end = numeric(0), thin = 1) 
  # function copied from coda package, otherwise doparallel does not see it 
  {
    if (is.matrix(data)) {
        niter <- nrow(data)
        nvar <- ncol(data)
    }
    else if (is.data.frame(data)) {
        if (!all(sapply(data, is.numeric))) {
            stop("Data frame contains non-numeric values")
        }
        data <- as.matrix(data)
        niter <- nrow(data)
        nvar <- ncol(data)
    }
    else {
        niter <- length(data)
        nvar <- 1
    }
    thin <- round(thin)
    if (length(start) > 1) 
        stop("Invalid start")
    if (length(end) > 1) 
        stop("Invalid end")
    if (length(thin) != 1) 
        stop("Invalid thin")
    if (missing(end)) 
        end <- start + (niter - 1) * thin
    else if (missing(start)) 
        start <- end - (niter - 1) * thin
    nobs <- floor((end - start)/thin + 1)
    if (niter < nobs) 
        stop("Start, end and thin incompatible with data")
    else {
        end <- start + thin * (nobs - 1)
        if (nobs < niter) 
            data <- data[1:nobs, , drop = FALSE]
    }
    attr(data, "mcpar") <- c(start, end, thin)
    attr(data, "class") <- "mcmc"
    data
  }
  
  nseeds = length(seed.index)
  table01 = permutations(2,nseeds,c(1,0),repeats=TRUE)
  table01 = table01[-which(rowSums(table01)<2),]
  sumindex = as.numeric(apply(table01,1,paste,collapse=""))
  R = length(datafile[[1]]$LL_total)
  LMD= c(rep(0,nseeds))
  listfromdata = list()
  #r=1
  listfromdata = (
    foreach(i=seed.index) %dopar%  # parallel, no .combine as I want a list, not a table or vector
     { out  =  mcmc (data = cbind(
          #datafile[[i]]$paramdraw[,1:5], 
          datafile[[i]]$alphadraw[-1:-burnin,,1], datafile[[i]]$betaMdraw[-1:-burnin,],
          datafile[[i]]$alphadraw[-1:-burnin,,2], datafile[[i]]$gammabetaSdraw[-1:-burnin,],  # these lines are for Gibbs sampler
          datafile[[i]]$lambdadraw[-1:-burnin,] , datafile[[i]]$sigma2mdraw[-1:-burnin,],    datafile[[i]]$sigma2ydraw[-1:-burnin,]
          #datafile[[i]]$rhodraw[-1:-burnin] ,        datafile[[i]]$sigma2mdraw[-1:-burnin,],    datafile[[i]]$sigma2ydraw[-1:-burnin,]
          ), start = burnin+1, end = R, thin = 1 )
       out
     }
  )
  for(i in seed.index)
  { LMD[i] = logMargDenNR(datafile[[i]]$LL_total[-1:-burnin]) }
  #for(i in seed.index)
  #{  listfromdata[[r]] = mcmc (data = cbind(
  #      #datafile[[i]]$paramdraw[,1:5], 
  #      datafile[[i]]$alphadraw[-1:-burnin,,1], datafile[[i]]$betaMdraw[-1:-burnin,],
  #      datafile[[i]]$alphadraw[-1:-burnin,,2], datafile[[i]]$gammabetaSdraw[-1:-burnin,],  # these lines are for Gibbs sampler
  #      datafile[[i]]$lambdadraw[-1:-burnin,] , datafile[[i]]$sigma2mdraw[-1:-burnin,],    datafile[[i]]$sigma2ydraw[-1:-burnin,]
  #      #datafile[[i]]$rhodraw[-1:-burnin] ,        datafile[[i]]$sigma2mdraw[-1:-burnin,],    datafile[[i]]$sigma2ydraw[-1:-burnin,]
  #      ), start = burnin+1, end = R, thin = 1 )
  #   LMD[r] = logMargDenNR(datafile[[i]]$LL_total[-1:-burnin])
  #   r=r+1
  #}
  listfromdata=listfromdata[!sapply(listfromdata, is.null)] 
  
  Psrf = matrix(0,nrow=nrow(table01),ncol=ncol(listfromdata[[1]]))
  Mpsrf = matrix(0,nrow=nrow(table01),ncol=1)
  for(j in 1:nrow(table01)) 
      { index1 = c(which(table01[j,]==1))
        j_listfromdata = listfromdata[index1]
        list_forRhat = mcmc.list(j_listfromdata)
        Psrf[j,] = gelman.diag(list_forRhat,autoburnin=FALSE)[[1]][,1]
        Mpsrf[j,1] = gelman.diag(list_forRhat,autoburnin=FALSE)[[2]]
      }
 
  #Tried to make parallel and did not work. Had to copy functions from code and still can't see mcpar
  #ncol_psrt_mpsrt = ncol(listfromdata[[1]])+1
  #Psrf_Mpsrt = list()
  #Psrf = matrix(0,nrow=nrow(table01),ncol=ncol(listfromdata[[1]]))
  #Mpsrf = matrix(0,nrow=nrow(table01),ncol=1)
  #Psrf_Mpsrt = (
  #  foreach(j=1:nrow(table01)) %dopar%  # parallel, no .combine=c as I need a list, not a table or vector as mcmc.list does not seem to work otherwise
  #  {#for(j in 1:nrow(table01) )   # need to add parallel computing
  #   index1 = c(which(table01[j,]==1))
  #   j_listfromdata = listfromdata[index1]
  #   list_forRhat = mcmc.list(j_listfromdata)
  #  pout = gelman.diag(list_forRhat,autoburnin=FALSE)
  #   pout
  #  }
  #)
  #for(j in 1:nrow(table01))
  #    { #Psrf[j,] = gelman.diag(list_forRhat,autoburnin=FALSE)[[1]][,1]
  #      #Mpsrf[j,1] = gelman.diag(list_forRhat,autoburnin=FALSE)[[2]]
  #      Psrf[j,] = Psrf_Mpsrt[j,-ncol_psrt_mpsrt]
  #      Mpsrf[j,1] = Psrf_Mpsrt[j,ncol_psrt_mpsrt]
  #}
  
  temp = which(Mpsrf<RhatCutoff)
  Rhat_sol1 = temp[which.max(rowSums(table01[temp,]))]   # solution 1, gives the row number in table01 to pick the seeds with solution 1
  index_sol1 = which(table01[Rhat_sol1,]==1)             # gives index of seeds with solution 1
  RhatEstSol1 = Psrf[Rhat_sol1,]
  MVRhatEsSol1 = Mpsrf[Rhat_sol1,1]
  LMD_Max_Sol1 = max(LMD[index_sol1])
  LMD_mean_Sol1 = mean(LMD[index_sol1])
  LMD_seedMax_Sol1 = index_sol1[which.max(LMD[index_sol1])]
  
  if(max(rowSums(table01[temp,]))<(nseeds)){ 
    index_sol2 = c(which(table01[Rhat_sol1,]==0))          # gives index of seeds with solution 2
    if(length(index_sol2)>1)
    { temp2 = rep(0,nseeds)
      for(t in index_sol2)  { temp2[t]=1 }
      Rhat_sol2 = temp[which(temp==which(sumindex==as.numeric(paste(temp2,collapse=""))))]  # solution 2, gives the row number in table01 to pick the seeds with solution 2
      RhatEstSol2 = Psrf[Rhat_sol2,]
      MVRhatEsSol2 = Mpsrf[Rhat_sol2,1]
      LMD_Max_Sol2 = max(LMD[index_sol2])
      LMD_seedMax_Sol2 = index_sol2[which.max(LMD[index_sol2])]
      LMD_mean_Sol2 = mean(LMD[index_sol2])
    }
    else {
      Rhat_sol2 = index_sol2
      RhatEstSol2 = Psrf[Rhat_sol2,]
      MVRhatEsSol2 = Mpsrf[Rhat_sol2,1]
      LMD_Max_Sol2 = max(LMD[index_sol2])
      LMD_seedMax_Sol2 = index_sol2
      LMD_mean_Sol2 = mean(LMD[index_sol2])
    }
  }
  else {  
    Rhat_sol2 = "NA" 
    RhatEstSol2 = "NA"
    MVRhatEsSol2 = "NA"
    if(max(rowSums(table01[temp,]))==(nseeds-1)) 
    {  index_sol2 = c(which(table01[Rhat_sol1,]==0))
       LMD_Max_Sol2 = max(LMD[index_sol2])
       LMD_seedMax_Sol2 = index_sol2[which.max(LMD[index_sol2])]
       LMD_mean_Sol2 = mean(LMD[index_sol2])
    }
    else {
       index_sol2 = "NA"
       LMD_Max_Sol2 = "NA"
       LMD_seedMax_Sol2 = "NA"
       LMD_mean_Sol2 = "NA"
    }
  } 
  return(list(RhatEstAll = Psrf, MVRhatEstAll = Mpsrf, index_sol1=index_sol1, index_sol2=index_sol2,
              Rhat_sol1=Rhat_sol1,Rhat_sol2=Rhat_sol2,
              RhatEstSol1 = RhatEstSol1,  RhatEstSol2 = RhatEstSol2, MVRhatEsSol1 = MVRhatEsSol1, MVRhatEsSol2 = MVRhatEsSol2, MpsrfAllSeeds =  Mpsrf[nrow(table01)],
              LMD=LMD, LMD_mean = mean(LMD),  LMD_Max = max(LMD), LMD_seedMax = which.max(LMD), 
              LMD_Max_Sol1 = LMD_Max_Sol1,LMD_Max_Sol2 = LMD_Max_Sol2,LMD_seedMax_Sol1 = LMD_seedMax_Sol1,LMD_seedMax_Sol2 =LMD_seedMax_Sol2,
              LMD_mean_Sol1 = LMD_mean_Sol1,   LMD_mean_Sol2 = LMD_mean_Sol2,seed.index=seed.index,burnin=burnin) )
}



#  NOT DEVELOPED FOR MORE THAN 2 SOLUTIONS !!!

library(bayesm);   library(coda); library(gtools)
FUN_Mediation_LMD_RHat_MS_MH = function(datafile,seed.index,burnin,RhatCutoff)
{ 
  nseeds = length(seed.index)
  table01 = permutations(2,nseeds,c(1,0),repeats=TRUE)
  table01 = table01[-which(rowSums(table01)<2),]
  sumindex = as.numeric(apply(table01,1,paste,collapse=""))
  R = length(datafile[[1]]$rhodraw)
  LMD= c(rep(0,nseeds))
  listfromdata = list()
  r=1
  for(i in seed.index)
  {  listfromdata[[r]] = mcmc (data = cbind(
        datafile[[i]]$paramdraw[,1:5], 
        datafile[[i]]$rhodraw  #,        datafile[[i]]$sigma2mdraw,    datafile[[i]]$sigma2ydraw
        ), start = burnin+1, end = R, thin = 1 )
     LMD[r] = logMargDenNR(datafile[[i]]$LL_total[-1:-burnin])  
     r=r+1
  }
  listfromdata=listfromdata[!sapply(listfromdata, is.null)] 
  
  Psrf = matrix(0,nrow=nrow(table01),ncol=ncol(listfromdata[[1]]))
  Mpsrf = matrix(0,nrow=nrow(table01),ncol=1)
  for(j in 1:nrow(table01) )
  {  index1 = c(which(table01[j,]==1))
     j_listfromdata = listfromdata[index1]
     list_forRhat = mcmc.list(j_listfromdata)
     Psrf[j,] = gelman.diag(list_forRhat,autoburnin=FALSE)[[1]][,1]
     #Mpsrf[j,1] = gelman.diag(list_forRhat,autoburnin=FALSE)[[2]]
  }
  
  temp = which(Mpsrf<RhatCutoff)
  Rhat_sol1 = temp[which.max(rowSums(table01[temp,]))]   # solution 1, gives the row number in table01 to pick the seeds with solution 1
  index_sol1 = which(table01[Rhat_sol1,]==1)             # gives index of seeds with solution 1
  RhatEstSol1 = Psrf[Rhat_sol1,]
  MVRhatEsSol1 = Mpsrf[Rhat_sol1]
  LMD_Max_Sol1 = max(LMD[index_sol1])
  LMD_mean_Sol1 = mean(LMD[index_sol1])
  LMD_seedMax_Sol1 = index_sol1[which.max(LMD[index_sol1])]
  
  if(max(rowSums(table01[temp,]))<(nseeds-1)){ 
    index_sol2 = c(which(table01[Rhat_sol1,]==0))          # gives index of seeds with solution 2
    temp2 = rep(0,nseeds)
    for (t in index_sol2)  { temp2[t]=1}
    Rhat_sol2 = temp[which(temp==which(sumindex==as.numeric(paste(temp2,collapse=""))))]  # solution 2, gives the row number in table01 to pick the seeds with solution 2
    RhatEstSol2 = Psrf[Rhat_sol2,]
    MVRhatEsSol2 = Mpsrf[Rhat_sol2]
    LMD_Max_Sol2 = max(LMD[index_sol2])
    LMD_seedMax_Sol2 = index_sol2[which.max(LMD[index_sol2])]
    LMD_mean_Sol2 = mean(LMD[index_sol2])
  }
  else {  
    Rhat_sol2 = "NA" 
    RhatEstSol2 = "NA"
    MVRhatEsSol2 = "NA"
    Rhat_sol2 = "NA"
    if(max(rowSums(table01[temp,]))==(nseeds-1)) 
    {  index_sol2 = c(which(table01[Rhat_sol1,]==0))
       LMD_Max_Sol2 = max(LMD[index_sol2])
       LMD_seedMax_Sol2 = index_sol2[which.max(LMD[index_sol2])]
       LMD_mean_Sol2 = mean(LMD[index_sol2])
       LMD_Max_Sol2 = max(LMD[index_sol2])
    }
    else {
       index_sol2 = "NA"
       LMD_Max_Sol2 = "NA"
       LMD_seedMax_Sol2 = "NA"
       LMD_mean_Sol2 = "NA"
       LMD_Max_Sol2 = "NA"       
    }
  } 
  return(list(RhatEstAll = Psrf, MVRhatEstAll = Mpsrf, index_sol1=index_sol1, index_sol2=index_sol2,Rhat_sol1=Rhat_sol1,Rhat_sol2=Rhat_sol2,
              RhatEstSol1 = RhatEstSol1,  RhatEstSol2 = RhatEstSol2, MVRhatEsSol1 = MVRhatEsSol1, MVRhatEsSol2 = MVRhatEsSol2, MpsrfAllSeeds =  Mpsrf[nrow(table01)],
              LMD=LMD, LMD_mean = mean(LMD),  LMD_Max = max(LMD), LMD_seedMax = which.max(LMD), 
              LMD_Max_Sol1 = LMD_Max_Sol1,LMD_Max_Sol2 = LMD_Max_Sol2,LMD_seedMax_Sol1 = LMD_seedMax_Sol1,LMD_seedMax_Sol2 =LMD_seedMax_Sol2,
              LMD_mean_Sol1 = LMD_mean_Sol1,   LMD_mean_Sol2 = LMD_mean_Sol2,seed.index=seed.index,burnin=burnin) )
}