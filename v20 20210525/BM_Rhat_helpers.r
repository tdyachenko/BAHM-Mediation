
#  NOT DEVELOPED FOR MORE THAN 2 SOLUTIONS !!!

#library(bayesm);   library(coda); library(gtools)


# testing
#output_list_b = list(output_list_b[[1]],output_list_b[[1]],output_list_b[[1]],output_list_b[[1]])
#testtemp = FUN_Mediation_LMD_RHat_MS(datafile=output_list_b, seed.index=seq(1:4), burnin=10, RhatCutoff=1.05)

# this will need to be run before some of the other results for the binary model.. since we want to know 
# which seeds are best for the two soltuions (if there are two solutions).
# also want to print some of these results to screen

# Function is updated 2020-02-28 (TD) to accomodate moderate mediation case
# lambda is now the parameter which is the transformation of previous rho

# TO DO (TD - 20200228): automatically find MCMC chains that did not converge to one mode and eliminate them from analysis

FUN_Mediation_LMD_RHat_MS_cov = function(inputdata, datafile,seed.index,burnin,RhatCutoff)
{ 
  nseeds = length(seed.index)
  table01 = permutations(2,nseeds,c(1,0),repeats=TRUE)  # table of all possible combinations of the seeds as groups in two groups
  table01 = table01[-which(rowSums(table01)<2),]        # only select permutations with 2 or more members in each group
  table01 = table01[-which(rowSums(table01)==(nseeds-1)),]        # only select permutationa with 2 or more members in each group
  sumindex = as.numeric(apply(table01,1,paste,collapse=""))
  R = length(datafile[[1]]$wdraw[1,])
  LMD= c(rep(0,nseeds))
  listfromdata = list()  # generate list for gelman.diag function
  r=1
  for(i in seed.index)
  {  listfromdata[[r]] = mcmc (data = cbind(
        #datafile[[i]]$paramdraw[,1:5], 
        datafile[[i]]$alphadraw[-1:-burnin,,1],datafile[[i]]$betaMdraw[-1:-burnin,],
        datafile[[i]]$alphadraw[-1:-burnin,,2], datafile[[i]]$gammabetaSdraw[-1:-burnin,],  # these lines are for Gibbs sampler
        datafile[[i]]$lambdadraw[-1:-burnin,] , datafile[[i]]$sigma2mdraw[-1:-burnin,],    datafile[[i]]$sigma2ydraw[-1:-burnin,]
        #datafile[[i]]$rhodraw[-1:-burnin] ,        datafile[[i]]$sigma2mdraw[-1:-burnin,],    datafile[[i]]$sigma2ydraw[-1:-burnin,]
        ), start = burnin+1, end = R, thin = 1 )
     LMD[r] = logMargDenNR(datafile[[i]]$LL_total[-1:-burnin])
     r=r+1
  }
  listfromdata=listfromdata[!sapply(listfromdata, is.null)] 
  
  Psrf = matrix(0,nrow=nrow(table01),ncol=ncol(listfromdata[[1]]))   
       # nrow=number of possible combinations of seeds into two groups from table01 (possible groupings)
       # ncol= nvar, including Rho and sigmas
  Mpsrf = matrix(0,nrow=nrow(table01),ncol=1)
  for(j in 1:nrow(table01) )
  {  index1 = c(which(table01[j,]==1))  # select seeds that qre groups together, indexed 1 in table01
     j_listfromdata = listfromdata[index1]
     list_forRhat = mcmc.list(j_listfromdata)
     temp = gelman.diag(list_forRhat,autoburnin=FALSE)
     Psrf[j,] = temp[[1]][,1]  # point estimate only
     Mpsrf[j,1] = temp[[2]]
  }
  
  DIC = FUN_DIC_mediation(inputdata, datafile,burnin,ModelFlag=2)  #### ADDED, needs to be tested

  temp = which(Mpsrf<RhatCutoff)
  Rhat_sol1 = temp[which.max(rowSums(table01[temp,]))]   # solution 1, gives the row number in table01 to pick the seeds with solution 1
  index_sol1 = which(table01[Rhat_sol1,]==1)             # gives index of seeds with solution 1
  RhatEstSol1 = Psrf[Rhat_sol1,]
  MVRhatEsSol1 = Mpsrf[Rhat_sol1]
  LMD_Max_Sol1 = max(LMD[index_sol1])
  LMD_mean_Sol1 = mean(LMD[index_sol1])
  LMD_seedMax_Sol1 = index_sol1[which.max(LMD[index_sol1])]
  DIC_Max_Sol1 = max(DIC[index_sol1])
  DIC_mean_Sol1 = mean(DIC[index_sol1])
  DIC_seedMax_Sol1 = index_sol1[which.max(DIC[index_sol1])]
  
  if(max(rowSums(table01[temp,]))<(nseeds)){ 
    index_sol2 = c(which(table01[Rhat_sol1,]==0))          # gives index of seeds with solution 2
    if(length(index_sol2)>1)
    { temp2 = rep(0,nseeds)
      for(t in index_sol2)  { temp2[t]=1 }
      Rhat_sol2 = temp[which(temp==which(sumindex==as.numeric(paste(temp2,collapse=""))))]  # solution 2, gives the row number in table01 to pick the seeds with solution 2
      RhatEstSol2 = Psrf[Rhat_sol2,]
      MVRhatEsSol2 = Mpsrf[Rhat_sol2]
      LMD_Max_Sol2 = max(LMD[index_sol2])
      LMD_seedMax_Sol2 = index_sol2[which.max(LMD[index_sol2])]
      LMD_mean_Sol2 = mean(LMD[index_sol2])
      DIC_Max_Sol2 = max(DIC[index_sol2])
      DIC_mean_Sol2 = mean(DIC[index_sol2])
      DIC_seedMax_Sol2 = index_sol1[which.max(DIC[index_sol2])]
    }
    else {
      Rhat_sol2 = index_sol2
      RhatEstSol2 = Psrf[index_sol2,]
      MVRhatEsSol2 = Mpsrf[index_sol2]
      LMD_Max_Sol2 = max(LMD[index_sol2])
      LMD_mean_Sol2 = mean(LMD[index_sol2])
      LMD_seedMax_Sol2 = index_sol2
      DIC_Max_Sol2 = max(DIC[index_sol2])
      DIC_mean_Sol2 = mean(DIC[index_sol2])
      DIC_seedMax_Sol2 = index_sol2
    }
  }
  else {  
    Rhat_sol2 = 999 
    RhatEstSol2 = 999
    MVRhatEsSol2 = 999
    if(max(rowSums(table01[temp,]))==(nseeds-1)) 
    {  index_sol2 = c(which(table01[Rhat_sol1,]==0))
       LMD_Max_Sol2 = max(LMD[index_sol2])
       LMD_mean_Sol2 = mean(LMD[index_sol2])
       LMD_seedMax_Sol2 = index_sol2[which.max(LMD[index_sol2])]
       DIC_Max_Sol2 = max(DIC[index_sol2])
       DIC_mean_Sol2 = mean(DIC[index_sol2])
       DIC_seedMax_Sol2 = index_sol2[which.max(DoC[index_sol2])]
    }
    else {
       index_sol2 = 999
       LMD_Max_Sol2 = 999
       LMD_seedMax_Sol2 = 999
       LMD_mean_Sol2 = 999
       DIC_Max_Sol2 = 999
       DIC_mean_Sol2 = 999
       DIC_seedMax_Sol2 = 999
    }
  } 
  
  if(sum(is.finite(RhatEstSol2))>0)
   {      Rhat = cbind(matrix(RhatEstSol1,ncol=1), matrix(RhatEstSol2,ncol=1)) }
   else { Rhat = cbind(matrix(RhatEstSol1,ncol=1), rep(999,length(RhatEstSol1)))    }

  table_forShiny = matrix(0,ncol=2,nrow=8)
  table_forShiny[1,] = c(format(LMD_seedMax_Sol1,nsmall=0),format(LMD_seedMax_Sol2,nsmall=0))
  table_forShiny[2,1] = ifelse(length(index_sol1)>0,toString(index_sol1),999)
  table_forShiny[2,2] = ifelse(length(index_sol1)>0,toString(index_sol2),999)
  table_forShiny[3,1] = round(mean(Rhat[,1],na.rm = TRUE),digits=2)
  table_forShiny[3,2] = ifelse(mean(Rhat[,2],na.rm = TRUE)>0,round(mean(Rhat[,2],na.rm = TRUE),digits=2),999)
  table_forShiny[4,] = c(round(MVRhatEsSol1,digits=2),    round(MVRhatEsSol2,digits=2))
  table_forShiny[5,] = c(round(LMD_mean_Sol1,digits=1),   round(LMD_mean_Sol2,digits=1))
  table_forShiny[6,] = c(round(LMD_Max_Sol1,digits=1),    round(LMD_Max_Sol2,digits=1))
  table_forShiny[7,] = c(round(DIC_mean_Sol1,digits=1),   round(DIC_mean_Sol2,digits=1))
  table_forShiny[8,] = c(round(DIC_Max_Sol1,digits=1),    round(DIC_Max_Sol2,digits=1))
  
  colnames(table_forShiny) = c("Solution 1","Solution 2")
  rownames(table_forShiny) = c("Seed number selected",
                               "Seeds in this group",
                               "Mean of the variable potential scale reduction factors' (psrf) point est.",       #Psrf point est.
                               "Multivariate potential scale reduction factor (mpsrf) point est.",      #Mpsrf point est.
                               "Mean LMD NR",
                               "Max LMD NR",
                               "Mean DIC",
                               "Max DIC")
  
 # if(LMD_Max_Sol1>LMD_Max_Sol2){
    colnames(Rhat) = c("Solution 1","Solution 2")
    output = list(table_forShiny = table_forShiny,
                  RhatEstAll = Psrf,
                  RhatEst = Rhat,
                  MVRhatEstAll = Mpsrf,
                  index_sol1=index_sol1,
                  index_sol2=index_sol2,
                  RhatEstSol1 = RhatEstSol1,
                  RhatEstSol2 = RhatEstSol2, 
                  MpsrfAllSeeds =  Mpsrf[nrow(table01)],
                  LMD=LMD,
                  LMD_mean = mean(LMD),
                  LMD_Max = max(LMD),
                  LMD_seedMax = which.max(LMD), 
                  DIC_mean = mean(DIC),
                  DIC_Max = max(DIC),
                  DIC_seedMax = which.max(DIC),
                  seed.index=seed.index,
                  burnin=burnin)
 # }
 # else{
 #   table_forShiny = cbind(table_forShiny[,2],table_forShiny[,1])
 #   colnames(table_forShiny) = c("Solution 1","Solution 2")
 #   Rhat = cbind(Rhat[,2],Rhat[,1])
 #   colnames(Rhat) = c("Solution 1","Solution 2")
 #   output = list(table_forShiny = table_forShiny,
 #                 RhatEstAll = Psrf,
 #                 RhatEst = Rhat,
 #                 MVRhatEstAll = Mpsrf,
 #                 index_sol1=index_sol2,
 #                 index_sol2=index_sol1,
 #                 RhatEstSol1 = RhatEstSol2,
 #                 RhatEstSol2 = RhatEstSol1, 
 #                 MpsrfAllSeeds =  Mpsrf[nrow(table01)],
 #                 LMD=LMD,
 #                 LMD_mean = mean(LMD),
 #                 LMD_Max = max(LMD),
 #                 LMD_seedMax = which.max(LMD), 
 #                 DIC_mean = mean(DIC),
 #                 DIC_Max = max(DIC),
 #                 DIC_seedMax = which.max(DIC),
 #                 seed.index=seed.index,
 #                 burnin=burnin)
 # }
  return( output  )
}



# --- Generate the best seed ----------------------------------------------------
FUN_BestSeed_ForShiny  = function(filename)
  { BestSeed <- as.numeric(filename$table_forShiny[1,1])
    return(BestSeed)
}

#FUN_BestSeed_ForShiny(testtemp)

