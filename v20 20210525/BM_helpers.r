


#---- (Mixture MS) Proportion of joint (alpha, beta) in each quadrant for both M and S segments MULTIX ---------------------------------------------------------------------------------
# this ONLY runs for ONE BEST seed, which will be taken from the results of RHat function

FUN_PDF_Mediation_AlphaBetaProportion_MSmixture_forShiny = function(filenamelist,seed.selected,x_vars,burnin)
  # seed.selected is ONE selected seed for analysis, not all seeds
{  
   nvarX = ncol(filenamelist[[1]]$alphadraw)
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
  ProportionsM <- round(QuadrantsCountsM/length(DrawsAnalysis),4)
  ProportionsS <- round(QuadrantsCountsS/length(DrawsAnalysis),4)
  Proportions <- cbind(ProportionsM,ProportionsS)
  
  temp <- c(rep(0,(nvarX-1)*2))
  for(n in 2:nvarX)
  {
     temp[n-1] <-paste(x_vars[n - 1], " M",sep = "")
     temp[(nvarX-1)+n-1] <- paste(x_vars[n - 1], " G",sep = "")
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
# used to output PDF plots of the effects

FUN_PDF_Mediation_FinalPlots_MSmixture_forShiny_Plot = function(dataset,filenamelist,seed.selected,burnin){
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
             main=paste("Segment M. Seed=",seed.selected[i]),
           xlab=bquote(alpha[M][.(p-1)]),ylab=expression(beta[M]), xlim=ylimA,ylim=ylimB)
        abline(h=0,v=0,col="gray")
      }
      for(p in 2:nvarX) {   
        plot(filename$alphadraw[-1:-burnin,p,2], filename$gammabetaSdraw[-1:-burnin,nvarX+1], main=paste("Segment G. Seed=",seed.selected[i]),
           xlab=bquote(alpha[S][.(p-1)]),ylab=expression(beta[S]), xlim=ylimA,ylim=ylimB)
        abline(h=0,v=0,col="gray")
      }
      for(p in 2:nvarX){
         hist(filename$gammabetaSdraw[-1:-burnin,p], xlab=bquote(gamma[S][.(p-1)]), xlim=ylimG, main=paste("Segment G. Seed=",seed.selected[i]),
           breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,p]),(max(filename$gammabetaSdraw[-1:-burnin,p])+breaksCalc),breaksCalc )))
         abline(v=0,col="darkred",lwd=2)
      }
      hist(filename$rhodraw[-1:-burnin], main=paste("Rho. Seed=",seed.selected[i]), xlab=bquote(rho), xlim=c(0,1),
           breaks=50)
    }
  #title(main=paste(dataset," Final Results. Binary mixture"),outer=T)
  #dev.off()
}

#-------------- BM model scatterplots (updated version for the App)--------------------
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
         #main=bquote("Scatterplot of " ~ alpha[.(p-1)] ~ " and " ~ beta ~ " for variable " ~ .(x_var[p-1]) ~ ". Segment M"),
         main=bquote("Segment M"),
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
         #main=bquote("Scatterplot of " ~ alpha[.(p-1)] ~ " and " ~ beta ~ ". Variable " ~ .(x_var[p-1]) ~ ". Segment S"),
         main=bquote("Segment S"),
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
         #main=bquote("Histogram of " ~ gamma[.(p-1)] ~ ". Variable " ~ .(x_var[p-1]) ~ ". Segment S"),
         main=bquote("Segment S"),
         #main=paste(x_var[p-1],"\n Segment S"),
         xlab=bquote(gamma[S][.(p-1)]), xlim=ylimG,
         breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,p]),(max(filename$gammabetaSdraw[-1:-burnin,p])+breaksCalc),breaksCalc )))
    abline(v=0,col="darkred",lwd=2)
  }
 
}

# testing
#debug(FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Effects)
#FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Effects(dataset="",filenamelist=output_list_b,seed.list,seed.selected=c(1),burnin=10)


FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Rho = function(dataset,filenamelist,seed.list,seed.selected,burnin){
  #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
  filename = filenamelist[[seed.selected]]
  hist((filename$rhodraw[-1:-burnin,]), main="", xlab=bquote(rho), xlim=c(0,1),col = "#75AADB",
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

FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Lambda = function(dataset,filenamelist,seed.list,seed.selected,burnin){
  #pdf(paste(dataset,"_FinalPlots.pdf", sep = ""), width=pdfW, height=pdfH)
  filename = filenamelist[[seed.selected]]
  numCharts = ncol(filename$lambdadraw)
  par(oma=c(0,0,2,0));   par(mfrow=c(numCharts,1))
  for(i in 1:numCharts)
  {  hist(filename$lambdadraw[-1:-burnin,i], main="", xlab=bquote(lambda[i-1]),
          xlim=c(min(filename$lambdadraw[-1:-burnin,i]),max(filename$lambdadraw[-1:-burnin,i])),col = "#75AADB",
       breaks=30)
  }
}


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








