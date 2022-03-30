
#-------------------- run model: 1 loop ------------------#


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

#---- (Mixture MS)  MCMC plots  MULTIX ----------------------------------------------------------------------------
# note: this will just be available as pdf download (no need to print to screen)
# this can be run without the RHat results
# title: MCMC plots
#  NOT DONE???

# TRANSFERRED TO output_helper (7-21-2021 TD)
# FUN_PDF_MCMC_Mediation_forShiny = function(dataset,filenamelist,seed.index,seed.list,burnin)

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
         main=bquote("Scatterplot of " ~ alpha[.(p-1)] ~ " and " ~ beta ~ ". Variable " ~ .(x_var[p-1]) ~ ". Segment S"),
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
         main=bquote("Histogram of " ~ gamma[.(p-1)] ~ ". Variable " ~ .(x_var[p-1]) ~ ". Segment S"),
         #main=paste(x_var[p-1],"\n Segment S"),
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
  
  for(i in 1:numCharts)
  {  hist(filename$lambdadraw[-1:-burnin,i], main="", xlab=bquote(lambda[i]), xlim=c(min(filename$lambdadraw[-1:-burnin,i]),max(filename$lambdadraw[-1:-burnin,i])),col = "#75AADB",
       breaks=50)
  }
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
# outputs HPD intervals

FUN_PDF_Mediation_Parameters_MSmixture_forShiny  = function(filenamelist, x_vars, m_var, z_vars, seed.list, burnin)
  { filename = filenamelist[[seed.list]]
    nvarX = ncol(filename$alphadraw[-1:-burnin,,1])
    tempCIs_M = rbind(
       cbind(colMeans(filename$alphadraw[-1:-burnin,,1]), 
         t( apply(filename$alphadraw[-1:-burnin,,1],2,hdi,credMass = 0.95))),  # aM
       cbind(colMeans(filename$betaMdraw[-1:-burnin,]),
         t( apply(filename$betaMdraw[-1:-burnin,],2,hdi,credMass = 0.95)))     # bM
       #cbind(colMeans(filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),
         #t(apply((filename$alphadraw[-1:-burnin,,1]*filename$betaMdraw[-1:-burnin,2]),2, hdi,credMass = 0.95)))  # abM
    )
    tempCIs_S = rbind(
       cbind(colMeans(filename$alphadraw[-1:-burnin,,2]), 
         t( apply(filename$alphadraw[-1:-burnin,,2],2,hdi,credMass = 0.95))),    # aS
       cbind(colMeans(filename$gammabetaSdraw[-1:-burnin,]),
         t( apply(filename$gammabetaSdraw[-1:-burnin,],2,hdi,credMass = 0.95))) # gbS
       #cbind(colMeans(filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),
         #t(apply((filename$alphadraw[-1:-burnin,,2]*filename$gammabetaSdraw[-1:-burnin,3]),2, hdi,credMass = 0.95))) # abS
    )
    tempCIs_Rho =  matrix(c(mean(filename$rhodraw[-1:-burnin]),
                     t(hdi(filename$rhodraw[-1:-burnin], credMass = 0.95))),ncol=3,nrow=1)
    tempCIs_lambda = cbind(colMeans(filename$lambdadraw[-1:-burnin, , drop = FALSE]),
         t(apply(filename$lambdadraw[-1:-burnin, , drop = FALSE],2,hdi,credMass = .95)))
    #templambda = t(apply(filename$lambdadraw[-1:-burnin, , drop = FALSE],2,hdi,credMass = .95))
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








