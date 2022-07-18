
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



#------(Aggregate model)----------------------#
# function to run aggregate model: returns a list of results (8 variables)


#---- (Aggregate model Results Page) --------------------------#
# Proportion of joint (alpha, beta) in each quadrant MULTIX
#TRANSFERRED TO output_helper (7-13-2021)

# FUN_PDF_Mediation_AlphaBetaProportion_Aggregate_forShiny  = function(filename, burnin, x_vars)
#   { 
#   nvarX = ncol(filename$alphadraw) - 1
#   QuadrantsCounts = array(0, dim=c(4, nvarX))
#   DrawsAnalysis   = c(seq(from = burnin+1, to = length(filename$LL_total), by = 1))
#   pp=pn=np=nn=0
#   for(j in 2:(nvarX+1)){  
#     for(r in DrawsAnalysis){ 
#       pp = pp + ifelse(filename$alphadraw[r,j,1]>0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
#       pn = pn + ifelse(filename$alphadraw[r,j,1]>0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
#       np = np + ifelse(filename$alphadraw[r,j,1]<0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
#       nn = nn + ifelse(filename$alphadraw[r,j,1]<0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
#     }
#     QuadrantsCounts[,(j-1)]=c(pp,pn,np,nn)
#     pp=pn=np=nn=0
#   }
#   Proportions <- QuadrantsCounts/length(DrawsAnalysis)
#   colnames(Proportions) <- x_vars
#   rownames(Proportions) <- c("I (++)","II (+-)","III (-+)","IV (--)")
#   # return(list(Proportions = QuadrantsCounts/length(DrawsAnalysis)))
#   return(Proportions)
# }


#---- (Aggregate model Results Page) --------------------------#
# HDP MULTIX
#TRANSFERED TO output_helper 7-13-2021

# FUN_PDF_Mediation_HDPI_Aggregate_forShiny  = function(filename, burnin, x_vars)
#   { nvarX = ncol(filename$alphadraw)
#   outputTable = NULL
#   for(j in 1:(nvarX)){  
#     tempa = hdi(filename$alphadraw[-1:-burnin,j,1], credMass = 0.95)
#     outputTable = rbind(outputTable,c(mean(filename$alphadraw[-1:-burnin,j,1]),tempa))
#   }
#   tempb = hdi(filename$betaMdraw[-1:-burnin,1], credMass = 0.95)
#   outputTable = rbind(outputTable,c(mean(filename$betaMdraw[-1:-burnin,1]),tempb))
#   for(j in 1:(nvarX)){  
#     tempg = hdi(filename$gammabetaSdraw[-1:-burnin,j], credMass = 0.95)
#     outputTable = rbind(outputTable,c(mean(filename$gammabetaSdraw[-1:-burnin,j]),tempg))
#   }
#   
#   #data.frame(outputTable)
#   
#   rownames_list = c(rep(0,nrow(outputTable)))
#   rownames_list[1] = "alpha_{0}"
#   rownames_list[nvarX+1] = "beta"
#   rownames_list[nvarX+2] = "gamma_{0}"
#   for(i in 2:(nvarX)) {
#     rownames_list[i]=paste0("alpha_{", i - 1, "}")
#     rownames_list[nvarX+1+i]=paste0("gamma_{", i - 1, "}")
#   }
#   #outputTable = cbind(rownames_list,outputTable)
#   
#   colnames(outputTable) <- c("Mean","Lower limit","Upper limit")
#   
#   # rownames(outputTable) <- rownames_list
#   # colnames(outputTable) <- c("Mean","Lower limit","Upper limit")
#   # rownames_list = c(rep(0,nrow(outputTable)))
#   # rownames_list[1] = bquote(alpha[0]) 
#   # rownames_list[nvarX+1] = paste(expression(beta))
#   # rownames_list[nvarX+2] = paste(expression(gamma[0]))
#   # for(i in 2:(nvarX)) {
#   #   rownames_list[i]=paste(expression(alpha[i-1]))
#   #   rownames_list[nvarX+1+i]=paste(expression(gamma[i-1]))
#   # }
#   rownames(outputTable) <- rownames_list
#   # return(list(Proportions = QuadrantsCounts/length(DrawsAnalysis)))
#   #return(grid.table(outputTable))
#   return(outputTable)
#   #return(htmlTable(outputTable))
# }


#-----(Aggregate model)-----------------------#
# function to generate PDF scatterplots of postreior draws of alphas and beta for aggregate model

#FUN_PDF_Mediation_ScatterPlots_Aggregate_forShiny_Plot = function(dataset, filename, burnin, x_vars) {
#    # dataset = the name (string) of the data file that needs to be displayed (i.e. "Loyalty data")
#    # filename = output file of MCMC (list of objects/tables)
#    # burnin = number of saved draws to exclude from plotting
#    # x_vars = vector of names of variables in X, excluding intercept
#  nvarX = ncol(filename$alphadraw) - 1
#  #pdf(paste(dataset, "Posterior Draws.pdf", sep = " "), width=pdfW, height=pdfH)
#  par(oma=c(0,0,2,0));
#  QuadrantsCounts = array(0, dim=c(4, nvarX))
#  #print(burnin)
#  #print(length(filename$LL_total))
#  DrawsAnalysis   = c(seq(from = burnin+1, to = length(filename$LL_total), by = 1))
#  pp=pn=np=nn=0
#  for(j in 1:nvarX){  
#    for(r in DrawsAnalysis){ 
#      pp = pp + ifelse(filename$alphadraw[r,j+1,1]>0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
#      pn = pn + ifelse(filename$alphadraw[r,j+1,1]>0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
#      np = np + ifelse(filename$alphadraw[r,j+1,1]<0,ifelse(filename$betaMdraw[r,1]>0,1,0),0)
#      nn = nn + ifelse(filename$alphadraw[r,j+1,1]<0,ifelse(filename$betaMdraw[r,1]<0,1,0),0)
#    }
#    QuadrantsCounts[,j]=c(pp,pn,np,nn)
#    pp=pn=np=nn=0
#  }
#  Proportions <- QuadrantsCounts/length(DrawsAnalysis)
#  
#  par(mfrow=c(nvarX,2))
#  ylimB = c( min(filename$betaMdraw[-1:-burnin])-0.1*min(filename$betaMdraw[-1:-burnin,1]),
#           max(filename$betaMdraw[-1:-burnin])+0.1*max(filename$betaMdraw[-1:-burnin,1]))
#  for(j in 1:nvarX){
#    ylimA = c( min(filename$alphadraw[-1:-burnin,j+1,1])-0.1*min(filename$alphadraw[-1:-burnin,j+1,1]),
#             max(filename$alphadraw[-1:-burnin,j+1,1])+0.1*max(filename$alphadraw[-1:-burnin,j+1,1]))
#    ylimG = c( min(filename$gammabetaSdraw[-1:-burnin,j+1])-0.1*min(filename$gammabetaSdraw[-1:-burnin,j+1]),
#             max(filename$gammabetaSdraw[-1:-burnin,j+1])+0.1*max(filename$gammabetaSdraw[-1:-burnin,j+1]))
#    breaksCalc = (ylimG[2]-ylimG[1])/25
#    plot(filename$alphadraw[-1:-burnin,j+1,1], filename$betaMdraw[-1:-burnin,1], 
#         main=bquote("Scatterplot of " ~ alpha[.(j)] ~ " and " ~ beta ~ " for variable " ~ .(x_vars[j])),
#         xlab=bquote(alpha[.(j)]), ylab=expression(beta), xlim=ylimA,ylim=ylimB)
#    abline(h=0,v=0,col="gray")
#    if(ylimA[2]>0 & ylimB[2]>0){
#      text(x=c(ylimA[2]),  y=c(ylimB[2]), 
#         labels = c(paste( "I (",round(Proportions[1,j],4),")")), pos=2, font=2,cex=0.9)
#    }
#    if(ylimA[2]>0 & ylimB[1]<0){
#      text(x=c(ylimA[2]),  y=c(ylimB[1]), 
#         labels = c(paste( "II (",round(Proportions[2,j],4),")")), pos=2, font=2,cex=0.9)
#    }
#    if(ylimA[1]<0 & ylimB[1]<0){
#      text(x=c(ylimA[1]),  y=c(ylimB[1]), 
#         labels = c(paste("IV (",round(Proportions[4,j],4),")")), pos=4, font=2,cex=0.9)
#    }
#    if(ylimA[1]<0 & ylimB[2]>0){
#      text(x=c(ylimA[1]),  y=c(ylimB[2]), 
#         labels = c(paste("IV (",round(Proportions[3,j],4),")")), pos=4, font=2,cex=0.9)
#    }
#    hist(filename$gammabetaSdraw[-1:-burnin,j+1], main=bquote(paste("Histogram of ",gamma[.(j)])),
#         xlab=bquote(gamma[.(j)]), xlim=ylimG,
#         breaks=c(seq(min(filename$gammabetaSdraw[-1:-burnin,j+1]),(max(filename$gammabetaSdraw[-1:-burnin,j+1])+breaksCalc),breaksCalc )))
#    abline(v=0,col="darkred", lwd=2)
#  }
#  #title(main=paste(dataset, "Aggregate Model.  Posterior Draws of parameters."),outer=T)
## dev.off()
#}
