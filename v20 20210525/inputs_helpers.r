
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

get_inputs_agg <-function(df, x_vars, y_var, m_var, z_vars,
                          Aa_var, Abg_var, nu_var, qm_var, qy_var,
                          R_var, seed_var, keep_var){
  
  df = as.matrix(df)
  
  X   = df[,x_vars]
  nhh = nrow(df)
  Z   = df[,z_vars]
  
  y = df[,y_var, drop=FALSE]
  m = df[,m_var, drop=FALSE]
  
  nvarX = 1 + length(x_vars) # number of X variables selected on the menu. intercept is added here
  nvarM = 1
  
  Data = list(
    X = cbind(rep(1, nhh), X), # datafile with only variables selected as X (independent variable(s))
    y = y, # datafile with only one variable selected as y (dependent variable)
    m = m,  # datafile with only one variable selected as m (mediator))
    Z = cbind(rep(1, nhh), Z)
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

#------(BM model)--------------------------------#
# prep inputs for BM model: returns a list of inputs for Data, Mcmc, Prior

get_inputs_binary <-function(df, x_vars, y_var, m_var, z_vars,
                             Aa_var, Abg_var, Al_var, nu_var, qy_var, qm_var,
                             R_var, keep_var, slambda1, slambda2)
{
  df  = as.matrix(df)
  X   = df[,x_vars]
  Z   = df[,z_vars]  #### ERIC: If no Z is selected, we should have something here, right? 
  
  nhh = nrow(df)
  y   = df[,y_var]
  m   = df[,m_var]
  nvarX = 1 + length(x_vars) # number of X variables selected on the menu. intercept is added here
  nvarM = 1
  nvarZ = 1 + length(z_vars) # number of z variables selected on the menu. intercept is added here
  
  Data = list(X = cbind(rep(1, nhh), X), 
              y = y, 
              m = matrix(m, ncol=1),
              Z = cbind(rep(1, nhh), Z)
  )
  Prior = list(ma  = c(rep(0, nvarX)),
               Aa  = Aa_var * diag(nvarX),
               mgb = c(rep(0, nvarX+1)), 
               Agb = Abg_var * diag(nvarX+1),
               ml = c(rep(0, nvarZ+1)),
               Al = Al_var * diag(nvarZ+1), #"input Al parameter, default is 0.01"
               nu  = nu_var,
               qy  = c(qy_var, qy_var),
               qm  = c(qm_var, qm_var)
  )
  Mcmc = list(Rep=R_var, keep=keep_var, slambda1=slambda1, slambda2=slambda2)
  
  return(list(Data=Data, Prior=Prior, Mcmc=Mcmc))
}

update_inputs_binary <- function(inputs_binary, seed_var) {
  inputs_binary$Mcmc$seed <- seed_var
  
  return(inputs_binary)
}


# --- Generate list of seeds ----------------------------------------------------
FUN_Mediation_SeedList_ForShiny  = function(seednum_var, seed_var)
{ set.seed(seed_var)
  seed.list <- c(sort(round(runif(seednum_var, 0, 10000))))
  #seed.index <- seq(1, seednum_var, 1)
  #num_seeds <- length(seed.index)  # should be the same as input$seednum_varalize list
  return(seed.list)
}
