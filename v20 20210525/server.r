#-#
# 
# This is the server logic for a Shiny web application.
#
#-#

#setwd("D:/Dropbox/Mediation/Shiny APP dev/R code/Shiny Code/v5")
#setwd("/Users/kgedney/Dropbox/Mediation/Shiny APP dev/R code/Shiny Code/v4")

# import libraries
library('shiny')
library('bayesm')
library('HDInterval')
library('coda')
library('gtools')
library(dplyr)
library(readr)
library(DT)
library(future)
library(promises)
library(ipc)
library(future.callr)
library(future.apply)

plan(list(
  tweak(callr, workers = max(1, availableCores() %/% 4)),
  tweak(multisession, workers = max(1, availableCores() %/% 2))
))

# load sample data to use as default # moved to see if this helps on the server
sample_df <- read.csv("sample_data_Loyalty.csv")

source('inputs_helpers.r')
source('output_helpers.r')
source('agg_helpers.r')
source('BM_helpers.r')
source('BM_Rhat_helpers.r')
source("FUN_Mediation_LCRM_2class_MS_Gibbs_Moderated_forShinyApp.R")



shinyServer(function(input, output, session) {
  
  # load helper functions
  
  
  
  #------------------------- get data --------------------------------#
  # 1. define df - the reactive function is explained here: 
  # https://shiny.rstudio.com/tutorial/written-tutorial/lesson6/
  df <- reactive({
    
    # check for input
    inFile <- input$file1
    
    # if no file is uploaded by user, then return the sample file
    if (is.null(inFile)){
      return(sample_df)
    }
    
    # otherwise, use the file uploaded by user
    else {
      validate(need(sum(is.na(read.csv(inFile$datapath))) == 0, "Error: please remove missing values in your file and upload again"))
      return(read.csv(inFile$datapath))
    }
  })
  
  output$obs <- renderText({
      return(paste0("Your file has ", nrow(df()), " observations and ", ncol(df()), " variables."))
  })
  
  df_column_list <- reactive({
    return(names(df()))
  })
  
  #---- add renderUI functionality to create dropdowns based on file uploaded -----#
  
  # Data
  output$radio_y <- renderUI({
    selectInput("radio_y", label = h5(em("Dependent Variable (y). Select one.")),
                choices = df_column_list(), selected = "y1")
  })
  
  output$radio_m <- renderUI({
    selectInput("radio_m", label = h5(em("Mediator (m). Select one.")),
                choices = df_column_list()[df_column_list() != input$radio_y], 
                selected = 'avg_m')
  })
  
  output$checkGroup_x <- renderUI({
    input$covariates_z  # ????
    df()
    isolate({
      selectizeInput("checkGroup_x",
                     label    = h5(em("Independent Variables (X). Select all that apply. You must select at least one.
                       For the aggregate model, generate the interactions variables BEFORE uploding the file. 
                       Then, select those variables as X variables.")), 
                     choices  = setdiff(df_column_list(), c(input$radio_y, input$radio_m, input$covariates_z)),
                     selected = if (!is.null(input$checkGroup_x)[1]) input$checkGroup_x else NA,
                     options = list(
                       placeholder = "Choose"
                     ),
                     multiple = TRUE)
    })
  })
  
  # NEED to have some warning pop up to not proceed unless y,m, and at least one x is selected
  
  output$covariates_z <- renderUI({
    input$checkGroup_x # ????
    df()
    isolate({
      tagList(
      selectizeInput("covariates_z",
                     label    = h5(em("Covariates. Leave BLANK if no covariates. 
                                      ")), 
                     choices  = setdiff(df_column_list(), c(input$radio_y, input$radio_m, input$checkGroup_x)),
                     selected = if (!is.null(input$covariates_z)[1]) input$covariates_z else NA,
                     options = list(
                       placeholder = "None"
                     ),
                     multiple = TRUE),
        
      # checkboxInput("interactions", "Would you like to create (an) interaction term(s) for the aggregate model?
       #             If the interactions for the aggregate model have been generated, select them as independent variables above and do not check this box.")
      )
    })
  })
  
  
  # Insert a check mark.
  # "Would you like to create an interaction term for the aggregate model?"
  # for the AGGREGATE model only:
  # Z variables become additional X variables regardless of the response
  # If checked interaction box, in the input file, generate new variables for each X and Z and add to other X variable
  
  # Priors
  output$select_Aa <- renderUI({
    numericInput("select_Aa", label = NULL, value = 0.01, step = 0.1)
  })
  output$select_Abg <- renderUI({
    numericInput("select_Abg", label = NULL, value = 0.01, step = 0.1)
  })
  output$select_Al <- renderUI({
    numericInput("select_Al", label = NULL, value = 0.01, step = 0.1)
  })
  output$select_nu <- renderUI({
    numericInput("select_nu", label = NULL, value = 5, step=1)
  })
  output$select_qy <- renderUI({
    if (is.null(input$radio_y) || !(input$radio_y %in% names(df()))) return(NULL)
    
    numericInput("select_qy", label = NULL, value = round(var(df()[,input$radio_y]), 3), step = 0.5)
  })
  output$select_qm <- renderUI({
    if (is.null(input$radio_m) || !(input$radio_m %in% names(df()))) return(NULL)
    
    numericInput("select_qm", label = NULL, value = round(var(df()[,input$radio_m]), 3), step = 0.5)
  })
  
  # MCMC
  output$select_R <- renderUI({
    #numericInput("select_R", label=("R (# of draws)"), value = 10000, step=100)
    numericInput("select_R", label=NULL, value = 10000, step=100)
  })
  output$select_seed <- renderUI({
    #numericInput("select_seed", label=('Random seed/starting value'), value = 123, step=1)
    numericInput("select_seed", label=NULL, #value = round(runif(1,0,10^5)))
                 value = 125                   )  # FOR TESTING ONLY
  })
  
  output$select_keep <- renderUI({
    #numericInput("select_keep", label=('Keep (thinning parameter to minimize autocorrelation in draws (default is R/1000))'), 
    numericInput("select_keep", label=NULL, value = 10, step=10)
  })
  
  output$select_burnin <- renderUI({
    #numericInput("select_burnin", label=("Burnin (number of saved MCMC draws to remove before analysis (default is 0.5*R/keep))"), 
    numericInput("select_burnin", label=NULL, 
                 value = round(input$select_R/20,0), step=1)
  })
  output$select_seednum <- renderUI({
    numericInput("select_seednum", label=NULL,   value = 10, step=1)
  })
  output$select_slambda <- renderUI({
    numericInput("select_slambda", label = NULL, value = 0.5, step=0.01)
  })
  
  output$ySelection <- renderText({
    paste0('Dependent Variable: ', input$radio_y)
  })
  
  output$mSelection <- renderText({
    paste0('Mediator: ', input$radio_m)
  })
  
  output$XSelection <- renderText({
    paste0('Independent Variables: ', paste(input$checkGroup_x, collapse = ", "))
  })
  
  output$ZSelection <- renderText({
    paste0('Covariates: ', paste(input$covariates, collapse = ", "))
  })
  
  # Forces output elements to initialize without tab being click
  outputOptions(output, 'radio_y', suspendWhenHidden=FALSE)
  outputOptions(output, 'radio_m', suspendWhenHidden=FALSE)
  outputOptions(output, 'checkGroup_x', suspendWhenHidden=FALSE)
  outputOptions(output, 'covariates_z', suspendWhenHidden=FALSE)
  outputOptions(output, 'select_Aa', suspendWhenHidden=FALSE)
  outputOptions(output, 'select_Abg', suspendWhenHidden=FALSE)
  outputOptions(output, 'select_Al', suspendWhenHidden=FALSE)
  outputOptions(output, 'select_nu', suspendWhenHidden=FALSE)
  outputOptions(output, 'select_qy', suspendWhenHidden=FALSE)
  outputOptions(output, 'select_qm', suspendWhenHidden=FALSE)
  outputOptions(output, "select_R", suspendWhenHidden=FALSE)
  outputOptions(output, "select_seed", suspendWhenHidden=FALSE)
  outputOptions(output, "select_keep", suspendWhenHidden=FALSE)
  outputOptions(output, "select_burnin", suspendWhenHidden=FALSE)
  outputOptions(output, "select_seednum", suspendWhenHidden=FALSE)
  outputOptions(output, "select_slambda", suspendWhenHidden=FALSE)
  
  
  #------------------ show Raw data ----------------------------------#
  # 2. render Raw Data inFile_table "Input" table
  output$Raw_table <- DT::renderDataTable({
    df()
  }, options = list(scrollX = TRUE))
  
  output$keepAlive <- renderText({
    req(input$count)
    paste("keep alive ", input$count)
  })
  
  #------------------ User Selected Inputs ----------------------------------#
  # 3. Define Model Inputs and Display in Inputs tab
  
  # calculate inputs only when Run Model button is clicked.
  # ref: https://shiny.rstudio.com/articles/action-buttons.html
  input_listA <- reactive({
    
    #req(input$select_R)
    #req(input$select_seed)
    #req(input$select_keep)
    # validate(need(input$select_R, message="cannot run! Visit 'Inputs for Analysis' tab"))
    
    # Data
    x_vars <- input$checkGroup_x
    y_var  <- input$radio_y
    m_var  <- input$radio_m
    z_var  <- input$covariates_z   # need for both models now
    
    # Priors
    Aa_var  <- input$select_Aa
    Abg_var <- input$select_Abg
    nu_var  <- input$select_nu
    qy_var  <- input$select_qy
    qm_var  <- input$select_qm
    
    # MCMC
    R_var    <- input$select_R
    seed_var <- input$select_seed
    seednum_var <- input$select_seednum
    keep_var <- input$select_keep
    slambda_var <- input$select_slambda
    
    numx_vars <<- length(input$checkGroup_x)
    
    return(get_inputs_agg(df(), 
                          x_vars, y_var, m_var, z_var, # Data
                          Aa_var, Abg_var, nu_var, qy_var, qm_var, #Priors
                          R_var, seed_var, keep_var)) # MCMC
    
  })
  
  #-------------------- Aggregate: Run Model ---------------------------------------#
  #  Run Aggregate Model! 
  
  # add progress bar https://shiny.rstudio.com/articles/progress.html
  aggregate_outputs <- reactiveValues(checkGroup_x = NULL, select_burnin = NULL, output_listA = NULL)
  
  parallel_compute_A <- function(model_data, model_prior, model_mcmc, progress) {
    progress$inc(amount = 1)
    progress$set(message = paste0("Beginning model fitting."))
    
    model_run <- FUN_Mediation_LCRM_2class_MS_Gibbs_Moderated_forShinyApp(
      Model = 1,
      Data  = model_data,
      Prior = model_prior,
      Mcmc  = model_mcmc)
    
    progress$inc(amount = 2)
    progress$set(message = paste0("Beginning model post-processing..."))
    
    return(model_run)
  }
  
  observeEvent(input$runA, {
    aggregate_outputs$checkGroup_x <<- input$checkGroup_x
    aggregate_outputs$radio_m <<- input$radio_m
    aggregate_outputs$select_burnin <<- input$select_burnin
    
    model_data <- input_listA()$Data
    model_prior <- input_listA()$Prior
    model_mcmc <- input_listA()$Mcmc
    
    progress <- AsyncProgress$new(session, min = 1, max = 4, message = "Beginning Model Estimation", detail = "This may take several minutes...")
    
    my_future <- future_promise({
      parallel_compute_A(model_data, model_prior, model_mcmc, progress)
    }, seed=TRUE)
    
    then(
      my_future,
      onFulfilled = function(value) {
        progress$close()
        
        aggregate_outputs$output_listA <<- value
      },
      onRejected = NULL
    )
    
    return(NULL)
  })
  
  agg_prop <- reactive({
    if(is.null(aggregate_outputs$output_listA)) {return(NULL)}
    return(FUN_PDF_Mediation_AlphaBetaProportion_forShiny(model = 1,
                                                          filename = aggregate_outputs$output_listA,
                                                          x_vars   = aggregate_outputs$checkGroup_x,
                                                          burnin   = aggregate_outputs$select_burnin))
  })
  
  output$aggregation_results <- renderUI({
    if (is.null(aggregate_outputs$output_listA)) {
      return(HTML("Click Run Aggregate Model to Begin!"))
    }
    
    prop_agg <- agg_prop()
    maxloc <- which.max(prop_agg)
    colnum <- maxloc %/% 4 + 1
    max_quad <- rownames(prop_agg)[maxloc %% 4]
    
    if (any(prop_agg >= .95)) {
      ## Mediation
      tagList(
        strong("According to the aggregate model, mediation is present in the sample for variable(s)", paste0(aggregate_outputs$checkGroup_x[colnum]),". "),
        paste0(scales::percent(max(prop_agg)), " of the joint posterior distirbution of parameters is in quadrant ", max_quad, " for variable(s) ",
               paste0(aggregate_outputs$checkGroup_x[colnum]),". "),
        br(),br(),
        "The independent variables are: ", paste0(aggregate_outputs$checkGroup_x, collapse = ", "), ".",
        br(),br(),
        "You can download the chart on the right as PDF by clicking on the button below."
      )
    } else {
      tagList(
        strong("According to the aggregate model, mediation is not present in the sample as proposed for variable(s)."),
        paste0("This is based on the analysis of the joint posterior of parameters
                          where the largest propotion of the distribution is ", scales::percent(max(prop_agg)), ".")
      )
    }
  })
  
  # calculate proportion of posterior draws in each quadrant
  output_proportions <- reactive({
    return(agg_prop())
  })
  
  
  # calculate HPD intervals 95%
  output_HDPI_A <- reactive({ 
    if(is.null(aggregate_outputs$output_listA)) {return(NULL)}
    clean_table(FUN_PDF_Mediation_HDPI_forShiny(model=1,
                                                filename = aggregate_outputs$output_listA,
                                                x_vars   = aggregate_outputs$checkGroup_x,
                                                m_var   = aggregate_outputs$radio_m,
                                                burnin   = aggregate_outputs$select_burnin,
                                                CIband=0.95) # Need to change to a variable/input 7-13-2021
    )
  })
  
  # calculate model fit (LMD NR)
  output_fitLMD <- reactive({ 
    if(is.null(aggregate_outputs$output_listA)) {return(NULL)}
    return(FUN_PDF_Mediation_LMD_NR_Aggregate_forShiny(filename = aggregate_outputs$output_listA,
                                                       burnin   = aggregate_outputs$select_burnin))
  })
  # calculate model fit (DIC)
  output_fitLMD_DIC <- reactive({ 
    if(is.null(aggregate_outputs$output_listA)) {return(NULL)}
      
      mydat <- list(
          y = input_listA()$Data$y,
          X = input_listA()$Data$X,
          m = input_listA()$Data$m
      )
    
    return(FUN_DIC_mediation(mydat,   ### NEED TO SEND THE ORIGINAL DATA that was used for estimation
                             McmcOutput = list(aggregate_outputs$output_listA),
                             burnin   = aggregate_outputs$select_burnin,
                             ModelFlag = 1))
  })
  # calculate model fit (LOO) later
  
  
  plotCountA <- reactive({
    if (is.null(aggregate_outputs$checkGroup_x)) return(2)
    
    length(aggregate_outputs$checkGroup_x)
  })
  plotHeightA <- reactive(350 * max(plotCountA(), 1))
  
  # generate scatterplots of posterior draws for each parameter
  output_scatterplots <- reactive({
    if(is.null(aggregate_outputs$output_listA)) {return(NULL)}
    return(FUN_PDF_Mediation_ScatterPlots_forShiny(model=1,
                                                   dataset  = "",
                                                   filename = aggregate_outputs$output_listA,
                                                   burnin   = aggregate_outputs$select_burnin,
                                                   x_vars   = aggregate_outputs$checkGroup_x))
  })
  
  #--------------- Aggregate: Format Results Page -------------------------------------#  
  
  output$plotA <- renderPlot({
    if (is.null(output_scatterplots())) return(NULL)
    
    output_scatterplots()
  })
  
  #renderPlot
  #renderUI
  # display plots
  output$plotA.ui <- renderUI({  # the height function only works with renderUI, not renderPlot
      plotOutput("plotA",  height = plotHeightA()
    )
  })
  #withProgress(expr = output_scatterplots(), message='Running scatterplots...')
  
  # download button for pdf
  output$downloadPDF <- downloadHandler(filename=function() { 'posterior_draws_A.pdf'} ,
                                        content= function(file){
                                          pdf(file=file)
                                          FUN_PDF_Mediation_ScatterPlots_forShiny(
                                            model=1,
                                            dataset="",
                                            filename=aggregate_outputs$output_listA,
                                            burnin=input$select_burnin,
                                            x_vars = aggregate_outputs$checkGroup_x
                                          )
                                          dev.off()
                                        })
  
  
  # print model fit
  output$fitA <- renderTable(expr = {
      y <- output_fitLMD()
      x <- rbind(y, as.matrix(output_fitLMD_DIC()))
      
      rownames(x) <- c("LMD NR", "DIC")
      
      return(x)
  }, rownames=TRUE, colnames=FALSE, bordered=TRUE)
  
  # print HDPI of posterior draws
  output$hdpiA_tbl <- renderTable(expr=output_HDPI_A(), colnames=TRUE, bordered=TRUE) 
  
  
  ###############
  #  FIX the rownames
  ###############
  
  # print proportions of posterior draws by quadrant
  output$proportionsA <- renderTable(expr = output_proportions(), 
                                     rownames=TRUE, colnames=TRUE, bordered=TRUE)
  
  
  
  #-------------------- Binary: Run Model ---------------------------------------#
  
  # Store outputs for the model as reactive values
  model_inputs <- reactiveValues(inputs = NULL, burnin = 500, checkGroup_x = NULL)
  model_outputs <- reactiveValues(output_listBM = NULL, output_RhatcalcBM = NULL, best.seed = NULL)
  
  my_inputs <- reactive({
      
    model_inputs$inputs <- "SET"
    
    # define inputs outside loop (since these will be the same)
    x_vars <- input$checkGroup_x
    y_var  <- input$radio_y
    m_var  <- input$radio_m
    z_var  <- input$covariates_z
    Aa_var  <- input$select_Aa
    Abg_var <- input$select_Abg
    Al_var <- input$select_Al
    nu_var  <- input$select_nu
    qy_var  <- input$select_qy
    qm_var  <- input$select_qm
    R_var    <- input$select_R
    keep_var <- input$select_keep
    slambda_var <- input$select_slambda
    
    inputs_binary <- get_inputs_binary(df(), x_vars, y_var, m_var, z_var,
                                       Aa_var, Abg_var, Al_var, nu_var, qy_var, qm_var,
                                       R_var, keep_var, slambda_var, slambda_var)
    
    return(inputs_binary)
  })
  
  seed_list <- reactive({
    FUN_Mediation_SeedList_ForShiny(seednum_var = input$select_seednum,
                                    seed_var = input$select_seed)
  })
  
  seed_index <- reactive({
    return(1:length(seed_list()))
  })
  
  model_results <- reactive({
    if (is.null(model_outputs$output_listBM)) return(NULL)
    
    return(model_outputs$output_listBM)
  })
  
  model_mediation <- reactive({
    if (is.null(model_results())) return(NULL)
      
      mydat <- list(
          y = input_listA()$Data$y,
          X = input_listA()$Data$X,
          m = input_listA()$Data$m,
          Z = input_listA()$Data$Z
      )
    
    model_outputs$output_RhatcalcBM <- FUN_Mediation_LMD_RHat_MS_cov(inputdata=mydat,
                                                                     datafile=model_results(),
                                                                     seed.index = seed_index(),
                                                                     burnin   = model_inputs$burnin,
                                                                     RhatCutoff   = 1.25)
    
    model_outputs$best.seed <- as.numeric(model_outputs$output_RhatcalcBM$table_forShiny[1,1])
    
    return(model_outputs$output_RhatcalcBM)
  })
  
  parallel_compute <- function(all_seeds, my_inputs, progress) {
    progress$set(message = "Beginning parallel computation")
     
    in_list <- my_inputs

    x <- future_lapply(1:length(all_seeds), function(j) {
      my_inputs$Mcmc$seed <- all_seeds[j]

      progress$inc(amount = .25)
      progress$set(message = paste0("Seed ", j, " of ", length(all_seeds), " Inputs Created."))
      
      model_run <- FUN_Mediation_LCRM_2class_MS_Gibbs_Moderated_forShinyApp(
        Model = 2,
        Data = my_inputs$Data,
        Prior = my_inputs$Prior,
        Mcmc  = my_inputs$Mcmc)
      
      progress$inc(amount = .25)
      progress$set(message = paste0("Seed ", j, " of ", length(all_seeds), " Model Complete."))
      
      return(model_run)
    }, future.seed=TRUE)
    
    progress$set(message = paste0("Beginning model post-processing..."))
    
    return(x)
  }
  
  observeEvent(input$runBM, {
    all_seeds <- seed_list()
    my_inputs <- my_inputs()
    
    progress <- AsyncProgress$new(session, min = 1, max = length(all_seeds), message = "Beginning Model Estimation", detail = "This may take several minutes...")
    
    my_future <- future_promise({
      parallel_compute(all_seeds, my_inputs, progress)
    }, seed=TRUE)
    
    then(
      my_future,
      onFulfilled = function(value) {
        model_outputs$output_listBM <<- value
        
        progress$close()
      },
      onRejected = NULL
    )
    
    return(NULL)
  })
  
  # download button for csv
  output$download_csv <- downloadHandler(
    filename = function(){"bm_w_means.csv"}, 
    content = function(fname){
      theseed <- as.numeric(model_outputs$output_RhatcalcBM$table_forShiny[1,1])
      
      model_res <- model_outputs$output_listBM[[theseed]]
      ws <- rowMeans(model_res$wdraw[,-1:-model_inputs$burnin])
      
      mytbl <- tibble(
        w = ws
      )
      
      write_csv(mytbl, fname)
    }
  )
  outputOptions(output, "download_csv", suspendWhenHidden=FALSE)
  
  output$mediation_result <- renderUI({
    if (is.null(model_inputs$inputs)) {
      return(HTML("Click Run Heterogeneous (BM) Model to Begin!"))
    } else if (is.null(output_proportionsBM())) {
      return(HTML("Model is now executing, please wait..."))
    }
    
    prop_bm <- output_proportionsBM()[,-1] %>% as_tibble() %>% mutate_all(parse_number) %>% as.matrix()
    
    maxloc <- which.max(prop_bm)
    colnum <- maxloc %/% 4 + 1
    max_quad <- rownames(prop_bm)[maxloc %% 4]
    
    if (any(prop_bm >= .95)) {
      ## Mediation
      
      tagList(
        strong("Mediation is present in the sample."),
        paste0(scales::percent(max(prop_bm)), " of the joint posterior distirbution of parameters α and ß for at least one of the segments is in quadrant ", max_quad, ". 
                          The model estimates that the average probability to mediate in the sample is ", scales::percent(round(output_HDPI()[[3]]$Mean, digits = 4)), 
                          ", which can be also interpreted as a percent of the sample mediating through the proposed mediator."),
        br(),br(),
        "The independent variables are: ", paste0(model_inputs$checkGroup_x, collapse = ", "), ".",
        #br(),br(),
        #"The mean of ρ is ", round(output_HDPI()[[3]]$Mean, digits = 4),
        br(),br(),
        "You can download individual specific proabilities to mediate by clicking the button below.
                          The probabilities are calculated as posterior means of individual parameters w.",
        br(),
        downloadButton("download_csv","Download CSV")
      )
    } else {
      tagList(
        strong("Mediation is not present in the sample as proposed."),
        paste0("This is based on the analysis of the joint posterior of parameters α and ß
                          where the largest propotion of the distribution is ", scales::percent(max(prop_bm)), ".
                          Even if there is a separation of respondents between the segments, none of these segments mediates through the proposed mechanism,
                          (the joint distirbution of α and ß is not away from 0)")
      )
    }
  })
  
  # print Rhat metrics
  output$test  <- renderTable({
    model_outputs$output_RhatcalcBM$table_forShiny  #  Should include DIC now
  }, rownames=TRUE, colnames=TRUE, bordered=TRUE)
  output$RhatEst  <- renderTable({
    model_outputs$output_RhatcalcBM$RhatEst
  }, rownames=TRUE, colnames=TRUE, bordered=TRUE)
  
  
  # select the solution to display
  
  #----------------------------------------------------  
  # calculate HPD intervals 95% for the BEST seed
  
  clean_table <- function(tbl) {
    if (!is.list(tbl)) tbl <- list(tbl)
    
    lapply(tbl, function(x) {
      mystr <- gsub("(alpha|beta|gamma|alpha|rho|lambda)(_\\{[a-zA-Z0-9_]+})?", "%%\\1\\2%%", rownames(x))
      
      y <- as.data.frame(x)
      y[] <- lapply(y, function(column) parse_guess(as.character(column)))
      y$Parameter <- mystr
      
      if ("Variable" %in% names(y)) {
          y <- y %>% select(Variable, Parameter, everything())
          y <- y[,c("Variable", "Parameter", setdiff(names(y), c("Variable", "Parameter")))]
      } else {
          y <- y[,c("Parameter", setdiff(names(y), c("Parameter")))]
      }
      
      return(y)
    })
  }
  
  output_HDPI <- reactive({
    if (is.null(model_mediation())) return(list(NULL, NULL, NULL))
    
    clean_table(FUN_PDF_Mediation_Parameters_MSmixture_forShiny(model_outputs$output_listBM,
                                                                x_vars = input$checkGroup_x,
                                                                m_var = input$radio_m,
                                                                z_var = input$covariates_z,
                                                                seed.list=model_outputs$best.seed, burnin = model_inputs$burnin))
  })
  
  # display ON SCREEN all non-rho HDPIs
  output$hdpiBM_M_tbl <- renderTable(expr = output_HDPI()[[1]], colnames=TRUE, bordered=TRUE, sanitize.text.function = function(x) x)
  output$hdpiBM_S_tbl <- renderTable(expr = output_HDPI()[[2]], colnames=TRUE, bordered=TRUE, sanitize.text.function = function(x) x)
  
  # display ON SCREEN HDPIs for rho and lambda
  output$hdpiRho_tbl <- renderTable(expr = output_HDPI()[[3]], colnames=TRUE, bordered=TRUE, sanitize.text.function = function(x) x)
  output$hdpiLambda_tbl <- renderTable(expr = output_HDPI()[[4]], colnames=TRUE, bordered=TRUE, sanitize.text.function = function(x) x)
  
  #----------------------------------------------------  
  # calculate proportion of posterior draws in each quadrant for selected seeds
  output_proportionsBM <- reactive({
    if (is.null(model_mediation())) return(NULL)
    
    test <- FUN_PDF_Mediation_AlphaBetaProportion_MSmixture_forShiny(model_outputs$output_listBM,
                                                                     seed.list=model_outputs$best.seed,  # hard coded for now, need an input function here later
                                                                     x_vars   = model_inputs$checkGroup_x,
                                                                     burnin   = model_inputs$burnin)
    
    test2 <- as.data.frame(test)
    test2$Var <- rownames(test)
    test2 <- test2[,c("Var", setdiff(colnames(test2), "Var"))]
    test2 <- as.matrix(test2)
    
    return(test2)
  })
  # display ON SCREEN proportions of posterior draws by quadrant
  
  container_dt <- reactive({
    tags$table(
      class = 'table table-bordered',
      tags$thead(
        tags$tr(
          tags$th(''),
          tags$th(class = 'dt-center', colspan = length(model_inputs$checkGroup_x), 'Segment M (mediating)'), ##colspan=length(model_inputs$checkGroup_x)
          tags$th(class = 'dt-center', colspan = length(model_inputs$checkGroup_x), 'Segment G (general)')),
        tags$tr(lapply(colnames(output_proportionsBM()), tags$th))
      ))
  })
  output$proportionsBM <- DT::renderDataTable({
    DT::datatable(output_proportionsBM(), container = container_dt(), rownames = FALSE, colnames = TRUE, class = "",
                  options = list(autoWidth = TRUE,
                                 searching = FALSE,
                                 paging = FALSE,
                                 info = FALSE,
                                 ordering = FALSE,
                                 columnDefs = list(list(className = "dt-center", targets = "_all"),
                                                   list(targets = "_all"))))
  })
  
  #----------------------------------------------------  
  # generate plots of posterior draws for each parameter for selected seeds
  output_scatterplotsBM_effects <- reactive({
    if (is.null(model_mediation())) return(NULL)
    
    FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Effects(dataset  = "",  
                                                                filenamelist = model_outputs$output_listBM,
                                                                seed.list = seed_list(),
                                                                seed.selected = model_outputs$best.seed, # hard coded for now, need an input function here later
                                                                burnin   = model_inputs$burnin,
                                                                x_var   = model_inputs$checkGroup_x)
    
  })    
  output_scatterplotsBM_rho <- reactive({
    if (is.null(model_mediation())) return(NULL)
    
    FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Rho(dataset  = "",  
                                                            filenamelist = model_outputs$output_listBM,
                                                            seed.list = seed_list(),
                                                            seed.selected = model_outputs$best.seed, # hard coded for now, need an input function here later
                                                            burnin   = model_inputs$burnin
                                                            #x_var   = model_inputs$checkGroup_x 
    )                                 
  })
  
  output_scatterplotsBM_lambda <- reactive({
      if (is.null(model_mediation())) return(NULL)
      
      FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_Lambda(dataset  = "",  
                                                                 filenamelist = model_outputs$output_listBM,
                                                                 seed.list = seed_list(),
                                                                 seed.selected = model_outputs$best.seed, # hard coded for now, need an input function here later
                                                                 burnin   = model_inputs$burnin
                                                                 #x_var   = model_inputs$checkGroup_x 
      )                                 
  })
  
  output_scatterplotsBM_w <- reactive({
    if (is.null(model_mediation())) return(NULL)
    
    FUN_PDF_Mediation_ParameterPlots_MSmixture_forShiny_meanW(dataset  = "",  
                                                              filenamelist = model_outputs$output_listBM,
                                                              seed.list = seed_list(),
                                                              seed.selected = model_outputs$best.seed, # hard coded for now, need an input function here later
                                                              burnin   = model_inputs$burnin
                                                              #x_var   = model_inputs$checkGroup_x 
    )
  })    
  # display ON SCREEN scatterplots of posterior draws for each parameter for selected seeds
  output$plotsBM_effects <- renderPlot({
    req(plotCountBM())
    if (plotCountBM() == 0){
      plot.new()
      return()
    }
    expr = output_scatterplotsBM_effects()
  })
  
  plotCountBM <- reactive({
    if (is.null(model_inputs$checkGroup_x)) return(2)
    
    length(model_inputs$checkGroup_x)
  })
  plotHeightBM <- reactive(350 * max(plotCountBM(), 1))
  
  #renderPlot
  #renderUI
  # display plots
  output$plotsBM_effects.ui <- renderUI({  # the height function only works with renderUI, not renderPlot
    plotOutput("plotsBM_effects",  height = plotHeightBM()
    )
  })
  #withProgress(expr = output_scatterplots(), message='Running scatterplots...')
  
  
  #output$plotsBM_effects <- renderPlot({
  #  withProgress(expr = output_scatterplotsBM_effects(), message='Running scatterplots...')
  #})
  output$plotBM_rho <- renderPlot({
    output_scatterplotsBM_rho()
  })
  output$plotBM_lambda <- renderPlot({
      output_scatterplotsBM_lambda()
  })
  output$plotBM_w <- renderPlot({
    output_scatterplotsBM_w()
  })
  
  
  #----------------------------------------------------  
  # MCMC plots
  output_MCMC_BM <- reactive({
    if (is.null(model_mediation())) return(NULL)
    
    FUN_PDF_MCMC_Mediation_forShiny(model = 2, dataset  = "",
                                    filenamelist = model_outputs$output_listBM,
                                    seed.list = seed_list(),
                                    seed.index = seed_index(),
                                    burnin   = model_inputs$burnin)
  })
  # display ON SCREEN scatterplots of posterior draws for each parameter for selected seeds
  output$plot_MCMC_BM <- renderPlot({
    withProgress(expr = output_MCMC_BM(), message='Running MCMC traceplots...')
  })
  
  # download button for pdf
  # This works well only in Browser not in RStudio
  
  output$downloadMCMC_BM <- downloadHandler(
    filename = function() { paste("MCMC for BM model.pdf") },
    content = function(file) {
      pdf(file=file,width=50, height=20)
      FUN_PDF_MCMC_Mediation_forShiny(dataset  = "",
                                      filenamelist = model_results(),
                                      seed.list = seed_list(),
                                      seed.index = seed_index(),
                                      burnin   = model_inputs$burnin)
      dev.off()
    } 
  )
  
  
  
  #--------------- Binary: Format Results Page -------------------------------------#  
  
  # need to calcualte and print on the page
  # LMD_agg = logMargDenNR(output_list()$LLtotal[-1:-burnin])
  
  ###############
  #  FIX the LMD calculation
  ###############
  
  #output$fitLMD <- renderPrint({
  #   logMargDenNR(output_list()$LLtotal[-1:-100])[[1]]
  #})
  
  # print list of model results
  # output$Results <- renderPrint({
  #   output_list()
  # })
  
  
  ###############
  #  FIX the rownames
  ###############
  
  # print proportions of posterior draws by quadrant
  #output$proportionsM <- renderTable(expr = output_proportionsBM()[[1]]$ProportionsM,
  #                                   rownames=TRUE, colnames=TRUE, bordered=TRUE)
  #output$proportionsS <- renderTable(expr = output_proportionsBM()[[2]]$ProportionsS,
  #                                   rownames=TRUE, colnames=TRUE, bordered=TRUE)
  
  
  
  # download button for pdf
   output$downloadPDF_b <- downloadHandler(filename='posterior_draws.pdf',
                                         content= function(file){
                                           pdf(file=file)
                                           FUN_PDF_Mediation_FinalPlots_MSmixture_forShiny_Plot(
                                              dataset="",
                                              filenamelist = model_results(),
                                              seed.list = seed_list(),
                                              seed.selected = c(1,2),
                                              burnin   = model_inputs$burnin
                                              #x_vars   = model_inputs$checkGroup_x
                                            )
                                            dev.off()
                                          })
  
  
  
  
  # Run Binary Model with Parallelization and 4 loops - NOT DONE
  
  #  # setup for loop
  #  num_seeds = 2
  #  seed.list = sort(round(runif(num_seeds, 0, 10000)))
  #  seed.index = seq(1, num_seeds, 1)
  #  
  #  avail_cores <- detectCores()
  #  num_workers <- avail_cores - 1
  #  registerDoParallel(cores=num_workers)
  #  
  #  # get outputBM  (PARALLEL VERSION - NOT DONE)
  #  outputBM <- eventReactive(input$run, {
  #   if(input$model_button == "Heterogeneous (Binary Mixture)"){
  #     
  #     # define inputs outside loop (since these will be the same)
  #     x_vars <- model_inputs$checkGroup_x
  #     y_var  <- input$radio_y
  #     m_var  <- input$radio_m
  #     Aa_var  <- input$select_Aa
  #     Abg_var <- input$select_Abg
  #     nu_var  <- input$select_nu
  #     qy_var  <- input$select_qy
  #     qm_var  <- input$select_qm
  #     R_var    <- input$select_R
  #     keep_var <- input$select_keep
  #     g_var <- input$select_g # only BM
  #     
  #     foreach(i=1:length(seed.list)) %dopar% { 
  #       # different seed each loop
  #       input_list <- get_inputs_binary(df(), x_vars, y_var, m_var,
  #                                       Aa_var, Abg_var, nu_var, qy_var, qm_var, #Priors
  #                                       R_var, seed_var=seed.list[i], keep_var, g_var)
  #       #cat("i=", i, fill=T)
  #       
  #       # run model
  #       return(FUN_Mediation_LCRM_2class_MS_Gibbs_forShinyApp(Data  = input_list()$Data,
  #                                                             Prior = input_list()$Prior,
  #                                                             Mcmc  = input_list()$Mcmc))
  #     }
  #   }
})



