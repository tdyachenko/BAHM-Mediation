#-#
#
# This is the user-interface definition of a Shiny web application.
# tutorial links:
# https://shiny.rstudio.com/tutorial/written-tutorial/lesson2/
#
#-#


# import libraries
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinyjs)
library(DT)

# ??? what does this do?
ketex_js <- " 
$(document).on('shiny:value', function(event) {
  if(event.name.endsWith('_tbl')) {
    var matches = event.value.match(/(%%+[^%]+%%)/g);
    var newvalue = event.value;
    for(var i=0; i<matches.length; i++){
      var code = '\\\\' + matches[i].slice(2,-2);
      newvalue = newvalue.replace(matches[i], katex.renderToString(code));
    }
    event.value = newvalue;
  }
})
"

# define UI
dashboardPage(
  
  dashboardHeader(title = "Bayesian Analysis of Heterogeneous Mediation (BAHM)",titleWidth=550),
  
  # item on the side bar
  dashboardSidebar(
    
    tags$head(
      useShinyjs(),
      tags$link(rel="stylesheet", href="https://cdn.jsdelivr.net/npm/katex@0.10.0-beta/dist/katex.min.css", integrity="sha384-9tPv11A+glH/on/wEu99NVwDPwkMQESOocs/ZGXPoIiLE8MU/qkqUcZ3zzL+6DuH", crossorigin="anonymous"),
      tags$script(src="https://cdn.jsdelivr.net/npm/katex@0.10.0-beta/dist/katex.min.js", integrity="sha384-U8Vrjwb8fuHMt6ewaCy8uqeUXv4oitYACKdB0VziCerzt011iQ/0TqlSlv8MReCm", crossorigin="anonymous"),
      tags$script(HTML(ketex_js))
    ),
    
    sidebarMenu(id = "mysidebar",
      menuItem("How to use", tabName = "how", icon = icon("question",lib = "font-awesome"),selected = TRUE),
      menuItem("Input", tabName = "input", icon = icon("th",lib = "font-awesome")),
      menuItem("Models", tabName = "input", icon = icon("bar-chart-o",lib = "font-awesome"),
          menuSubItem("Heterogeneous (BM)", tabName = "BM",icon = icon("angle-double-right")),
          menuSubItem("Aggregate", tabName = "aggregate",icon = icon("angle-double-right"))
      ),
      menuItem("Resources", tabName = "ref", icon = icon("book",lib = "font-awesome"))
      )
  ), 
  
  dashboardBody(
    
    tabItems(
        
        tabItem(tabName = "how",
                h4(strong("Welcome!")),
        				br(),
                helpText("This application is currently under continuous development. While most functionality works, you may encounter unexpected behavior, please contact <email>,
                         and we will do our best to look into it and correct it if needed. Thank you!"),
                helpText("More information will be provided shortly. Please go to Input Tab to perform the analysis."),
        ),
                
      
      # SIDE tab content - # TAB Input
      tabItem(tabName = "input",
        
        tabsetPanel(type = 'tabs', id = 'input_tabs',
                   
          # Raw Tab - show data that user has inputted
          tabPanel("Data", value='Raw Data', 
                   
                  #fluidRow(
                  #h4('INSTRUCTIONS: Upload your data and select variables for analysis'),
                  hr(),
                            # add widgets for user to select variables
                            # this is dynamically done based on the file they upload using the RenderUI function in the server
                  h4(strong("Upload Data:")),
          				h5(strong("This may take several minutes in this BETA version.")),
          				h5(strong("You should see the table with your data below when it is done.")),
          				br(),
                  fileInput('file1', 
                           em('Upload CSV file or use default data'),
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv",
                                      ".xls")
                  ),
                  h4(strong("Select Variables:")),
                  helpText("Note: only unique variable selections will be allowed"),
                            fluidRow(
                              column(width = 4,
                                     uiOutput("radio_y"),
                                     uiOutput("radio_m")
                              ),
                              column(width = 8,
                                     uiOutput("checkGroup_x"),
                                     conditionalPanel(condition = "input.checkGroup_x.length > 1",
                                        helpText("WARNING: Selecting more than one independent variable will change the interpretation of the results and will assume simultaneous mediation using both variables.")
                                     ),
                                     uiOutput("covariates_z")
                              )
                            ),
                   
                   #column(width = 4,
                   #),
                   #column(width = 8,
                  
                   #show number of observations and variables in the uploaded file
                   textOutput("obs"),

                   box(width = NULL,
                        solidHeader = FALSE,
                        title = "Table: Raw Data",
                        status = "primary",
                        DT::dataTableOutput("Raw_table")
                        )
            #)
                   #)
          ),
            
          # Input Tab - show model settings and default parameters, allow user to change them
          
          tabPanel("Settings for Analysis", value = 'Settings for Analysis', 
                  
             h5("INSTRUCTIONS: Review Settings for Analysis"),
             br(),
             h5(strong('Selected Variables:')), 
             p(textOutput("XSelection")), 
             p(textOutput("ySelection")),  
             p(textOutput("mSelection")),
             p(textOutput("ZSelection")),
             br(),
             
             h5(strong('Selected Parameters:')),
             hr(),
             
             h5(em("MCMC:")),
             fluidRow(
               column(6, "R (# of draws)"),
               column(3, uiOutput("select_R"))
             ),
             #uiOutput("select_R"),
             fluidRow(
               column(6, "Random seed/starting value"),
               column(3, uiOutput("select_seed"))
             ),
             #uiOutput("select_seed"),
             fluidRow(
               column(6, "Keep (thinning parameter to minimize autocorrelation in draws (default is R/1000))"),
               column(3, uiOutput("select_keep"))
             ),
             #uiOutput("select_keep"),
             fluidRow(
               column(6, "Burnin (number of saved MCMC draws to remove before analysis (default is 0.5*R/keep))"),
               column(3, uiOutput("select_burnin"))
             ),
             #uiOutput("select_burnin"),
             fluidRow(
               column(6, "Number of MCMC chains for the Binary Mixture model (default is 12)"),
               column(3, uiOutput("select_seednum"))
             ),
             fluidRow(
               column(6, "Initial MH step size for lambda (default is 0.5)"),
               column(3, uiOutput("select_slambda"))
             ),
             fluidRow(
                 column(6, "Confidence Interval Band (default is 0.95)"),
                 column(3, uiOutput("select_ciband"))
             ),
             br(),
             hr(),
             
             h5(em("Priors: (any changes in the settings might impact impact the results of the estimation)")),
             fluidRow(
               column(6, HTML("<i>A<sub>&alpha;</sub> </i> (precision parameter for <i>&alpha;</i> parameters)")),
               column(3, uiOutput("select_Aa"))
             ),
             #uiOutput("select_Aa"),
             fluidRow(
               column(6, HTML("<i>A<sub>&beta;&gamma;</sub> </i> (precision parameter for <i>&beta;</i> and <i>&gamma;</i> parameters)")),
               column(3, uiOutput("select_Abg"))
             ),
             #uiOutput("select_Abg"),
             fluidRow(
               column(6, HTML("<i>A<sub>&lambda;</sub> </i> (precision parameter for <i>&lambda;</i> parameters)")),
               column(3, uiOutput("select_Al"))
             ),
             #uiOutput("select_Al"),
             fluidRow(
               column(6,  HTML("<i>&nu;</i>   (degrees of freedom for the variance parameters)")),
               column(3, uiOutput("select_nu"))
             ),
             #uiOutput("select_nu"),
             fluidRow(
               column(6,  HTML("<i>q<sub>y</sub></i> (scaling parameter for variance <i>&sigma;<sub>y</sub><sup>2</sup></i>, default is <i>var(y))</i>")),
               column(3, uiOutput("select_qy"))
             ),
             #uiOutput("select_qy"),
             fluidRow(
               column(6,  HTML("<i>q<sub>m</sub></i> (scaling parameter for variance <i>&sigma;<sub>m</sub><sup>2</sup></i>, default is <i>var(m))</i>")),
               column(3, uiOutput("select_qm"))
             ),
             #uiOutput("select_g"),
             # this block does not work foe me (TD))
             #if (radioButtons$model_button != 'Aggregate'){
                 #fluidRow(
                #          column(6,  HTML("<i>g</i> (proportion parameter - Beta prior disturbution is symmetric)")),
                #          column(3, uiOutput("select_g"))
                #       ),
                    # },
             br()
          ) # end of Data tabPanel "Settings"
        ) # end of tabsetPanel "Data"
      ),  # end of tabItem "input"
      
      # SIDE tab content - # TAB Models (aggregate)
      tabItem(tabName = "aggregate",
 
       # button which will re-run model
        actionButton("runA", "Run Aggregate Model", style="color: #fff; background-color: #FF0000; border-color: #DC143C; width:30%"),
        br(),
        br(),
        fluidRow(
          column(width = 4, 
                  box(width=NULL,
                      solidHeader = TRUE,
                      title = "Results",
                      status = "primary",
                      uiOutput("aggregation_results")
                  ), 
                  box(width=NULL,
                      solidHeader = FALSE,
                      title = "Aggregate Model. Proportions of Posterior Draws by Quadrant",  #table
                      status = "primary",
                      tableOutput("proportionsA")
                  ),
                  box(width=NULL,
                      solidHeader = FALSE,
                      title = "Model fit",   #table
                      status = "primary",
                      tableOutput("fitA"),
                      tableOutput("fitA_DIC")
                  )
            ),
           column(width = 8, 
                  box(width=NULL,
                      status = "primary",
                      solidHeader = FALSE,
                      title = "Aggregate Model. Plots of the Posterior Draws of Parameters",    #Figure
                      "Numbers in the corners are proportions of the posterior distribution in that quadrant.",
                      br(),br(),
                      #plotOutput("scatterplotsA"),
                      uiOutput("plotA.ui"),
                      h5('Downloads:', align='left'),
                      downloadButton('downloadPDF', 'PDF of Plots' ),
                      br()
                  ),
                  box(width=NULL,
                      solidHeader = FALSE,
                      title = "Aggregate Model. Parameter 95% HPDIs",   #table
                      status = "primary",
                      tableOutput("hdpiA_tbl")
                  )
           )
        )
     ),
 #
      # SIDE tab content - # TAB Models (BM)
      tabItem(tabName = "BM",         
        
       # button which will re-run model
       actionButton("runBM", "Run Heterogeneous (BM) Model", style="color: #fff; background-color: #FF0000; border-color: #DC143C; width:30%"),
       br(),
       br(),
       
       tags$head(
        HTML(
          "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
        )
      ),
      # textOutput("keepAlive"),
       
       # Results Tab - Binary Mixture - this only appears when 'model_button' is selected
       tabsetPanel(type = 'tabs', id = 'BM_tabs',
                   
        # Summary
        tabPanel(title="Summary: Binary Mixture Model estimation", value='BMpar_het', 
                
                br(),      
                fluidRow(
                  column(width = 8, 
                         box(width=NULL,
                             title = "Results (Estimation may take 5-10 minutes)",
                             status = "primary",
                             solidHeader = TRUE,
                             uiOutput("mediation_result")
                         )
                  ),
                )
        ),
                         
           # Parameters (Step 1)
           tabPanel(title="Step 1. Parameters (BM)", value='BMpar_het', 
              
              br(),      
              fluidRow(
                  column(width = 8,
                          box(width = NULL,solidHeader = FALSE,
                              title = "Proportions of posterior distribution by Quadrant",  #table
                              status = "primary",
                              DT::dataTableOutput("proportionsBM")
                          ),
                          box(width = NULL,solidHeader = FALSE,
                              title = "Plots of the Posterior Draws of Parameters",  #Figure
                              "Numbers in the corners are proportions of the posterior distribution in that quadrant",
                              status = "primary",
                              uiOutput("plotsBM_effects.ui"),
                              h5('Downloads:', align='left'),
                              downloadButton('downloadPDF_b', 'PDF of Plots' ),
                          )
                         ),
                  column(width = 8,
                         box(width = NULL,
                             solidHeader = FALSE,
                             title = "Parameter 95% HPDIs for the Mediating Segment (M)", #table
                             status = "primary",
                             tableOutput("hdpiBM_M_tbl")
                         ),
                         box(width = NULL,
                             solidHeader = FALSE,
                             title = "Parameter 95% HPDIs for the General Segment (S)", #table  # TODO change to M* if segmentFlag ==2
                             status = "primary",
                             tableOutput("hdpiBM_S_tbl")
                         )
                  )
                )
              ),
                  
          # Heterogeneity Parameters (Steps 2a and 2b)
          tabPanel(title="Steps 2a and 2b. Heterogeneity in Mediation", value='BMpar',   
                 
            br(),
            fluidRow(    
                column(width = 4,
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Distribution of the individual probabilities to mediate",  #Figure
                          status = "primary",
                          plotOutput("plotBM_w", height = "400px")  # LATER make it a function of number of variables (Maybe)
                      ),
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Histogram of the posterior distribution of parameter rho",  #Figure
                          status = "primary",
                          plotOutput("plotBM_rho", height = "400px")
                      ),
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Histogram of the posterior distributions(s) of parameter(s) Lambda",  #Figure (THIS MIGHT BE MORE THAN ONE!)
                          status = "primary",
                          plotOutput("plotBM_lambda", height = "400px")
                      )
                ),
                column(width = 4,
                      h5('Downloads:', align='left'),
                      downloadButton('downloadPDF_BM', 'PDF of Plots' ),
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Posterior distribution of parameter RHO (95% HPDI)",  #table
                          status = "primary",
                          tableOutput("hdpiRho_tbl")
                      ),
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Posterior distribution(s) of parameter(s) Lambda (95% HPDI)",  #table with multiple rows - # of covariates in Z + intercept
                          status = "primary",
                          tableOutput("hdpiLambda_tbl")
                      )
                  ),
              )      
           ),
        
           # Fit and diagnostics for BM
           tabPanel(title="Fit (BM)", value='BMfit', 
            
              # button which will re-run model
              #actionButton("runBMRhat", "Run Convergence Diagnostic", style="color: #fff; background-color: #FF0000; border-color: #DC143C; width:30%"),
       
              h5(strong('Model fit and Diagnostics', align='left')),  #table
              tableOutput("fitBM"),
              br(),
              h5(strong('Rhat for each variable (point est.)', align='left')),  #table
              tableOutput("RhatEst"),
              
              br(),
              h5(strong('Lambda Rejection rate (shoud be around 0.65-0.75)', align='left')),  #table
              tableOutput("rejectRate"),
              
              #h4('Table: Convergence Diagnostics', align='left'),
              #textOutput("BestSeed"),
              #tableOutput("Rhat_sol1"),
              #tableOutput("Rhat_sol2"),
              br()
           ),
           
           # MCMC 
           tabPanel(title="MCMC (BM)", value='BMmcmc', 
                    
              downloadButton(outputId = 'downloadMCMC_BM', label = 'Download'),
             
              box(width=12,
                  solidHeader = TRUE,
                  title = "MCMC traceplots",  #Figure
                  status = "primary",
                  plotOutput("plot_MCMC_BM", height = "800px")  
                )
              
           )
       )  # end of tabsetPanel for BM
     ), # end of tabItem "BM"
 #                  
     # SIDE tab content - # TAB Models (aggregate) 
     tabItem(tabName = "ref",
             h4(strong("Resources")),
     				 br(),
             "This app is based on the paper ",
             em("Is Your Sample Truly Mediating? Bayesian Analysis of Heterogeneous Mediation (BAHM)"),
                "by Tatiana L. Dyachenko (University of Georgia, Terry College of Business) and Greg M. Allenby
                (The Ohio State University, Fisher College of Business). The paper is forthcoming in the Journal of Consumer Research.",  
     				 br(),
     				 br(),
             " The paper with a short tutorial will be avaiable shortly using the link below.",
             br(),
             br(),
             #p("This link will be removed when submitting for review:"),
             a(href="https://doi.org/10.1093/jcr/ucac041", "link to the paper"),
             br(),
             br(),
             a(href="https://www.rdocumentation.org/packages/coda/versions/0.19-2/topics/gelman.diag", "References for function gelman.diag"),
             br(),
             p("Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511."),
             p("Brooks, S P. and Gelman, A. (1998) General Methods for Monitoring Convergence of Iterative Simulations.
                  Journal of Computational and Graphical Statistics, 7, 434-455."),
             br(),
<<<<<<< HEAD
     				 "Development assistance for this application was provided by ",
=======
     				 "Development assistance for this was provided by ",
>>>>>>> 3089b22639f81e7a88a8da8b6337f08694766c56
     				 a(href="https://omnianalytics.org", "Omni Analytics Group"),
     				 br()

     )
  )  # end of tabItems
 
  ) # end of dashboardBody
 
  
 
) # end of dashboardPage
