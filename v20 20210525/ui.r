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
library(shinycssloaders)
library(DT)

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
  
  dashboardHeader(title = "Bayesian Analysis of Heterogeneous Mediation",titleWidth=450),
  
  dashboardSidebar(
    
    tags$head(
      tags$link(rel="stylesheet", href="https://cdn.jsdelivr.net/npm/katex@0.10.0-beta/dist/katex.min.css", integrity="sha384-9tPv11A+glH/on/wEu99NVwDPwkMQESOocs/ZGXPoIiLE8MU/qkqUcZ3zzL+6DuH", crossorigin="anonymous"),
      tags$script(src="https://cdn.jsdelivr.net/npm/katex@0.10.0-beta/dist/katex.min.js", integrity="sha384-U8Vrjwb8fuHMt6ewaCy8uqeUXv4oitYACKdB0VziCerzt011iQ/0TqlSlv8MReCm", crossorigin="anonymous"),
      tags$script(HTML(ketex_js))
    ),
    
    sidebarMenu(id = "mysidebar",
      menuItem("Input", tabName = "input", icon = icon("th",lib = "font-awesome"),selected = TRUE),
      menuItem("Models", tabName = "input", icon = icon("bar-chart-o",lib = "font-awesome"),
          menuSubItem("Heterogeneous (BM)", tabName = "BM",icon = icon("angle-double-right")),
          menuSubItem("Aggregate", tabName = "aggregate",icon = icon("angle-double-right"))
      ),
      menuItem("Resources", tabName = "ref", icon = icon("book",lib = "font-awesome"))
      )
  ), 
  
  dashboardBody(
    
    tabItems(
      
      # First tab content
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
                                     uiOutput("covariates_z")
                              )
                            ),
                   
                   #column(width = 4,
                   #),
                   #column(width = 8,
                  
                   # NEED TO HAVE THIS CODED
                   #helpText("Your file has " " observations and " "variables."),
                   helpText("Your file has *** observations and *** variables."),  # placeholder
                  
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
               column(6, "Number of MCMC chains for the Binary Mixture model (default is 10)"),
               column(3, uiOutput("select_seednum"))
             ),
             fluidRow(
               column(6, "Initial MH step size for lambda (default is 0.5)"),
               column(3, uiOutput("select_slambda"))
             ),
             br(),
             hr(),
             
             h5(em("Priors:")),
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
      
      # Second tab content
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
                      withSpinner(type = 8, uiOutput("aggregation_results"))
                  ), 
                  box(width=NULL,
                      solidHeader = FALSE,
                      title = "Aggregate Model. Parameter 95% HPDIs",   #table
                      status = "primary",
                      tableOutput("hdpiA_tbl")
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
                      tableOutput("fitA")
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
                  )
           )
        )
     ),
 #
      # Third tab content
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
      textOutput("keepAlive"),
       
       # Results Tab - Binary Mixture - this only appears when 'model_button' is selected
       tabsetPanel(type = 'tabs', id = 'BM_tabs',
                   
           # Heterogenity Parameters
           tabPanel(title="Heterogeneity Parameters (BM)", value='BMpar_het', 
              
              br(),      
              fluidRow(
                  column(width = 4, 
                      box(width=NULL,
                          title = "Results (Estimation may take 5-10 minutes)",
                          status = "primary",
                          solidHeader = TRUE,
                          withSpinner(type = 8, uiOutput("mediation_result"))
                      ),
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Posterior distribution parameter RHO (95% HPD)",  #table
                          status = "primary",
                          tableOutput("hdpiRho_tbl")
                      )
                  ),
                  column(width = 8,
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Histogram of the posterior distribution of parameter rho",  #Figure
                          status = "primary",
                          plotOutput("plotBM_rho", height = "400px")
                      ),
                      box(width=NULL,
                          solidHeader = FALSE,
                          title = "Distribution of the individual posterior means of the probability to mediate",  #Figure
                          status = "primary",
                          plotOutput("plotBM_w", height = "400px")  # TO DO make it a function of number of variables
                      )
                  )
              )      
           ),
           
           # Effects Parameters
           tabPanel(title="Parameters (BM)", value='BMpar',   

              br(),
              fluidRow(
                column(width = 4,
                       box(width = NULL,
                           solidHeader = FALSE,
                           title = "Parameter 95% HPDIs for the Mediating Segment", #table
                           status = "primary",
                           tableOutput("hdpiBM_M_tbl")
                       ),
                       box(width = NULL,
                           solidHeader = FALSE,
                           title = "Parameter 95% HPDIs for the General Segment", #table
                           status = "primary",
                           tableOutput("hdpiBM_S_tbl")
                       )
                ),
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
                      uiOutput("plotsBM_effects.ui")  # TO DO: size to make a function of number of parameters
                  )
                )
              )
           ),
           
           # Fit and diagnostics for BM
           tabPanel(title="Fit (BM)", value='BMfit', 
            
              # button which will re-run model
              #actionButton("runBMRhat", "Run Convergence Diagnostic", style="color: #fff; background-color: #FF0000; border-color: #DC143C; width:30%"),
       
              h5(strong('Model fit and Diagnostics', align='left')),  #table
              tableOutput("test"),
              br(),
              h5(strong('Rhat for each variable (point est.)', align='left')),  #table
              tableOutput("RhatEst"),
              
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
     tabItem(tabName = "ref",
             h3("Resources"),
             "This app is based on the working paper:",
             em("Bayesian Analysis of Heterogeneous Mediation,"),
             "which is currently under review process in the Journal of Consumer Research.",  br(),
             "  Authors' names are temporarily hidden to adhere to the blind review process.",
             br(),
             br(),
             p("This link will be removed when submitting for review:"),
             a(href="https://ssrn.com/abstract=2600140", "link to the paper"),
             br(),
             br(),
             a(href="https://www.rdocumentation.org/packages/coda/versions/0.19-2/topics/gelman.diag", "References for function gelman.diag"),
             br(),
             p("Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511."),
             p("Brooks, S P. and Gelman, A. (1998) General Methods for Monitoring Convergence of Iterative Simulations.
                  Journal of Computational and Graphical Statistics, 7, 434-455."),
             br()
             

     )
  )  # end of tabItems
 
  ) # end of dashboardBody
 
  
 
) # end of dashboardPage
