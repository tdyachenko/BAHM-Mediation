
# The first line is only if this is your first time to load the application.

## Install needed libraries 
install.packages(c("shiny", "bayesm", "HDInterval", "coda", "gtools","dplyr", "readr", "DT", "future", "promises", "ipc", "future.callr",
                   "future.apply", "data.table", "shinyWidgets", "devtools","htmltools", "shinydashboard","shinyjs"))

## Getting the package from github
devtools::install_github("tdyachenko/BAHM-Mediation", ref = "r-pkg",force = TRUE)

# select 1 if you'd like to update the packages (recommended)
# select 3 is no updates are requested

library(bahm)
run_bahm()
