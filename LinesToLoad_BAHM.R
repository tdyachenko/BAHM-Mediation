
# The first line is only if this is your first time to load the application.

# In RStudio:   Move the cursor to the command line and click "Run" bottom 
# ??? In R console: copy and paste the line and click Enter

devtools::install_github("tdyachenko/BAHM-Mediation", ref = "r-pkg")

# select 1 if you'd like to update the packages (recommended)
# select 3 is no updates are requested

# Select the lines below and click "Run"

library(bahm)
run_bahm()
