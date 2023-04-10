# Bayesian Analysis of Heterogeneous Mediation (BAHM)
 
This project hosts code for Bayesian heterogneous meditaion described in
 
_Tatiana L Dyachenko, Greg M Allenby (2022 in press), Is Your Sample Truly Mediating? Bayesian Analysis of Heterogeneous Mediation (BAHM), 
Journal of Consumer Research, 	https://doi.org/10.1093/jcr/ucac041_

(new code for categorical outcome (y) and mediator (m) variables is coming soon)
  
    
## You can get an R package "bahm" following instructions below. This will open the application on your device.

1. You need to have R installed on your computer:
	- go to https://cran.r-project.org/mirrors.html
	- select a mirror/location and download and install R software
	- you can also install R and RStudio (which provide some convenience) from one place:
  https://posit.co/download/rstudio-desktop/

2. Open file called LinesToLoad_BAHM.R and follow the instructions in the file (3 lines)

```
## Getting the package from github 
devtools::install_github("tdyachenko/BAHM-Mediation", ref = "r-pkg")

## Loading the package 
library(bahm)

## Starting the app on your device
run_bahm()
```
Alternatively, the app can be run as a web application, but it might be much slower.
https://bayesianmediationanalysis.shinyapps.io/BAHM/

