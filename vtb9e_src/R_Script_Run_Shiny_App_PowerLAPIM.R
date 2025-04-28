# Users can download the app and run locally on their computer by executing the following commands in R or Rstudio.

# Check if R packages are installed

list.of.packages = c("nlme","MASS","tidyverse","future.apply","gridExtra","formattable","htmltools",
"shiny","DT","ggplot2","gridExtra","data.table","plyr","dplyr","tidyr","shinyjs")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(nlme)
library(MASS)
library(tidyverse)
library(future.apply)
library(gridExtra)
library(formattable)
library(htmltools)
library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(shinyjs)

library(devtools)
devtools::install_github("ginettelafit/PowerLAPIM", force = T)

library(PowerLAPIM)

# Using Gist: users can launch this app with:
shiny::runGist('1d186b6d9bc76f5e41871ce40e5cee47')
