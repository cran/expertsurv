## ---- include = FALSE,echo = FALSE--------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
library("expertsurv")

options(rmarkdown.html_vignette.check_title = FALSE)



## -----------------------------------------------------------------------------
#A param_expert object; which is a list of 
#length equal to the number of timepoints
param_expert_example1 <- list()

#If we have 1 timepoint and 2 experts
#dist is the names of the distributions
#wi is the weight assigned to each expert (usually 1)
#param1, param2, param3 are the parameters of the distribution
#e.g. for norm, param1 = mean, param2 = sd
#param3 is only used for the t-distribution and is the degress of freedom.
#We allow the following distributions:
#c("normal","t","gamma","lognormal","beta") 


param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
                                         wi = c(0.5,0.5), # Ensure Weights sum to 1
                                         param1 = c(0.1,0.12),
                                         param2 = c(0.005,0.005),
                                         param3 = c(NA,3))
param_expert_example1

#Naturally we will specify the timepoint for which these probabilities where elicited

timepoint_expert <- 14


#In case we wanted a second timepoint -- Just for illustration

# param_expert_example1[[2]] <- data.frame(dist = c("norm","norm"),
#                                          wi = c(1,1),
#                                          param1 = c(0.05,0.045),
#                                          param2 = c(0.005,0.005),
#                                          param3 = c(NA,NA))
# 
# timepoint_expert <- c(timepoint_expert,18)


## ----echo = FALSE, fig.cap = "Expert prior distributions", out.width="75%"----
knitr::include_graphics("plots/Vignette_Example_1_Expert_Opinion.png")


## ----echo = FALSE, fig.cap = "Model Comparison", out.width="75%"--------------
knitr::include_graphics("plots/Vignette_Example_1_DIC.png")

## ----echo = FALSE, fig.cap = "Survival function with Expert prior", out.width="75%"----
knitr::include_graphics("plots/Vignette_Example_1.png")

## -----------------------------------------------------------------------------
#Check the coding of the arm variable
#Comparator is 0, which is our id_St
unique(data$arm)


## ----echo = FALSE, fig.cap = "Survival function with Expert prior (left) and Vague prior (right)", out.width="75%"----
knitr::include_graphics("plots/Vignette_Example_2.png")

## ----echo = FALSE, fig.cap = "Survival difference", out.width="75%"-----------
knitr::include_graphics("plots/Vignette_Example_3.png")

