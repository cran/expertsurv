## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'center',
  out.width = '70%',
  error = TRUE,
  message = TRUE,
  warning = TRUE
)

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


## ----echo = FALSE, fig.cap = "Expert prior distributions"---------------------
img_temp <- "Vignette_Example_1_Expert_Opinion.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))

## ----echo = FALSE, fig.cap = "Model Comparison"-------------------------------

img_temp <- "Vignette_Example_1_DIC.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))


## ----echo = FALSE, fig.cap = "Survival function with Expert prior"------------
img_temp <- "Vignette_Example_1.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))

## ----echo = FALSE, fig.cap = "Survival function with Expert prior (left) and Vague prior (right)"----

img_temp <- "Vignette_Example_2.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))


## ----echo = FALSE, fig.cap = "Survival difference"----------------------------

img_temp <- "Vignette_Example_3.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))


## ----include = FALSE, eval=FALSE----------------------------------------------
#  library("expertsurv"); library("dplyr"); library("ggplot2")
#  data(bc)
#  set.seed(236236)
#  ## Age at diagnosis is correlated with survival time. A longer survival time
#  ## gives a younger mean age
#  bc$age <- rnorm(dim(bc)[1], mean = 65 - scale(bc$recyrs, scale=F), sd = 5)
#  ## Create age at diagnosis in days - used later for matching to expected rates
#  bc$agedays <- floor(bc$age * 365.25)
#  ## Create some random diagnosis dates between 01/01/1984 and 31/12/1989
#  bc$diag <- as.Date(floor(runif(dim(bc)[1], as.Date("01/01/1984", "%d/%m/%Y"),
#                                 as.Date("31/12/1989", "%d/%m/%Y"))),
#                     origin="1970-01-01")
#  ## Create sex (assume all are female)
#  bc$sex <- factor("female")
#  ## 2-level prognostic variable
#  bc$group2 <- ifelse(bc$group=="Good", "Good", "Medium/Poor")
#  head(bc)
#  
#  
#  ## reshape US lifetable to be a tidy data.frame, and convert rates to per person-year as flexsurv regression is in years
#  survexp.us.df <- as.data.frame.table(survexp.us, responseName = "exprate") %>%
#    mutate(exprate = 365.25 * exprate)
#  survexp.us.df$age <- as.numeric(as.character(survexp.us.df$age))
#  survexp.us.df$year <- as.numeric(as.character(survexp.us.df$year))
#  
#  ## Obtain attained age and attained calendar year in (whole) years
#  bc <- bc %>% mutate(attained.age.yr = floor(age + recyrs),
#                      attained.year = lubridate::year(diag + rectime))
#  
#  ## merge in (left join) expected rates at event time
#  bc <- bc %>% left_join(survexp.us.df, by = c("attained.age.yr"="age",
#                                               "attained.year"="year",
#                                               "sex"="sex"))
#  
#  param_expert_example_gpm <- list()
#  param_expert_example_gpm[[1]] <- data.frame(dist = c("norm"),
#                                           wi = c(1), # Ensure Weights sum to 1
#                                           param1 = c(0.40),
#                                           param2 = c(0.05),
#                                           param3 = c(NA))
#  
#  timepoint_expert <- 20
#  
#  example_gpm  <- expertsurv::fit.models.expert(formula=Surv(recyrs,censrec)~as.factor(group2),data=bc,
#                                 distr=c("gomp"),
#                                 method="mle",
#                                 opinion_type = "survival",
#                                 times_expert = timepoint_expert,
#                                 param_expert = param_expert_example_gpm,
#                                 id_St = min(which(bc$group2 =="Good")))
#  
#  expert_opinion_flex <- example_gpm$misc$flex_expert_opinion[[1]]
#  expert_opinion_flex$bhazard_par <- c(0.5)
#  
#  model.gomp.sep.rs <- expertsurv:::flexsurvreg(Surv(recyrs, censrec)~as.factor(group2),
#                                       data=bc, dist="gompertz",
#                                       anc = list(shape = ~ as.factor(group2)),
#                                       bhazard=exprate,expert_opinion = expert_opinion_flex)
#  
#  
#  
#  ## All-cause survival
#  ss.gomp.sep.rs.surv <- flexsurv::standsurv(model.gomp.sep.rs,
#                                      type = "survival",
#                                      at = list(list(group2 = "Good"),
#                                                list(group2 = "Medium/Poor")),
#                                      t = seq(0,30,length=50),
#                                      rmap=list(sex = sex,
#                                                year = diag,
#                                                age = agedays
#                                      ),
#                                      ratetable = survexp.us,
#                                      scale.ratetable = 365.25,
#                                      newdata = bc)
#  #> Marginal all-cause survival will be calculated
#  #> Calculating marginal expected survival and hazard
#  
#  # All-cause hazard
#  ss.gomp.sep.rs.haz <- flexsurv::standsurv(model.gomp.sep.rs,
#                                     type = "hazard",
#                                     at = list(list(group2 = "Good"),
#                                               list(group2 = "Medium/Poor")),
#                                     t = seq(0,30,length=50),
#                                     rmap=list(sex = sex,
#                                               year = diag,
#                                               age = agedays
#                                     ),
#                                     ratetable = survexp.us,
#                                     scale.ratetable = 365.25,
#                                     newdata = bc)
#  
#  plot(example_gpm, add.km = T, t = 0:50,plot_opinion  = TRUE)
#  ggsave("inst/image/Vignette_Example_4.png", units = "in", width = 8, height = 6)
#  
#  
#  plot(ss.gomp.sep.rs.surv, expected = T)+theme_bw()+ylab("Survival")+labs(color = "")+
#    scale_color_manual( values = c("Expected" = "Black", "group2=Good" = "red","group2=Medium/Poor" = "blue"),
#      labels = c("Expected" = "Expected", "group2=Good" ="Group - Good" , "group2=Medium/Poor"="Group - Medium/Poor" ))
#  ggsave("inst/image/Vignette_Example_5.png", units = "in", width = 8, height = 6)
#  
#  
#  plot(ss.gomp.sep.rs.haz, expected = T)+theme_bw()+ylab("Hazard")+labs(color = "")+
#    scale_color_manual( values = c("Expected" = "Black", "group2=Good" = "red","group2=Medium/Poor" = "blue"),
#      labels = c("Expected" = "Expected", "group2=Good" ="Group - Good" , "group2=Medium/Poor"="Group - Medium/Poor" ))
#  ggsave("inst/image/Vignette_Example_6.png", units = "in", width = 8, height = 6)
#  
#  

## ----echo = FALSE, fig.cap = "Survival Estimates without General Population Mortality"----

img_temp <- "Vignette_Example_4.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))


## ----echo = FALSE, fig.cap = "All-cause Survival Estimates including General Population Mortality"----
img_temp <- "Vignette_Example_5.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))

## ----echo = FALSE, fig.cap = "All-cause Hazard Functions including General Population Mortality"----
img_temp <- "Vignette_Example_6.png"
knitr::include_graphics(system.file(paste0("image/",img_temp), package = "expertsurv"))

## ----eval=FALSE, include=FALSE------------------------------------------------
#  
#  #Uniform prior on lambda implies non uniform St*
#  
#  t_expert <- 10
#  n <- 1000000
#  a <- 0
#  b <- .3
#  lambda<-runif(n,a,b)
#  St <-exp(-lambda*t_expert)
#  St_vec <- seq(0,1, by = 0.001)
#  
#  St_1 <- exp(-Inf*t_expert)
#  St_2<- exp(-b*t_expert)
#  k <- -log(St_2)/t_expert
#  
#  plot(density(St))
#  #lines(St_vec, -log(St_vec)/t_expert)
#  lines(St_vec, (1/(St_vec*t_expert))/k, col = "red")
#  
#  
#  
#  #Uniform St* implies non uniform lambda
#  
#  t_expert <- 10
#  n <- 1000000
#  a <- 0.01
#  b <- 0.999
#  St<-runif(n,a,b)
#  lambda <- -log(St)/t_expert
#  lambda_vec <- seq(0,1, by = 0.001)
#  
#  lambda_1 <- -log(b)/t_expert
#  lambda_2<- -log(a)/t_expert
#  
#  k <- -log(St_2)/t_expert
#  
#  plot(density(lambda))
#  #lines(St_vec, -log(St_vec)/t_expert)
#  lines(lambda_vec, t_expert*exp(-lambda_vec*t_expert), col = "red")
#  
#  

