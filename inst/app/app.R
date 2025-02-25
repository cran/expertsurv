  source("shiny-helper-funs.R")
  
  # Access compile_mods from the global environment
  compile_mods <- .GlobalEnv$compile_mods
  
  options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)
  
 
  ui = fluidPage(shinyjs::useShinyjs(),
                 shinybusy::add_busy_spinner(spin = "semipolar"),			 
				tags$head(tags$title("ShinyExpertsurv")),  
                 titlePanel(
                   div(style = "display: flex; align-items: center;",
                       img(src = "hexsticker.png", height = "150px", style = "margin-left: 20px;"),
                       h1("ShinyExpertsurv", style = "flex-grow: 1;")
                   )
                 ),
                 tags$head(
                   tags$style(HTML(".fixed-height {min-height: 100px; /* Adjust this value as needed */ }")),
                   tags$style(HTML("#run_analysis {
                                      font-size: 18px;
                                      color: #FFFFFF;
                                      background-color: #5bc0de;
                                      border: 2px solid #000000;
                                      width: 100%;
                                      height: 60px;}"))),
                 sidebarPanel(
                   wellPanel(
                     h3("Reference Documents"),
                     tags$a(href = "README.html", "expertsurv Package README"),
                     tags$br(),
                     tags$a(href = "ShinyExpertsurv-Vignette.html", "ShinyExpertsurv App README")),
                   wellPanel(
                     fileInput('df_upload', 'Choose .csv data file to upload',
                               accept = c(".csv")),
                     varSelectInput("variables", "Variable:", data.frame(NULL), multiple = TRUE),
                     p("Data should have the following columns: time and status. If your data has two treatment arms please include an arm column."),
                     numericInput("n_expert", "Number of Experts", value = 1, min = 1),
                     numericInput("n_timepoint", "Number of Timepoints", value = 1,min = 1,  max = 2),
                     numericInput("xlim", "Limit of x-axis on Survival Plot", value = round(10,digits = 0)),
                     checkboxInput(inputId ="expert_opt", label = "Show Advanced options for expert opinion", value = FALSE),
                     checkboxInput(inputId ="MLV_opt", label = "Include Most Likely Values (MLV)", value = FALSE),
                     selectInput(inputId ="pool_type_eval", label = "Pooling approach for experts", 
                                 choices = c("Linear Pool" = "linear pool","Logarithmic Pool"= "log pool"),selected = "linear pool"),
                     selectInput(inputId ="dist_select", label = "Select the best fitting distribution for Expert Pooling", 
                                 choices = c("Best Fitting" = "best",
                                             "Normal"= "normal",
                                             "T-distribution" = "t",
                                             "Gamma" = "gamma",
                                             "Log-Normal" = "lognormal",
                                             "Beta" = "beta"), 
                                 selected = "best"),
                     actionButton(paste0('update_expert'), "Plot/Update Survival Curves and Expert Opinions")),
                   hr(),
                   tabsetPanel(id = "Timepoints",
                               tabPanel("Timepoints1",
                                        numericInput(paste0("time1"), label= "Timepoint", value= 1),
                                        textInput('quant_vec1', 'Enter a Vector of Quantiles', "0.025,0.5,0.975"),
                                        helpText("Enter the judgements in the table below,
                                                              one column per expert. Enter quantile values corresponding to the cumulative probabilities. 
                                                             Enter Most Likely Values (i.e. Mode) for each expert if included."),
                                        shinyMatrix::matrixInput(
                                          inputId = "matrix1",
                                          value = m_default_gen(),
                                          class = "numeric",
                                          cols = list(names = TRUE,
                                                      editableNames = FALSE),
                                          rows = list(names = FALSE,
                                                      editableNames = FALSE)),
                                        shinyMatrix::matrixInput(
                                          inputId = "matrix1_mode",
                                          value = m_default_gen2(),
                                          class = "numeric",
                                          cols = list(names = TRUE,
                                                      editableNames = FALSE),
                                          rows = list(names = TRUE,
                                                      editableNames = FALSE)),
                                        plotOutput(paste0("expert_plot1"))),
                               tabPanel("Timepoints2",
                                        numericInput(paste0("time2"), label= "Timepoint", value= 1),
                                        textInput('quant_vec2', 'Enter a Vector of Quantiles', "0.025,0.5,0.975"),
                                        helpText("Enter the judgements in the table below, one column per expert. Enter quantile values corresponding to the cumulative probabilities."),
                                        shinyMatrix::matrixInput(
                                          inputId = "matrix2",
                                          value = m_default_gen(),
                                          class = "numeric",
                                          cols = list(names = TRUE,
                                                      editableNames = FALSE),
                                          rows = list(names = FALSE,
                                                      editableNames = FALSE)),
                                        shinyMatrix::matrixInput(
                                          inputId = "matrix2_mode",
                                          value = m_default_gen2(),
                                          class = "numeric",
                                          cols = list(names = TRUE,
                                                      editableNames = FALSE),
                                          rows = list(names = TRUE,
                                                      editableNames = FALSE)),
                                        
                                        plotOutput(paste0("expert_plot2"))))),
                 mainPanel(
                   h3("Kaplan-Meier Survival Plot"),
                   plotOutput("plot_km_expert1"),
                   
                   fluidRow(
                     column(width = 3, actionButton("run_analysis", "Run Analysis")),
                     column(width = 3, 
                            selectInput("opinion_type", label = "Choose opinion type",
                                        choices = c("Survival at timepoint(s)" = "survival",
                                                    "Mean difference between survival" = "mean",
                                                    "No expert opinion" = "no_expert"),
                                        selected = "survival")),
                     column(width = 3, class = "fixed-height", 
                            selectInput("stat_type", label = "Choose statistical approach",
                                        choices = c("Frequentist" = "mle", "Bayesian" = "bayes"),
                                        selected = "mle"),
                            conditionalPanel(
                              condition = "input.stat_type == 'bayes'",
                              numericInput("iterations", "Number of iterations:", value = 2000)
                            )),
                     column(width = 3, 
                            shinyWidgets::pickerInput(
                              inputId = "param_mod",
                              label = "Choose models:",
                              choices = c("Exponential" = "exp", "Weibull" = "wei", "Gompertz" = "gomp", 
                                          "Log-Logistic"= "llo", "Log-normal" = "lno", "Generalized-Gamma" = "gga",
                                          "Royston-Parmar (1 knot)" = "rps"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 3"
                              ),
                              multiple = TRUE,
                              selected  = c("exp", "wei")
                            )),
                     column(width = 3, 
                            selectInput("id_trt", label = "Select name of treatment corresponding to expert opinion",
                                        choices =  character(0)))
                   ),
                   
                   fluidRow(
                     column(width = 3, 
                            selectInput("gof_type", label = "Choose goodness of fit measure",
                                        choices = c("AIC" = "aic", "BIC" = "bic"),
                                        selected = "AIC")),
                     column(width = 3, 
                            selectInput("incl_psa", label = "Include Statistical Uncertainty in Plots",
                                        choices = c("Yes" = "yes", "No" = "no"),
                                        selected = "no"))
                   ),
                   
                   plotOutput("plot_gof"),
                   
                   fluidRow(
                     column(width = 3, textInput('file_name', 'Enter filename for saved output', "Output-File")),
                     column(width = 3, 
                            selectInput("outFormat", label = "Report format",
                                        choices = list(html = "html_document", pdf = "pdf_document", Word = "word_document")))
                   ),
                   
                   fluidRow(
                     column(width = 3, downloadButton("save_output", "Download R objects")),
                     column(width = 3, downloadButton("report", "Download report"))
                   )
                 )
                 
  )
  
  
  create_server <- function(compile_mods){function(input, output, session){
    
    #shinyjs::hide("incl_psa") 
    
    value <- reactiveValues(
      m_default = m_default_gen(),
      n_expert_prev = 1,
      quant_vec2 = NULL,
      id_trt = NULL) #Up to max timepoints
    
    
    observeEvent(input$stat_type,{
      if(input$stat_type == "mle"){
        updateSelectInput(session,"gof_type",choices = c("AIC" = "aic", "BIC" = "bic"))
      }else{
        updateSelectInput(session,"gof_type",choices = c("WAIC" = "waic", "PML" = "pml"))
      }
      
    })
    
    
    observeEvent(input$MLV_opt,{
      if(input$MLV_opt){
        shinyjs::show("matrix1_mode")
        shinyjs::show("matrix2_mode")
      }else{
        shinyjs::hide("matrix1_mode")
        shinyjs::hide("matrix2_mode")
        
      }
      
    })
    
    observeEvent(input$df_upload,{
      inFile <- input$df_upload
      df_upload <- utils::read.csv(inFile$datapath)
      value$df_upload <- df_upload
      #varSelectInput("variables", "Variable:", df_upload, multiple = TRUE),
      vars <- names(df_upload)
      # Update select input immediately after clicking on the action button. 
      updateSelectInput(session, "variables","Variable:", choices = vars)
      
    })
    observeEvent(input$expert_opt,{
      #browser()
      if(input$expert_opt){
        shinyjs::show("pool_type_eval")
        shinyjs::show("dist_select")
        shinyjs::show("MLV_opt")
      }else{
        shinyjs::hide("pool_type_eval")
        shinyjs::hide("dist_select")
        shinyjs::hide("MLV_opt")
      }
    })
    
    observeEvent(input$Timepoints, {
      i <- as.numeric(gsub("Timepoints", "", input$Timepoints))
      output$expert_plot <- renderPlot({
        value[[paste0("expert_plot", i)]]
      })
      
    })
    
    
    observeEvent({input$variables},{
    
      if(length(input$variables) == 2 & exists("value$df_work")){
        df_work$arm <- 1
      }
    })
    
    
    observeEvent({input$update_expert},{
      if(length(input$variables) >= 2 ){
        df_work <- dplyr::select(value$df_upload, !!!input$variables)
        
        #browser()    
        if(length(input$variables) == 2){
          df_work$arm <- 1
        }
        
        colnames(df_work) <- c("time", "status", "arm")
        trt_vec <- unique(df_work[["arm"]])
        
        prev_input <- input$opinion_type
        
        
        if(length(trt_vec) == 1){
          result.km <- survfit(Surv(time, status) ~ 1, data = df_work, conf.type="log-log")
          km.data <- data.frame(cbind(result.km[[c("time")]],
                                      result.km[[c("surv")]],
                                      result.km[[c("upper")]],
                                      result.km[[c("lower")]],
                                      arm = 1))
          
          
          
          if(!any(prev_input %in% c("survival","no_expert"))){
            prev_input <- character(0)
          }
          
          updateSelectInput(session,"opinion_type",choices = c("Survival at timepoint(s)" = "survival",
                                                               "No expert opinion" = "no_expert"), selected  = prev_input)
          shinyjs::hide("id_trt") #hide id_trt panel
          value$id_trt <- NULL
          df_work$arm <- 1
        }else{
          
          shinyjs::show("id_trt") #hide id_trt panel
          trt_vec_char <- as.character(trt_vec)
          
          updateSelectInput(inputId = "id_trt", choices = trt_vec_char,selected = head(trt_vec_char,1) )
          
          km.data <- NULL
          for(i in unique(df_work$arm)){
            df_temp <-  filter(df_work, arm == i)
            result.km_temp <- survfit(Surv(time, status) ~ 1, data = df_temp, conf.type="log-log")
            km.data_temp <- data.frame(cbind(result.km_temp[[c("time")]],
                                             result.km_temp[[c("surv")]],
                                             result.km_temp[[c("upper")]],
                                             result.km_temp[[c("lower")]],
                                             arm = i))
            
            km.data <-  rbind(km.data,km.data_temp)
          }
          
          updateSelectInput(session,"opinion_type",
                            choices = c("Survival at timepoint(s)" = "survival",
                                        "Mean difference between survival"= "mean",
                                        "No expert opinion" = "no_expert"),
                            selected = prev_input)
          
        }
        
        colnames(km.data) <- c("Time", "Survival", "upper", "lower", "arm")
        #browser()
        value$km.data <- km.data
        value$df_work <- df_work
        value$id_trt <- input$id_trt
        #Need to adjust for arm
        
        plot_fit <- ggplot2::ggplot(value$km.data, aes(x = Time,y =Survival, col = factor(arm)))+
          geom_step()+
          ylim(0,1)+
          xlim(0, input$xlim)+
          geom_step(aes(x  = Time, y =upper, col = factor(arm)))+
          geom_step(aes(x  = Time, y =lower, col = factor(arm)))+
          theme_light()#+
  
        if(!any(is.na(input[["matrix1"]][,2]))){ # If Expert opinions are not NA values
          
          times_expert_vec <- c()
          df.linear_all <- NULL
          param_expert <- list()
          df.linear <- list()
          scale_vec <- c()
          final_scale <- 0.2 # We want the difference between the final time and density to be less than 0.2
          
          for(i in 1:input$n_timepoint){ #Update
            
            if(input$opinion_type == "survival"){
              St_opinion = 1
            }else{
              St_opinion = 0
            }
            if(input$MLV_opt){
              output_pool <- return_pooled_info(input[[paste0("matrix",i)]], St_indic = St_opinion,dist = input$dist_select, mode =input[[paste0("matrix",i,"_mode")]][1,])
            }else{
              output_pool <- return_pooled_info(input[[paste0("matrix",i)]], St_indic = St_opinion,dist = input$dist_select)
            }
            if(input$opinion_type == "survival"){
              output_pool[[2]] <- output_pool[[2]]+ xlim(c(0,1)) #If survival we want to truncate.
            }
            if(input$n_expert == 1){
              output_pool[[2]][["layers"]][[3]] <-NULL
              output_pool[[2]][["layers"]][[2]] <-NULL
            }  

            value[[paste0("expert_plot",i)]] <- output_pool[[2]]
            times_expert = input[[paste0("time",i)]]
            times_expert_vec <- c(times_expert_vec, times_expert)
            df.curr <-  subset(output_pool[[2]]$data, ftype == input$pool_type_eval)
            df.curr <- rename(df.curr,y = x)
            
            scale_vec <- c(scale_vec,final_scale*(input$xlim-times_expert)/max(df.curr$fx))
            df.linear[[i]] <- df.curr
            output_pool[[1]][,"dist"] <- gsub("normal", "norm", output_pool[[1]][,"dist"])
            param_expert[[i]] <- output_pool[[1]]
          }
          
          for(i in 1:input$n_timepoint){
            df.linear[[i]] <- mutate(df.linear[[i]],x = times_expert_vec[i] + fx*min(scale_vec), 
                     times_expert = times_expert_vec[i])
            df.linear_all <- rbind(df.linear_all, df.linear[[i]])
          }
          
          
          
          value$param_expert <- param_expert
          value$timepoint_expert <- times_expert_vec
          value$df.linear_all <- df.linear_all
          
          if(input$opinion_type == "survival"){
            plot_fit <- plot_fit+
              geom_ribbon(data = df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                          fill = "sky blue", alpha = 0.5, colour = "grey")
          }
          
          
        }
        
        output$plot_km_expert1<- renderPlot(plot_fit)
        
      }  
      
    })
    
    
    
    observeEvent(input$n_timepoint, {
      
      if(input$n_timepoint > 1){
        shiny::showTab(inputId = "Timepoints", target = "Timepoints1")
        shiny::showTab(inputId = "Timepoints", target = "Timepoints2")
      }
      if(input$n_timepoint == 1){
        shiny::showTab(inputId = "Timepoints", target = "Timepoints1")
        shiny::hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
    })
    
    observeEvent(input$opinion_type,{
      shinyjs::hide("plot_gof")
      if(input$opinion_type == "survival"){
        
        updateSelectInput(session,"id_trt", label = "Select treatment ID corresponding to expert opinion")
        
        
        shinyjs::show("time1")
        if(input$n_timepoint > 1){
          shiny::showTab(inputId = "Timepoints", target = "Timepoints1")
          shiny::showTab(inputId = "Timepoints", target = "Timepoints2")
        }
        if(input$n_timepoint == 1){
          shiny::showTab(inputId = "Timepoints", target = "Timepoints1")
          shiny::hideTab(inputId = "Timepoints", target = "Timepoints2")
        }
        
      }
      
      if(input$opinion_type == "mean"){
        
        updateSelectInput(session,"id_trt", label = "Select treatment ID corresponding to expert opinion - Mean difference of selected treatment vs other treatment")
        shiny::showTab(inputId = "Timepoints", target = "Timepoints1")
        shinyjs::hide("time1")
        shiny::hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
      
      if(input$opinion_type == "no_expert"){
        shiny::hideTab(inputId = "Timepoints", target = "Timepoints1")
        shiny::hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
      
      if(input$opinion_type == "survival" | input$opinion_type == "mean"){
        shinyjs::show("n_expert")
        if(input$opinion_type == "survival" ){
          shinyjs::show("n_timepoint") 
        }else{
          shinyjs::hide("n_timepoint")
        }
        
        shinyjs::show("expert_opt")
        #browser()
        if(input$expert_opt){
          shinyjs::show("pool_type_eval")
          shinyjs::show("dist_select")
          shinyjs::show("MLV_opt")
        }else{
          shinyjs::hide("pool_type_eval")
          shinyjs::hide("dist_select")
          shinyjs::hide("MLV_opt")
        }
        
        if(is.null(value$id_trt)){
          
          shinyjs::hide("id_trt")
        }else{
          shinyjs::show("id_trt")
        }
        
        
      }else{ #No Expert Opinion
        shinyjs::hide("n_expert")
        shinyjs::hide("n_timepoint")
        shinyjs::hide("expert_opt")
        shinyjs::hide("MLV_opt")
        shinyjs::hide("pool_type_eval")
        shinyjs::hide("dist_select")
        shinyjs::hide("id_trt")
        
      }
      
    })
    
    observeEvent(input$id_trt,{
      value$id_trt <- input$id_trt
    })
    
    observeEvent(input$n_expert, {
      # browser()
      if(input$n_expert == 1){
        #shinyjs::hideElement(id = "pool_type_eval")
        shinyjs::hide("pool_type_eval")
        
      }else{
        #shinyjs::showElement(id = "pool_type_eval")
        shinyjs::show("pool_type_eval")
      }
      
      for(i in 1:2){ #Modify this force it to me 2 which is the max number of timepoints
        mat_exist <- input[[paste0("matrix",i)]]
        #browser()
        mat_exist_mode <- input[[paste0("matrix",i,"_mode")]]
        
        if(input$n_expert > value$n_expert_prev){
          extra_cols <- input$n_expert - value$n_expert_prev 
          mat_bind <- matrix(nrow = nrow(mat_exist), ncol = extra_cols)
          mat_exist <- cbind(mat_exist,mat_bind)
          
          mat_bind_mode <- matrix(nrow = nrow(mat_exist_mode), ncol = extra_cols)
          mat_exist_mode <- cbind(mat_exist_mode,mat_bind_mode)
          
        }else if(input$n_expert == value$n_expert_prev){
        } else{
          mat_exist <- mat_exist[,1:(input$n_expert+1),drop = F]
          mat_exist_mode <- mat_exist_mode[,1:(input$n_expert),drop = F]
        }
        colnames(mat_exist) <- c("Cum Prob", paste0("Expert_", 1:input$n_expert))
        shinyMatrix::updateMatrixInput(session, paste0("matrix",i), value=mat_exist )
        
        colnames(mat_exist_mode) <-  paste0("Expert_", 1:input$n_expert)
        shinyMatrix::updateMatrixInput(session, paste0("matrix",i,"_mode"), value=mat_exist_mode )
      }
      value$n_expert_prev <- input$n_expert
      
    })
    
    
    toListen <- reactive({
      list(input$quant_vec1,input$quant_vec2)
    })
    
    observeEvent(toListen(),{
      
      for(i in 1:input$n_timepoint){#max number of quant_vec
        #browser()
        if(!is.null(input[[paste0("quant_vec",i)]])){
          quant_vec_temp <- input[[paste0("quant_vec",i)]]
          quant_num <- as.numeric(unlist(strsplit(quant_vec_temp,",")))
          if(length(quant_num)==0){
            new_mat <-   matrix(ncol = input$n_expert +1, nrow = 1) # Handle case when nothing is entered
            colnames(new_mat) <- c("Cum Prob", paste0("Expert_",1:input$n_expert ))
            
          }else{
            mat_exist <- input[[paste0("matrix",i)]]
            retain_quant_index <-which(mat_exist[,1] %in% quant_num)
            retain_quant <- mat_exist[retain_quant_index,1]
            change_quant_index <- which(quant_num %!in% retain_quant)
            new_mat <- matrix(ncol = input$n_expert +1, nrow = length(quant_num))
            colnames(new_mat) <- c("Cum Prob", paste0("Expert_",1:input$n_expert ))
            
            if(length(retain_quant_index)>0){
              new_mat[1:length(retain_quant_index),] <-mat_exist[retain_quant_index,]
              
            }
            if(length(change_quant_index)>0){
              new_mat[(length(retain_quant_index)+1):nrow(new_mat),1] <- quant_num[change_quant_index]
            }
          }
          shinyMatrix::updateMatrixInput(session, paste0("matrix",i), value=new_mat)
          
        }
      }})
    
    
    
    observeEvent(input$run_analysis, {
      
      if(input$id_trt == ""){
        id_trt_work<- 1
      }else{
        #browser()
        id_trt_work<- min(which(as.character(value$df_work[["arm"]]) == input$id_trt))
      }
      
      if(input$opinion_type == "mean"){
        id_comp_work<- min(which(as.character(value$df_work[["arm"]]) != input$id_trt))
        timepoint_expert_work <- NULL 
      }else{
        id_comp_work <- NULL
        timepoint_expert_work <- value$timepoint_expert
      }
      
      
      if(length(unique(value$df_work[["arm"]]))==1){
        formula_text <- "Surv(time,status)~1"
        id_trt_work<- id_comp_work<- id_St <- NULL
        
      }else{
        formula_text <- "Surv(time,status)~factor(arm)"
      }
      #browser()
      if(!is.null(value$param_expert)& input$opinion_type != "no_expert"){
        #browser()
        mod_fit  <- try({fit.models.expert(formula=as.formula(formula_text),data=value$df_work,
                                           distr=input$param_mod,
                                           method=input$stat_type,
                                           pool_type = input$pool_type_eval,#"log pool", 
                                           opinion_type = input$opinion_type,
                                           times_expert = timepoint_expert_work, 
                                           param_expert = value$param_expert,
                                           id_trt = id_trt_work,
                                           id_comp = id_comp_work,
                                           id_St  = id_trt_work,
                                           k = 1,
                                           iter = input$iterations,
                                           compile_mods = compile_mods)})
        value$mod_fit <- mod_fit
      }
      
      if(input$opinion_type == "no_expert"){
        
        mod_fit  <- fit.models.expert(formula=as.formula(formula_text),data=value$df_work,
                                      distr=input$param_mod,
                                      method=input$stat_type, k = 1)
        
        value$mod_fit <- mod_fit
      }
      
    })
    
    value$plot_km_expert1 <- eventReactive(value$mod_fit,{
      validate(need(class(value$mod_fit)!="try-error","Model estimation failed; Note Frequentist approach is typically much more fragile when expert opinion is in conflict with the data." ))
      
      if(input$incl_psa == "yes"){
        if(input$stat_type == "bayes"){
          nsim_eval <- input$iterations
        }else{
          nsim_eval <- 1000
        }
        plot_ci <- TRUE
        
      }else{
        nsim_eval <- 1
        plot_ci <- FALSE
      }
        
      if(input$opinion_type == "survival"){
          plot(value$mod_fit, add.km = TRUE,t = seq(0,input$xlim,length.out = 100),nsim  = nsim_eval, plot_ci = plot_ci)+
            geom_ribbon(data = value$df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                        fill = "sky blue", alpha = 0.5, colour = "grey")
          
        }else{
          plot(value$mod_fit, add.km = TRUE,t = seq(0,input$xlim,length.out = 100),nsim  = nsim_eval, plot_ci = plot_ci)
          
        } 
        
        
    })
    
    
    value$plot_gof <- eventReactive(value$mod_fit,{
      #browser()
      model.fit.plot(value$mod_fit,type = input$gof_type)
    })
    
    observeEvent(input$update_expert, {
      for (i in 1:input$n_timepoint) {
        local({
          j <- i  # Use a local variable to ensure correct value in the loop
          output[[paste0("expert_plot", j)]] <- renderPlot({
            value[[paste0("expert_plot", j)]]
          })
        })
      }
      shinyjs::hide("plot_gof")
    })
    
    
    observeEvent(input$run_analysis, {
      # browser()
      output$plot_km_expert1 <- renderPlot(
        value$plot_km_expert1())
      
      output$plot_gof <- renderPlot(
        value$plot_gof())
      
      shinyjs::show("plot_gof")
      
    })
    
    # observeEvent(input$save_output,{
    #   #browser()
    #   list_output <- list(model = value$mod_fit, surv_plt = value$plot_km_expert1(), gof_plt = value$plot_gof())
    #   
    #   
    #   #Need to fix this code chunk
    #   if(input$opinion_type != "no_expert"){
    #   for(i in 1:input$n_timepoint){
    #     list_output[[paste0("Timepoint",i)]] <- value[[paste0("expert_plot",i)]] 
    #   }
    #   }
    #   saveRDS(list_output,
    #           file = paste0(input$file_name,".rds"))
    #   
    #   #readRDS(file = paste0(input$file_name,".rds"))
    # })
    
    output$save_output <- downloadHandler(filename = paste0(input$file_name,".rds"),
                                          content = function(file){
                                            #browser()
                                            list_output <- list(model = value$mod_fit)
                                            #Need to fix this code chunk
                                            if(input$opinion_type != "no_expert"){
                                              for(i in 1:input$n_timepoint){
                                                list_output[[paste0("Timepoint",i)]] <- value[[paste0("expert_plot",i)]] 
                                              }
                                            }
                                            saveRDS(list_output,
                                                    file = file)
                                          })
    
    # observeEvent(input$save_output,{
    #   #browser()
    #   list_output <- list(model = value$mod_fit, surv_plt = value$plot_km_expert1(), gof_plt = value$plot_gof())
    #   
    #   
    #   #Need to fix this code chunk
    #   if(input$opinion_type != "no_expert"){
    #     for(i in 1:input$n_timepoint){
    #       list_output[[paste0("Timepoint",i)]] <- value[[paste0("expert_plot",i)]] 
    #     }
    #   }
    #   saveRDS(list_output,
    #           file = paste0(input$file_name,".rds"))
    #   
    #   #readRDS(file = paste0(input$file_name,".rds"))
    # })
    
    
    output$report <- downloadHandler(filename = function() {
      switch(input$outFormat, 
             html_document = paste0(input$file_name,"-report.html"), 
             pdf_document = paste0(input$file_name, "-report.pdf"),
             word_document = paste0(input$file_name,"-report.docx"))
    }, content = function(file) {
      
      #Will need these lines later
      # tempReport <- file.path(tempdir(), "elicitationShinySummary.Rmd")
      # file.copy(system.file("shinyAppFiles", "elicitationShinySummary.Rmd", 
      #                       package = "SHELF"), tempReport, overwrite = TRUE)
      
      list_output <- list(fit = value$mod_fit, 
                          surv_plt = value$plot_km_expert1(), 
                          gof_plt = value$plot_gof(),
                          xlim = input$xlim,
                          n_timepoint = input$n_timepoint,
                          opinion_type = input$opinion_type)
      
      if(input$opinion_type != "no_expert"){
        for(i in 1:input$n_timepoint){
          list_output[[paste0("expert_plot",i)]] <- value[[paste0("expert_plot",i)]] 
        }
      }
      
      
      
      file_pathway <- system.file("Report", "Generate-Report.Rmd", package = "expertsurv")
      if (file_pathway == "") {
        stop("File Generate-Report.Rmd not found.")
      }
      
      tempReport <- file.path(paste0(tempdir(), "\\Generate-Report.Rmd"))
      file.copy(file_pathway , tempReport, overwrite = TRUE)
      
      # rmarkdown::render(tempReport, output_file = file, #File Name
      #                   params = params, output_format = input$outFormat,
      #                   envir = new.env(parent = globalenv()))
          
      
      template <- use_parameters(tempReport, names(list_output),
                                 is.file = TRUE)
      writeLines(template,tempReport)
      
      rmarkdown::render(tempReport, output_file = file, #File Name
                        params = list_output, output_format = input$outFormat,
                        envir = new.env(parent = globalenv()))
      
      
    })
    
    
    
    
  }}
  

  shinyApp(ui = ui, server = create_server(compile_mods = compile_mods))
  



