#' Plot the survival curves using \code{ggplot2}
#' 
#' @param exArgs list of extra options passed to \code{plot.survHE}. These 
#' include whether the KM curve should be added \code{add.km} and whether
#' the user specifies a profile of covariates (in the list \code{newdata}).
#' Other possibilities are additional (mainly graphical) options. 
#' These are: \code{xlab} = a string with the label for the
#' x-axis (default = "time") \code{ylab} = a string with the label for the
#' y-axis (default = "Survival") \code{lab.profile} = a (vector of) string(s)
#' indicating the labels associated with the strata defining the different
#' survival curves to plot. Default to the value used by the Kaplan Meier
#' estimate given in \code{fit.models}. \code{newdata} = a list (of lists) 
#' providing the values for the relevant covariates If NULL, then will use 
#' the mean values for the covariates if at least one is a continuous variable, 
#' or the combination of the categorical covariates. \code{xlim} = a vector 
#' determining the limits for the x-axis \code{colors} = a vector of characters 
#' defining the colours in which to plot the different survival curves 
#' \code{lab.profile} = a vector of characters defining the names of the models fitted 
#' \code{add.km} = TRUE (whether to also add the Kaplan Meier estimates of the data) 
#' \code{annotate} = FALSE (whether to also add text to highlight the observed vs
#' extrapolated data)
#' \code{legend.position} = a vector of proportions to place the legend. Default
#' to 'c(.75,.9)', which means 75% across the x-axis and 90% across the y-axis
#' \code{legend.title} = suitable instructions to format the title of the legend;
#' defaults to 'element_text(size=15,face="bold")' but there may be other 
#' arguments that can be added (using 'ggplot' facilities)
#' \code{legend.text} = suitable instructions to format the text of the legend;
#' defaults to 'element_text(colour="black", size=14, face="plain")' but there 
#' may be other arguments that can be added (using 'ggplot' facilities)
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Parametric survival models
#' @import dplyr
#' @import ggplot2
#' @examples
#' #' 
#' data(bc)
#' 
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' plot(mle)
#' @noRd 
plot_ggplot_expertsurv <- function(exArgs) {
  
  # # First checks the class of the input
  # if(class(x)=="survHE") {
  #   # If x is a 'survHE' object, then there's only one object to deal with
  #   surv.curv=make_surv_curve_plot(x,mods,nsim=1,t,newdata,add.km)
  # }
  
  # Extracts the 'survHE' objects from the list 'exArgs. If there are none, then stop with an error message!
  w <- which(unlist(lapply(1:length(exArgs),function(i) class(exArgs[[i]])))=="expertsurv")
  if(length(w)==0){
    stop("Please give at least one 'survHE' object, generated by a call to 'fit.models(...)")
  } else {
    survHE_objs=lapply(1:length(w),function(i) exArgs[[w[i]]])
  }
  names(survHE_objs)=names(exArgs)[w]
  
  # Check some basic inputs - if not specified by the user (and stored in 'exArgs', then
  # sets up defaults to be used in the other functions)
  # t = time to plot on the x-axis (default = the times in the original data/survHE object)
  if (!exists("t",exArgs)) {t <- sort(unique(survHE_objs[[1]]$misc$km$time))} else {t <- exArgs$t}
  # newdata = possible list with a profile of covariates (default = NULL)
  if (!exists("newdata",exArgs)) {newdata <- NULL} else {newdata <- exArgs$newdata}
  
  # Do we want a single survival curve, or is this for the PSA?
  if (!exists("nsim",exArgs)) {nsim <- 1} else {nsim <- exArgs$nsim}

  # What model should be used from the object 'x'?
  if (!exists("mods",exArgs)) {
    mods=1:sum(unlist(lapply(survHE_objs,function(x) length(x$models))))
  } else {mods=exArgs$mods}
  # Should the KM curve be added to the plots?
  if (!exists("add.km",exArgs)) {add.km <- FALSE} else {add.km <- exArgs$add.km}
  
  # Should the graph be annotated with extrapolation vs observed data?
  if(exists("annotate",where=exArgs)){annotate=exArgs$annotate} else {annotate=FALSE}
  
  # Makes the dataframe with the data to plot
  # toplot = lapply(1:length(survHE_objs),function(i){
  #   make_data_surv(survHE_objs[[i]],
  #                  mods=1:length(survHE_objs[[i]]$models),
  #                  nsim=nsim,
  #                  t=t,
  #                  newdata=newdata,
  #                  add.km=add.km
  #   )[[1]] %>% mutate(object_name=as.factor(names(survHE_objs)[i]))
  # }) %>% bind_rows() %>%
  #   group_by(object_name,model_name) %>% mutate(mods_id=cur_group_id()) %>% ungroup() %>%
  #   filter(mods_id%in%mods)
  
  ##############################################################################################
  # Tries to only select the relevant models based on the choice indicated by the user
  # Makes a tibble with the *only* objects + the models selected in each of them
  # Initialises 'obj' and 'mod' to avoid binding issues
  obj <- mod <- NULL
  all_models=tibble(
    obj=unlist(
      lapply(1:length(survHE_objs),function(x) {
        rep(names(survHE_objs)[x],length(survHE_objs[[x]]$models))
      })
    ),
    mod=unlist(lapply(survHE_objs,function(x) 1:length(x$models)))
  ) %>% slice(mods) %>% arrange(obj)
  
  # Makes a vector with the index of *only* the objects selected
  sel_mods=unique(match(all_models$obj,names(survHE_objs)))
  
  # Makes the dataset to plot, including *only* the objects and models selected
  toplot = lapply(sel_mods,function(i){
    make_data_surv(survHE_objs[[i]],
                   mods=all_models %>% filter(obj==names(survHE_objs)[i]) %>% pull(mod), 
                   nsim=nsim,
                   t=t,
                   newdata=newdata,
                   add.km=add.km 
    )[[1]] %>% mutate(object_name=as.factor(names(survHE_objs)[i]))
  }) %>% bind_rows() %>% 
    group_by(object_name,model_name) %>% mutate(mods_id=cur_group_id()) %>% ungroup() 
  ##############################################################################################
  
  # If so, then builds the relevant data
  if(add.km==TRUE) {
    datakm = lapply(1:length(survHE_objs),function(i){
      make_data_surv(survHE_objs[[i]],
                     mods=1, #1:length(survHE_objs[[i]]$models), 
                     nsim=1,
                     t=t,
                     newdata=newdata,
                     add.km=add.km
      )[[2]] %>% mutate(object_name=as.factor(names(survHE_objs)[i]))
    }) %>% bind_rows() %>% 
      group_by(object_name,model_name) %>% mutate(mods_id=cur_group_id()) %>% ungroup()
  } else {
    datakm=NULL
  }

  # Now makes the plot using the helper function
  surv.curv=make_surv_curve_plot(toplot,datakm,mods)

  # Optional arguments
  if(exists("lab.profile",exArgs)){
    surv.curv=surv.curv+
      scale_linetype_manual(labels=exArgs$lab.profile,values=1:length(exArgs$lab.profile))
  }
  # if both colours & labels are specified for the models chosen
  if(exists("colour",exArgs) & exists("lab.model",exArgs)) {
    surv.curv=surv.curv+scale_color_manual(labels=exArgs$lab.model,values=exArgs$colour)
  }
  # if only the colours
  if(exists("colour",exArgs) & !exists("lab.model",exArgs)) {
    surv.curv=surv.curv+scale_color_manual(values=exArgs$colour)
  }
  # if only the labels
  if(exists("lab.model",exArgs) & !exists("colour",exArgs)) {
    surv.curv=surv.curv+scale_color_manual(values=1:length(exArgs$lab.model),labels=exArgs$lab.model)
  }
  if(exists("xlab",where=exArgs)){
    surv.curv=surv.curv+labs(x=exArgs$xlab)
  }
  if(exists("ylab",where=exArgs)){
    surv.curv=surv.curv+labs(y=exArgs$ylab)
  }
  if(exists("main",where=exArgs)) {
    surv.curv=surv.curv+labs(title=exArgs$main)+theme(plot.title=element_text(size=18,face="bold"))
  }
  if(annotate==TRUE){
    cutoff=max(survHE_objs[[1]]$misc$km$time)
    surv.curv=surv.curv + #geom_vline(xintercept=cutoff,linetype="dashed",size=1.5) +
      geom_segment(aes(x=cutoff,y=-Inf,xend=cutoff,yend=-.01),size=0.9) + 
      geom_segment(aes(x=cutoff,y=-.01,xend=cutoff*.85,yend=-.01),arrow=arrow(length=unit(.25,"cm"),type="closed"),size=1.1)+
      geom_segment(aes(x=cutoff,y=-.01,xend=cutoff*1.15,yend=-.01),arrow=arrow(length=unit(.25,"cm"),type="closed"),size=1.1)+
      annotate(geom="text",x=cutoff,y=-Inf,hjust=1.1,vjust=-1,label="Observed data",size=5) +
      annotate(geom="text",x=cutoff,y=-Inf,hjust=-0.1,vjust=-1,label="Extrapolation",size=5) +
      ylim(-0.01,1) + 
      geom_rect(data=data.frame(xmin=-Inf,xmax=cutoff,ymin=-Inf,ymax=Inf),
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey",alpha=.1)
  } else{
    surv.curv=surv.curv+ylim(0,1)
  }
  if(exists("legend.position",exArgs)){
    surv.curv=surv.curv+theme(legend.position=exArgs$legend.position)
  }
  if(exists("legend.title",exArgs)){
    surv.curv=surv.curv+theme(legend.title=exArgs$legend.title)
  }
  if(exists("legend.text",exArgs)){
    surv.curv=surv.curv+theme(legend.text=exArgs$legend.text)
  }
  # to remove the profiles legend
  #surv.curv=surv.curv+guides(linetype=FALSE)
  # to modify the profile legend
  #surv.curv=surv.curv+scale_linetype_discrete(name="XXX",label=c("XX","YY","ZZ"))
  # to remove the models legend
  #surv.curv=surv.curv+guides(colour=FALSE)
  # to modify the profile legend
  #surv.curv=surv.curv+scale_color_discrete(name="XXX",label=c("XX","YY","ZZ"))
  # +scale_linetype_manual(labels=c("Control","Treated"),values=c("dotdash","solid"))

  # Now prints the plot
  surv.curv
}


#' Make the dataset to be used by \code{ggplot2} to plot the survival curves
#' 
#' @param x The 'survHE' object
#' @param mods The models to be considered
#' @param nsim The number of simulations to generate
#' @param t The vector of times
#' @param newdata The list of "new" covariares proffiles
#' @param add.km Should the KM estimate be plotted too?
#' @return \item{surv.curv}{The \code{ggplot2} object with the graph}
#' @note Something will go here
#' @author Gianluca Baio
#' @import tibble
#' @import dplyr
#' @keywords Parametric survival models
#' @noRd 
make_data_surv <- function(x,mods=1:length(x$models),nsim=1,t=NULL,newdata=NULL,add.km=FALSE) {
  if(is.null(t)) {
    t <- sort(unique(x$misc$km$time))
  }
  #s=lapply(1:length(x$models),function(i) {
  s=lapply(mods,function(i) {
    make.surv(x,mod=i,t=t,nsim=nsim,newdata=newdata)
  })
  strata=lapply(1:length(s),function(i) {
    lapply(1:nrow(s[[i]]$des.mat),function(x){
      s[[i]]$des.mat %>% as_tibble() %>% select(-matches("(Intercept)",everything())) %>% slice(x) %>% 
        round(digits=2) %>% mutate(strata=paste0(names(.),"=",.,collapse=","))
    }) %>% bind_rows(.) %>% select(strata)
  })

  # toplot=lapply(1:length(mods),function(i) {
  #   lapply(1:length(s[[mods[i]]]$S),function(j) {
  #     s[[mods[i]]]$S[[j]] %>% bind_cols(strata=as.factor(strata[[mods[i]]][j,]),model_name=as.factor(names(x$models)[mods[i]]))
  #   })
  # }) %>% bind_rows(.)
  # out=list(toplot)
  toplot=lapply(1:length(mods),function(i) {
    lapply(1:length(s[[i]]$S),function(j) {
      s[[i]]$S[[j]] %>% bind_cols(strata=as.factor(as.character(strata[[i]][j,])),model_name=as.factor(names(x$models)[mods[i]]))
    })
  }) %>% bind_rows(.)
  out=list(toplot)
  
  # Add the data for the KM curve?
  if(add.km==TRUE) {
    # If the number of strata in the KM computed in 'fit.models' is not the same as the 
    # number of rows in the design matrix from 'make.surv', then re-do a KM with no covariates
    if(length(x$misc$km$strata)!=nrow(s[[1]]$des.mat)){
      x$misc$km=rms::npsurv(update(x$misc$formula,~1),data=x$misc$data)
      x$misc$km$call$formula=as.formula(deparse(update(x$misc$formula,~1)))
    }
    # Now uses info in the KM table in the survHE object to create a dataset to plot
    datakm=bind_cols(t=x$misc$km$time,n.risk=x$misc$km$n.risk,n.event=x$misc$km$n.event,
                     n.censor=x$misc$km$n.censor,S=x$misc$km$surv,lower=x$misc$km$lower,
                     upper=x$misc$km$upper) %>% mutate(model_name="Kaplan Meier")
    # If 'strata' is not in the KM object (will happen if there's only 1)
    if(is.null(x$misc$km$strata)) {
      datakm$strata=as.factor("all")
    } else {
      datakm$strata=as.factor(rep(1:length(x$misc$km$strata),x$misc$km$strata))
    }
    out$datakm=datakm
  }
  # Returns the output as a list with the dataset(s) to plot
  return(out)
}


#' Make the actual \code{ggplot2} plot with the survival curves
#' 
#' @param toplot The dataset with the relevant data
#' @param dataKM The dataset with the (optional) data for the KM estimate
#' @param mods The models to be plotted (a vector of numbers)
#' @return \item{out}{A list with the dataset to be plotted including the survival curves}
#' @import ggplot2
#' @note Something will go here
#' @author Gianluca Baio
#' @keywords Parametric survival models
#' @noRd 
make_surv_curve_plot <- function(toplot,datakm=NULL,mods) {
  # Does the model have covariates?
  if (all(toplot$strata=="=")) {
    # In this case not (intercept only), so remove the linetype as not needed
    linetype=NULL
  } else {
    # If it does have covariates then use 'strata' to plot a curve per profile
    linetype=toplot$strata
  }
  surv.curv=ggplot() 
  # Am I plotting a single 'survHE' object?
  if(length(levels(toplot$object_name))==1) {
    surv.curv=surv.curv+
      geom_line(data=toplot,aes(x=t,y=S,group=model_name:strata,col=model_name,linetype=linetype),size=.9) 
  } else {
    surv.curv=surv.curv+
      geom_line(data=toplot,aes(x=t,y=S,group=model_name:strata:object_name,col=object_name:model_name,linetype=linetype),size=.9)   
  }
  surv.curv=surv.curv +
    theme_bw() + 
    theme(axis.text.x = element_text(color="black",size=12,angle=0,hjust=.5,vjust=.5),
          axis.text.y = element_text(color="black",size=12,angle=0,hjust=.5,vjust=.5),
          axis.title.x = element_text(color="black",size=14,angle=0,hjust=.5,vjust=.5),
          axis.title.y = element_text(color="black",size=14,angle=90,hjust=.5,vjust=.5)) +
    theme(axis.line = element_line(colour = "black"),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size=18, face="bold")) +
    theme(legend.position=c(.75,.78),
          legend.title=element_text(size=15,face="bold"),
          #legend.title = element_blank(),
          legend.text = element_text(colour="black", size=14, face="plain"),
          legend.background=element_blank()) +
    labs(y="Survival",x="Time",title=NULL,
         color=ifelse(length(mods)==1,"Model","Models"),
         linetype="Profile") + 
    # This ensures that the model legend is always before the profile legend
    guides(color=guide_legend(order=1),
           linetype=guide_legend(order=2))
  # If uses more than 1 simulation from distribution of survival curves, then add ribbon
  if(any(grepl("low",names(toplot)))) {
    surv.curv=surv.curv+geom_ribbon(data=toplot,aes(x=t,y=S,ymin=low,ymax=upp,group=model_name:strata),alpha=.2)
  }
  
  # Add KM plot? 
  if(!is.null(datakm)) {
    surv.curv=surv.curv+geom_step(data=datakm,aes(t,S,group=as.factor(strata)),color="darkgrey") + 
      geom_ribbon(data=datakm,aes(x=t,y=S,ymin=lower,ymax=upper,group=as.factor(strata)),alpha=.2) 
  }
  surv.curv
}

