########Function to calculate AMME

# This function calls the appropriate algorithm for estimating the MEMS. The dependent variable is assumed to be a dyad or dyad equivalent (e.g., dyad-time-period, dyad-group)
#
#
# micro_model is the model or list of models used to represent the micro process. A dyad or dyad equivalent (dyad time period, dyad group) is assumed to be the dependent variable.
# macro_model is the model used to predict a network-level or node-level outcome variable taking.
#micro_process is the micro process of interest, provided as a character string that matches the model output.
#mediator is the mediating variable of interest provided as a character string. It should match the output from the macro_model exactly
#link_id is a required vector of IDs used to link the micro_model output to the macro_model input. If calculating a network-level mediator, this should be the network identifier or network-group/network-time period identifier. If calculating a node-level mediator, this should be the node ID or node-time-period/node-group identifier. Observations should correspond exactly to rows in the macro_model data matrix. If calculating multiple network statistics at different levels of analysis, link_id may be provided as an ordered list of identifier, where each entry in the list is a vector of link_ids. If provided as a list, the first entry should correspond to the mediator and the remaining entries should correspond to the network controls.

#macro_function is the function to be calculated on the network object that provides the mediating variable. The function can be provided as the name of a function in the sna package to implement default settings. To include specific parameters within that function, write a user function of the form macro_function<-function(x){gtrans(x,measure="weakcensus)}. To call functions from igraph, igraph must be loaded in the global environment or called within the user function and object_type must be specified as object_type=="igraph". MEMS currently only implements functions for igraph and network objects, but that behavior can be overridden by specifying a user function that takes an igraph or network object as its input and computes relevant statistics from there. For example, the following user function calculates statistics from the netseg package for R assuming that x is an igraph object: macro_function=function(x){netseg::assort(x,"type")}. Also note that MEMS only takes into account the attributes included in the model. To include node or edge attributes not specified in the model, write a user function that takes a network or igraph object as its input and assigns the attribute. For example, if we want to calculate race assortativity using igraph, but race is not included in the model, we would specify network_object="igraph" and the user function: macro_function=function(x,race){require(igraph); V(x)$race<-race;assortativity(x,V(x)$race)}. Whatever macro function is specified, the output should either be a numeric value or a vector of numeric values for the algorithm to proceed.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object. May be provided as a vector when calculating controls using distinct object types. If provided as a vector, the first entry should correspond to the mediator and the remaining entries should correspond to the controls.
#control is an optional character string or vector of character strings that provide the named entries for each network control variable in the macro_model. Entries should correspond exactly to macro_model output. If left NULL, not controls are computed.
#control_functions is a list of functions to calculate on the simulated macro networks. These may either be user functions or functions that inherent from igraph or statnet. The list should correspond exactly to the elements listed in "controls"
#the interval provides the values over which to vary the micro-process It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#algorithm tells netmediate which estimation routine to use. Parametric and nonparametric are available. Defaults to parametric.
#full_output tells R whether to return simulated distribution in addition to summary statistic
#time interval is an optional parameter to be used with rem.dyad models or GLM/GLMER models for REM data. Time interval tells which time points to compute the MEMS. Two options are available: aggregate or a continuous vector of numeric values. If "aggregate" is specified, a single MEMS is calculated by creating an aggregate cross-sectional representation of the event sequence treating the final observed time of event activity as the limiting condition for simulated event sequences. If a continuous interval is specified of the form c(0,2, 3), unique networks are created containing all ties that occur between those intervals. For example, if c(1,2,5) is specified, unique networks are constructed using all ties that occur between the 1 and 2 time periods and the 2 and 5 time periods, and the MEMS is computed over those two networks. If unspecified, defaults to all time periods.
#covar_list is a list of covariates used to estimate the MEMS when using a rem.dyad object.
#edgelist is an edgelist that acts as the dependent variable in rem.dyad when using a rem.dyad object. Otherwise left empty.
#net_logit_y and net_logit_x are the dependent and independent variables specified using network logistic regression with the netlogit function in sna. y must be a binary adjacency matrix and x stack of independent network variables.
#group_id is an optional vector of group identifiers to use when estimating a glm or glmer on grouped data (i.e., multiple time periods, multiple networks) on a dependent tie variable. When specified, the MEMS command will induce unique networks for each grouping factor. If left unspecified, all groups/time period are pooled. If using glmer, the grouping factor does not have to be provided as part of the model or used as a random effect.
#node numbers is a vector listing the number of nodes in each network when using GLM or GLMER. If estimating MEMS aggregated over all networks, this shoud be the total number of nodes in all networks. Required when using GLM or GLMER, ignored otherwise.

AMME_nonparam<-function(micro_model=micro_model,
                     macro_model=macro_model,
                     micro_process=micro_process,
                     mediator=mediator,
                     macro_function=macro_function,
                     link_id=link_id,
                     object_type=object_type,
                     controls=controls,
                     control_functions=control_functions,
                     interval=interval,
                     nsim=nsim,
                     algorithm=algorithm,
                     silent=silent,
                     full_output=full_output,
                     SAOM_data=SAOM_data,
                     time_interval=time_interval,
                     covar_list=covar_list,
                     edgelist=edgelist,
                     net_logit_y=net_logit_y,
                     net_logit_x=net_logit_x,
                     group_id=group_id,
                     node_numbers=node_numbers){

  if(class(macro_model)[1]%in%"plm"){
    stop("PLM objects not supported for nonparametric estimation. Try parametric estimation instead.")
  }

  if(class(macro_model)[1]%in%c("lnam")){
    link_fun<-list(family="gaussian",link="identity")
  }else{
    link_fun<-family(macro_model)
  }
  if(!link_fun$link%in%c("logit","probit","identity")){
    stop("netmediate currently only supports logit, probit, and linear GLMs for parametric AMME estimation.")
  }

  ##
  if(class(link_id)[1]!="list"){
    ID_vec<-link_id
    link_id<-list()
    for(i in 1:(length(controls)+1)){
      link_id[[i]]<-ID_vec
    }
  }

  if(length(link_id)<length(controls)){

    stop("Dimensions of link_id do not match number of controls specified.")

  }


  ####step 1: get MEMS output
  message("Getting output for micro-macro relationship.")
  MEMS_output<-MEMS(model=micro_model,
                    micro_process=micro_process,
                    macro_function=macro_function,
                    link_id=link_id,
                    object_type=object_type,
                    controls=controls,
                    control_functions=control_functions,
                    interval=interval,
                    nsim=nsim,
                    algorithm=algorithm,
                    silent=silent,
                    full_output=TRUE,
                    SAOM_data=SAOM_data,
                    time_interval=time_interval,
                    covar_list=covar_list,
                    edgelist=edgelist,
                    net_logit_y=net_logit_y,
                    net_logit_x=net_logit_x,
                    group_id=group_id,
                    node_numbers=node_numbers)



  #step 2: start bootstrap
  if(class(macro_model)[1]%in%c("glmerMod","lmerMod")){
    coef<-lme4::fixef(macro_model)

  }else{
    coef<-coef(macro_model)
  }
  if(!mediator%in%names(coef)){
    stop("Mediator does not match model output. Check spelling and respecify.")
  }
  if(!is.null(controls) & !any(controls%in%names(coef))){
    stop("One more controls do not match model output. Check spelling and respecify.")
  }


  if(any(is.na(coef))|
     any(is.infinite(coef))){
    stop("Infinite or missing values in parameter estimates. Algorithm cannot continue.")

  }
  interval<-sort(interval) #order from lowest to highest
  theta<-matrix(NA,nrow=nsim,ncol=length(coef))
  colnames(theta)<-names(coef)


  if(class(macro_model)[1]=="lnam"){

    model_frame<-cbind.data.frame(macro_model$y,
                                macro_model$x)
    colnames(model_frame)<-c("y",colnames(macro_model$x))
    if(!is.null(macro_model$W1)){
      model_frame$rho1.1<-macro_model$W1[,,]%*%macro_model$y
    }
    if(!is.null(macro_model$W2)){
      model_frame$rho2.1<-macro_model$W2[,,]%*%macro_model$residuals
    }


  }else{
    model_frame<-model.frame(macro_model)
  }
  if(class(macro_model)[1]=="plm"){
    index<-attributes(macro_model$model)$index
    model_frame<-cbind(model_frame,index)
  }

  message("Starting bootstrap estimation of AMME.")
  for(j in 1:nsim){

    bs_mat<-model_frame[sample(nrow(model_frame), replace=TRUE), ]

    sampling_error<-1
    tracer<-1
    class(sampling_error)<-"try-error"

    if(class(macro_model)[1]%in%c("glm","lm")){
      bs_model<-glm(bs_mat[,1] ~  1+., data = data.frame(bs_mat[,-c(1)]),family=link_fun$family)
    }
    if(class(macro_model)[1]%in%c("lnam")){
      bs_model<-glm(bs_mat[,1] ~  -1+., data = data.frame(bs_mat[,-c(1)]),family=link_fun$family)
    }

    if(class(macro_model)[1]%in%c("glmerMod","lmerMod")){
      while(class(sampling_error)[1]=="try-error"){

        bs_mat<-model_frame[sample(nrow(model_frame), replace=TRUE), ]
        bs_model<-try(lme4::glmer(formula(macro_model), data = data.frame(bs_mat),family=link_fun$family),silent=TRUE)
        class(sampling_error)<-class(bs_model)
        tracer<-tracer+1

        if(tracer>=100){
          stop("Unable to estimate GAM in a sequence of 1000 bootstrap samples. Estimation may not be feasible with these data.")
          }
        }
      }

    if(class(macro_model)[1]%in%"Gam"){


      while(class(sampling_error)[1]=="try-error"){
        bs_mat<-macro_model$data
        bs_mat<-bs_mat[sample(nrow(bs_mat),replace=TRUE),]
        bs_model<-try(gam::gam(formula(macro_model), data = data.frame(bs_mat),family=link_fun$family),silent = TRUE)
        class(sampling_error)<-class(bs_model)
        tracer<-tracer+1

        if(tracer>=100){
          stop("Unable to estimate GAM in a sequence of 1000 bootstrap samples. Estimation may not be feasible with these data.")
        }
      }

    }

    if(class(macro_model)[1]%in%"plm"){
      formula<-formula(macro_model)
      bs_model<-plm::plm(formula,data=bs_mat,
                         model=macro_model$args$model,
                         effect=macro_model$args$effect,
                         index=colnames(bs_mat)[c(ncol(bs_mat)-1,ncol(bs_mat))])
    }

    if(class(macro_model)[1]%in%c("glmerMod","lmerMod")){
      theta[j,colnames(theta)%in%names(lme4::fixef(bs_model))]<-lme4::fixef(bs_model)
    }else{
      theta[j,colnames(theta)%in%names(coef(bs_model))]<-coef(bs_model)
    }


    if(silent==FALSE){

      message("Getting bootstrap parameters ", j, " of ", nsim, " for macro_model.")

    }


  }


theta[is.na(theta)]<-0
  ##step 3 set up parameters to estimate AMME from macro model

  if(class(macro_model)[1]=="lnam"){

    #no intercept with LNAM, so DV can be ignored
    model_mat<-data.frame(macro_model$x)
    if(!is.null(macro_model$W1)){
      model_mat$rho1.1<-macro_model$W1[,,]%*%macro_model$y
    }
    if(!is.null(macro_model$W2)){
      model_mat$rho2.1<-macro_model$W2[,,]%*%macro_model$residuals
    }
    model_mat<-as.matrix(model_mat)

  }else{

    model_mat<-stats::model.matrix(macro_model)

  }

  if(class(macro_model)[1]%in%c("glmerMod","lmerMod")){
    model_mat<-model_mat[,!colnames(model_mat)%in%names(lme4::ranef(macro_model))]
  }


   ##step 4: start simulation

  ##now need to loop over parameter vector and predict
  #so, need to be able to match MEMS output to parameter iteration
  message("Beginning AMME estimation.")
  #get change stat values
  X_vals<-length(interval)
  M_vals<-length(interval)
  AMME_vec<-vector(length=nsim)# output data
  prop_vec<-vector(length=nsim) #proportion explained data
  crosswalk<-list() #only used for nested data
  for(i in 1:length(link_id)){
    crosswalk[[i]]<-data.frame(link_id=link_id[[i]])
  }

  for(j in 1:nsim){
    sim_model_mat<-model_mat

    #get predictions at each X by M interval
    for(i in 1:X_vals){

      k_vec<-matrix(NA,nrow=nrow(sim_model_mat),ncol=M_vals)

      for(m in 1:M_vals){ #loop over m intervals

        #use merges if observations don't align (i.e., nested data)
        if(nrow(sim_model_mat)!=nrow(MEMS_output$out_dat_main[[j]])){
          sim_dat<-data.frame(link_id[[1]]==rownames(MEMS_output$out_dat_main[[j]]),
                              M=MEMS_output$out_dat_main[[j]][,m])
          sim_dat<-suppressMessages(plyr::join(sim_dat,crosswalk[[1]]))
          sim_model_mat[,mediator]<-sim_dat$M

        }else{
          sim_model_mat[,mediator]<-MEMS_output$out_dat_main[[j]][,m]
        }

        if(!is.null(controls)){
          for(k in 1:length(MEMS_output$out_dat_controls)){
            control<-names(MEMS_output$out_dat_controls)[k]

            if(nrow(sim_model_mat)!=nrow(MEMS_output$out_dat_controls[[k]][[j]])){
              sim_dat<-data.frame(link_id[[k+1]]==rownames(MEMS_output$out_dat_main[[j]]),
                                  M=MEMS_output$out_dat_main[[j]][,m])
              sim_dat<-suppressMessages(plyr::join(sim_dat,crosswalk[[k+1]]))
              sim_model_mat[,control]<-sim_dat$M

            }else{
              sim_model_mat[,control]<-MEMS_output$out_dat_controls[[k]][[j]][,i] #for control k, simulation j, interval i
            }
          }
        }

        #get prediction
        lp<-as.matrix(sim_model_mat)%*%theta[j,]

        if(link_fun$link=="logit"){
          lp <- 1/(1 + exp(-lp))
        }

        if(link_fun$link=="probit"){
          lp<-VGAM::probitlink(lp,inverse=TRUE)
        }

        k_vec[,m]<-lp


      }#close M_vals loop

      #calculate differences in columns
      change_mat<-as.matrix(k_vec[,-c(1)])
      prop_mat<-change_mat

      for(k in 1:ncol(change_mat)){
        upper<-k+1
        lower<-k
        change_mat[,k]<-k_vec[,upper]-k_vec[,lower]
        prop_mat[,k]<-change_mat[,k]/k_vec[,upper]

      }

      #store all change stats
      if(i ==1){
        x_vec<-change_mat
      }else{
        x_vec<-cbind(x_vec,change_mat)
      }


    } #close x_vals loop

    #compute mean and store result
    AMME_vec[j]<-mean(x_vec,na.rm=T)
    prop_mat<-prop_mat[!is.infinite(prop_mat),]
    prop_vec[j]<-mean(prop_mat,na.rm=T)

    if(silent==FALSE){
      message("Simulation ", j, " of ",nsim," complete.")
    }
  }#close j loop



  ##package output
  summary_dat<-matrix(NA,nrow=2,ncol=5)
  rownames(summary_dat)<-c("AMME","Prop. Change in M")
  colnames(summary_dat)<-c("Estimate","Std. Dev.","lower 95% CI","Upper 95% CI","BS p-val")

  summary_dat[1,1]<-mean(AMME_vec,na.rm=TRUE)
  summary_dat[1,2]<-sd(AMME_vec,na.rm=TRUE)

  if(summary_dat[1,1]<0){
    summary_dat[1,5]<-length(AMME_vec[which(AMME_vec>=0)])/nsim
  }else{
    summary_dat[1,5]<-length(AMME_vec[which(AMME_vec<=0)])/nsim

  }
  summary_dat[1,3]<-quantile(AMME_vec,.025,na.rm=TRUE)
  summary_dat[1,4]<-quantile(AMME_vec,.975,na.rm=TRUE)
  summary_dat[2,1]<-mean(prop_vec,na.rm=TRUE)


  if(full_output==FALSE){
    return(summary_dat)
  }else{
    out_dat<-list(summary_dat=summary_dat,
                  AMME_obs=AMME_vec,
                  prop_explained_obs=prop_vec)
    return(out_dat)

  }



}#close function
