#'Function to calculate AMME
#'This function calls the appropriate algorithm for estimating the MEMS. The dependent variable is assumed to be a dyad or dyad equivalent (e.g., dyad-time-period, dyad-group)
#'
#'
#'@micro_model is the model or list of models used to represent the micro process. A dyad or dyad equivalent (dyad time period, dyad group) is assumed to be the dependent variable.
#'@macro_model is the model used to predict a network-level or node-level outcome variable taking.
#'@micro_process is the micro process of interest, provided as a character string that matches the model output.
#'@mediator is the mediating variable of interest provided as a character string. It should match the output from the macro_model exactly
#'@link_id is a required vector of IDs used to link the micro_model output to the macro_model input. Link_id elements should correspond to unique observations in the macro_model data matrix. If calculating a network-level mediator, this should be the network identifier or network-group/network-time period identifier. If calculating a node-level mediator, this should be the node ID or node-time-period/node-group identifier. If calculating multiple network statistics at different levels of analysis, link_id may be provided as an ordered list of identifier, where each entry in the list is a vector of link_ids. If provided as a list, the first entry should correspond to the mediator and the remaining entries should correspond to the network controls.
#'@macro_function is the function to be calculated on the network object that provides the mediating variable. The function can be provided as the name of a function in the sna package to implement default settings. To include specific parameters within that function, write a user function of the form macro_function<-function(x){gtrans(x,measure="weakcensus)}. To call functions from igraph, igraph must be loaded in the global environment or called within the user function and object_type must be specified as object_type=="igraph". MEMS currently only implements functions for igraph and network objects, but that behavior can be overridden by specifying a user function that takes an igraph or network object as its input and computes relevant statistics from there. For example, the following user function calculates statistics from the netseg package for R assuming that x is an igraph object: macro_function=function(x){netseg::assort(x,"type")}. Also note that MEMS only takes into account the attributes included in the model. To include node or edge attributes not specified in the model, write a user function that takes a network or igraph object as its input and assigns the attribute. For example, if we want to calculate race assortativity using igraph, but race is not included in the model, we would specify network_object="igraph" and the user function: macro_function=function(x,race){require(igraph); V(x)$race<-race;assortativity(x,V(x)$race)}. Whatever macro function is specified, the output should either be a numeric value or a vector of numeric values for the algorithm to proceed.
#'@object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object. May be provided as a vector when calculating controls using distinct object types. If provided as a vector, the first entry should correspond to the mediator and the remaining entries should correspond to the controls.
#'@control is an optional character string or vector of character strings that provide the named entries for each network control variable in the macro_model. Entries should correspond exactly to macro_model output. If left NULL, not controls are computed.
#'@control_functions is a list of functions to calculate on the simulated macro networks. These may either be user functions or functions that inherent from igraph or statnet. The list should correspond exactly to the elements listed in "controls"
#'@interval provides the values over which to vary the micro-process It defaults to 0 and 1
#'@nsim is the number of Monte Carlo samples
#'@silent tells R whether to provide updates on the algorithm's progress
#'@algorithm tells netmediate which estimation routine to use. Parametric and nonparametric are available. Defaults to parametric.
#'@full_output tells R whether to return simulated distribution in addition to summary statistic
#'@SAOM_data is a siena object that contains the data to use for SAOM estimation. Only required when using SAOM to calculate the MEMS.
#'@SAOM_var is an optional list of varCovar and varDyadCovar to calculate on each simulated network. Only required when the macro statistic is a funciton of time varying node or dyad attributes. CoCovar and coDyadCovar (time invariant node covariates and dyad covariates) objects are handled internally. List must be provided with named entries if used in a call to AMME with a SAOM micro-model.
#'@time_interval is an optional parameter to be used with rem.dyad models or GLM/GLMER models for REM data. Time interval tells which time points to compute the MEMS. Two options are available: aggregate or a continuous vector of numeric values. If "aggregate" is specified, a single MEMS is calculated by creating an aggregate cross-sectional representation of the event sequence treating the final observed time of event activity as the limiting condition for simulated event sequences. If a continuous interval is specified of the form c(0,2, 3), unique networks are created containing all ties that occur between those intervals. For example, if c(1,2,5) is specified, unique networks are constructed using all ties that occur between the 1 and 2 time periods and the 2 and 5 time periods, and the MEMS is computed over those two networks. If unspecified, defaults to all time periods. When estimating the AMME, the time_interval is assumed to correspond to the temporal component of the macro_model.
#'@covar_list is a list of covariates used to estimate the MEMS when using a rem.dyad object.
#'@edgelist is an edgelist that acts as the dependent variable in rem.dyad when using a rem.dyad object. Otherwise left empty.
#'@net_logit_y and net_logit_x are the dependent and independent variables specified using network logistic regression with the netlogit function in sna. y must be a binary adjacency matrix and x stack of independent network variables.
#'@group_id is an optional vector of group identifiers to use when estimating a glm or glmer on grouped data (i.e., multiple time periods, multiple networks) on a dependent tie variable. When specified, the MEMS command will induce unique networks for each grouping factor. If left unspecified, all groups/time period are pooled. If using glmer, the grouping factor does not have to be provided as part of the model or used as a random effect.
#'@node_numbers is a vector listing the number of nodes in each network when using GLM or GLMER. If estimating MEMS aggregated over all networks, this shoud be the total number of nodes in all networks. Required when using GLM or GLMER, ignored otherwise.


AMME<-function(micro_model,
               macro_model,
               micro_process,
               mediator,
               macro_function,
               link_id,
               object_type=NULL,
               controls=NULL,
               control_functions=NULL,
               interval=c(0,1),
               nsim=500,
               algorithm="parametric",
               silent=FALSE,
               full_output=FALSE,
               SAOM_data=NULL,
               SAOM_var=NULL,
               time_interval=NULL,
               covar_list=NULL,
               edgelist=NULL,
               net_logit_y=NULL,
               net_logit_x=NULL,
               group_id=NULL,
               node_numbers=NULL){



  if(!class(macro_model)[1]%in%c("lm","glm","glmerMod","lmerMod","Gam","lnam","plm")){
    stop("Macro_model type not currently supported. AMME only takes lm, glm, glmerMod, lmerMod, gam, lnam, and plm type objects.")
  }

  if(length(link_id)<2 & class(link_id)[1]!="list"){
    stop("Link_id must have length>=1 to estimate the AMME.")
  }



  if(!is.null(controls) & length(controls)!=length(control_functions)){
    stop("Elements in controls do not match elements in control_functions exactly. Please respecify.")
  }

  if(class(control_functions)[1]=="function"){
    control_functions<-list(control_functions)
  }

  if(!is.null(controls)& class(control_functions)[1]!="list"){
    stop("control_functions must be provided as a list of functions.")
  }

  if(length(object_type)==1 & !is.null(controls)){
    object_type<-rep(object_type,length(controls)+1)
  }

####try with a single AMME function for each model type
  if(algorithm=="parametric"){
    if(class(macro_model)[1]=="Gam"){
      stop("Parametric estimation does not support generalized additive models. Use nonparametric estimation instead.")
    }

   results<-AMME_param(micro_model=micro_model,
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
                       SAOM_var=SAOM_var,
                       time_interval=time_interval,
                       covar_list=covar_list,
                       edgelist=edgelist,
                       group_id=group_id,
                       node_numbers=node_numbers)

  }else{

    results<-AMME_nonparam(micro_model=micro_model,
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
                        covar_list=covar_list,
                        edgelist=edgelist,
                        net_logit_y=net_logit_y,
                        net_logit_x=net_logit_x,
                        group_id=group_id,
                        node_numbers=node_numbers)

  }

  return(results)



}
