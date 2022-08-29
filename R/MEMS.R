#'
#' This function calls the appropriate algorithm for estimating the MEMS. The dependent variable is assumed to be a dyad or dyad equivalent (e.g., dyad-time-period, dyad-group)
#'
#'
#'@model is the model object or in the case of multiple models fit to multiple networks, a list of model objects.
#'@micro_process is the micro process of interest, provided as a character string that matches the model output.
#'@macro_function is the function to be calculated on the network object. The function can be provided as the name of a function in the sna package to implement default settings. To include specific parameters within that function, write a user function of the form macro_function<-function(x){gtrans(x,measure="weakcensus)}. To call functions from igraph, igraph must be loaded in the global environment or called within the user function and object_type must be specified as object_type=="igraph". MEMS currently only implements functions for igraph and network objects, but that behavior can be overridden by specifying a user function that takes an igraph or network object as its input and computes relevant statistics from there. For example, the following user function calculates statistics from the netseg package for R assuming that x is an igraph object: macro_function=function(x){netseg::assort(x,"type")}. Also note that MEMS only takes into account the attributes included in the model. To include node or edge attributes not specified in the model, write a user function that takes a network or igraph object as its input and assigns the attribute. For example, if we want to calculate race assortativity using igraph, but race is not included in the model, we would specify network_object="igraph" and the user function: macro_function=function(x,race){require(igraph); V(x)$race<-race;assortativity(x,V(x)$race)}. Whatever macro function is specified, the output should either be a numeric value or a vector of numeric values for the algorithm to proceed.
#'@object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#'@the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#'@nsim is the number of Monte Carlo samples
#'@silent tells R whether to provide updates on the algorithm's progress
#'@algorithm tells netmediate which estimation routine to use.
#'@full_output tells R whether to return simulated distribution in addition to summary statistic
#'@SAOM_data is a siena object that contains the data to use for SAOM estimation. Only required when using SAOM to calculate the MEMS.
#'@SAOM_var is an optional list of varCovar and varDyadCovar to calculate on each simulated network. Only required when the macro statistic is a funciton of time varying node or dyad attributes. CoCovar and coDyadCovar (time invariant node covariates and dyad covariates) objects are handled internally. List must be provided with named entries if used in a call to AMME with a SAOM micro-model. If providing a list of sienaFit objects, SAOM_var should be a list of lists the correspond to each entry sienaFit object.
#'@time_interval is an optional parameter to be used with rem.dyad models. Time interval tells which time points to compute the MEMS. Two options are available: aggregate or a continuous vector of numeric values. If "aggregate" is specified, a single MEMS is calculated by creating an aggregate cross-sectional representation of the event sequence treating the final observed time of event activity as the limiting condition for simulated event sequences. If a continuous interval is specified of the form c(0,2, 3), unique networks are created containing all ties that occur between those intervals. For example, if c(1,2,5) is specified, unique networks are constructed using all ties that occur between the 1 and 2 time periods and the 2 and 5 time periods, and the MEMS is computed over those two networks. If unspecified, defaults to all time periods. The same behavior may be obtained when using a GLM or GLMER to estimate a relational event model by assigning group memberships that align with the desired time intervals and then passing these memberships with the group_id parameter.
#'@covar_list is a list of covariates used to estimate the MEMS when using a rem.dyad object.
#'@edgelist is an edgelist that acts as the dependent variable in rem.dyad when using a rem.dyad object. Otherwise left empty.
#'@net_logit_y and net_logit_x are the dependent and independent variables specified using network logistic regression with the netlogit function in sna. y must be a binary adjacency matrix and x stack of independent network variables.
#'@group_id is an optional vector of group identifiers to use when estimating a glm or glmer on grouped data (i.e., multiple time periods, multiple networks) on a dependent tie variable. When specified, the MEMS command will induce unique networks for each grouping factor. If left unspecified, all groups/time period are pooled. If using glmer, the grouping factor does not have to be provided as part of the model or used as a random effect.
#'@node numbers is a vector listing the number of nodes in each network when using GLM or GLMER. If estimating MEMS aggregated over all networks, this shoud be the total number of nodes in all networks. Required when using GLM or GLMER, ignored otherwise.
#'@mediator, macro_function,link_id,controls, and control_functions are optional parameters intended for internal use with the AMME function. They are ignored when specified externally

MEMS<- function(model,
                     micro_process,
                     macro_function,
                     object_type=NULL,
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
                     node_numbers=NULL,
                     mediator=NULL,
                     link_id=NULL,
                     controls=NULL,
                     control_functions=NULL) {

###checking compatability
  if(!algorithm%in%c("parametric","nonparametric")){
    stop("Must specify either 'parametric' or 'nonparametric' estimation algorithm.")
  }

  if(length(algorithm)>1){
    stop("Must specify a single algorithm")
  }

  if(!is.null(object_type)&&
     !any(object_type%in%c("igraph","network"))){
    stop("Macro functions currently only applied using igraph and network objects. Adjust the function to comport with network or igraph objects.")
  }

  if(is.null(object_type)){
    object_type<-"network"
  }






####call correct function
  if(class(model)[1]=="list"){

    type_list<-list()
    for(i in 1:length(model)){
      type_list[[i]]<-class(model[[i]])
      if(!type_list[[i]]%in%c("ergm","sienaFit")){
        stop("Pooled MEMS estimation currently only implemented for ERGM and SAOM.")
      }

      if(i>1){
        if(type_list[[i]]!=type_list[[i-1]]){
          stop("All micro models must be the same type of model.")
        }
      }

    }


    if(type_list[[1]]=="ergm"){
      results<-MEMS_pooled_ergm(model=model,
                         micro_process=micro_process,
                         macro_function=macro_function,
                         object_type=object_type,
                         interval=interval,
                         nsim=nsim,
                         silent=silent,
                         full_output=full_output,
                         algorithm=algorithm,
                         mediator=mediator,
                         link_id=link_id,
                         controls=controls,
                         control_functions=control_functions)

    }

    if(type_list[[1]]=="sienaFit"){
      if(algorithm=="nonparametric"){
        stop("Nonparametric estimation not supported for stochastic actor-oriented models. Bootstrapping violates core model assumptions. Please re-specify parametric estimation.")
      }

      if(is.null(SAOM_data)){
        stop("SAOM MEMS calculation requires SAOM data object. Specify SAOM data source and try again.")
      }
      if(class(SAOM_data)[1]!="list"){
        stop("SAOM data must be provided as a list of siena objects for pooled estimatino to proceed.")
      }
      results<-MEMS_pooled_saom(model=model,
                         micro_process=micro_process,
                         macro_function=macro_function,
                         object_type=object_type,
                         interval=interval,
                         nsim=nsim,
                         silent=silent,
                         SAOM_data=SAOM_data,
                         SAOM_var=SAOM_var,
                         full_output=full_output,
                         algorithm=algorithm,
                         mediator=mediator,
                         link_id=link_id,
                         controls=controls,
                         control_functions=control_functions)

    }
    return(results)

  }else{


  if(!class(model)[1]%in%c("btergm","ergm","sienaFit","rem.dyad","glm","lm","glmerMod","lmerMod","netlogit")){
    stop("Model type currently not supported. Netmediate currently only supports btergm, ergm, sienaFit, rem.dyad, glm, lm, glmerMod, lmerMod, and netlogit type objects.")
  }

  if(class(model)[1]%in%c("btergm","ergm")){
    results<-MEMS_ergm(model=model,
                       micro_process=micro_process,
                       macro_function=macro_function,
                       object_type=object_type,
                       interval=interval,
                       nsim=nsim,
                       silent=silent,
                       full_output=full_output,
                       algorithm=algorithm,
                       mediator=mediator,
                       link_id=link_id,
                       controls=controls,
                       control_functions=control_functions)
  }



  if(class(model)[1]%in%"sienaFit"){

    if(algorithm=="nonparametric"){
      stop("Nonparametric estimation not supported for stochastic actor-oriented models. Bootstrapping violates core model assumptions. Please re-specify parametric estimation.")
    }
    if(is.null(SAOM_data)){
      stop("SAOM MEMS calculation requires SAOM data object. Specify SAOM data source and try again.")
    }

    results<-MEMS_saom(model=model,
                       micro_process=micro_process,
                       macro_function=macro_function,
                       object_type=object_type,
                       interval=interval,
                       nsim=nsim,
                       silent=silent,
                       SAOM_data=SAOM_data,
                       SAOM_var=SAOM_var,
                       full_output=full_output,
                       algorithm=algorithm,
                       mediator=mediator,
                       link_id=link_id,
                       controls=controls,
                       control_functions=control_functions)
  }


  if(class(model)[1]%in%c("rem.dyad")){
    results<-MEMS_rem(model=model,
                       micro_process=micro_process,
                       macro_function=macro_function,
                       object_type=object_type,
                       interval=interval,
                       nsim=nsim,
                       silent=silent,
                       full_output=full_output,
                      time_interval=time_interval,
                      covar_list=covar_list,
                       edgelist=edgelist,
                       algorithm=algorithm,
                      mediator=mediator,
                      link_id=link_id,
                      controls=controls,
                      control_functions=control_functions)
  }


  if(class(model)[1]%in%c("glm","lm","glmerMod","lmerMod")){
    results<-MEMS_glm(model=model,
                       micro_process=micro_process,
                       macro_function=macro_function,
                       object_type=object_type,
                       interval=interval,
                       nsim=nsim,
                       silent=silent,
                       full_output=full_output,
                       algorithm=algorithm,
                       group_id=group_id,
                       node_numbers=node_numbers,
                      mediator=mediator,
                      link_id=link_id,
                      controls=controls,
                      control_functions=control_functions)
  }




  if(class(model)[1]%in%c("netlogit")){
    if(algorithm=="parametric"){
       message("Parametric estimation not possible with QAP. Returning nonparametric estimates instead.")
    }
    if(is.null(net_logit_x)|is.null(net_logit_y)){
      stop("net_logit_x and net_logit_y must be specified when estimating MEMS with network logistic regression.")
    }
    results<-MEMS_QAP(model=model,
                       micro_process=micro_process,
                       macro_function=macro_function,
                       object_type=object_type,
                       interval=interval,
                       nsim=nsim,
                       silent=silent,
                       net_logit_y=net_logit_y,
                       net_logit_x=net_logit_x,
                       full_output=full_output,
                      mediator=mediator,
                      link_id=link_id,
                      controls=controls,
                      control_functions=control_functions)
  }



  return(results)
}
}
