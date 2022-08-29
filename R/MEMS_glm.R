#
# This function calls the appropriate GLM algorithm for estimating the MEMS
#
#
# model is the GLM object
#micro_process is the micro process of interest, provided as a character string
#macro_function is the function to be calculated on the network object. NOTE: currently only supports statistics calculated on network and igraph objects.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#full_output tells R whether to return simulated distribution in addition to summary statistic
#algorithm tells netmediate which estimation routine to use.
#group_id is an optional vector of group identifiers to use when estimating a glm or glmer on grouped data (i.e., multiple time periods, multiple networks) on a dependent tie variable. When specified, the MEMS command will induce unique networks for each grouping factor. If left unspecified, all groups/time period are pooled. If using glmer, the grouping factor does not have to be provided as part of the model or used as a random effect.
#node numbers is a vector listing the number of nodes in each network when using GLM or GLMER. If estimating MEMS aggregated over all networks, this shoud be the total number of nodes in all networks. Required when using GLM or GLMER, ignored otherwise.


MEMS_glm<- function(model=model,
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
                    control_functions=control_functions) {

  message("Global and node-level statistics that are a function of *exogenous* node-level and edge-level characteristics are not intrinsically available for glm and glmer objects. If the macro statistic of interest is a function of node or edge-level exogenous characteristics, adjust macro_function to assign exogenous node and edge attributes before computing macro statistic of interest. See help file for examples.")
  if(is.null(node_numbers)){
    stop("node_numbers must be provided when using GLM or GLMER.")
  }

  if(algorithm=="parametric"){
    results<-MEMS_glm_param(model=model,
                            micro_process=micro_process,
                            macro_function=macro_function,
                            object_type=object_type,
                            interval=interval,
                            nsim=nsim,
                            silent=silent,
                            full_output=full_output,
                            group_id=group_id,
                            node_numbers=node_numbers,
                            mediator=mediator,
                            link_id=link_id,
                            controls=controls,
                            control_functions=control_functions)

  }else{


    results<-MEMS_glm_nonparam(model=model,
                               micro_process=micro_process,
                               macro_function=macro_function,
                               object_type=object_type,
                               interval=interval,
                               nsim=nsim,
                               silent=silent,
                               full_output=full_output,
                               group_id=group_id,
                               node_numbers=node_numbers,
                               mediator=mediator,
                               link_id=link_id,
                               controls=controls,
                               control_functions=control_functions)

  }

  return(results)

}
