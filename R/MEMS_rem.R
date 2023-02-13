#
# This function calls the appropriate REM algorithm for estimating the MEMS
#
#
# model is the REM object
#micro_process is the micro process of interest, provided as a character string
#macro_function is the function to be calculated on the network object. NOTE: currently only supports statistics calculated on network and igraph objects.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#algorithm tells netmediate which estimation routine to use.
#full_output tells R whether to return simulated distribution in addition to summary statistic
#time interval is an optional parameter to be used with rem.dyad models or GLM/GLMER models for REM data. Time interval tells which time points to compute the MEMS. Two options are available: aggregate or a continuous vector of numeric values. If "aggregate" is specified, a single MEMS is calculated by creating an aggregate cross-sectional representation of the event sequence treating the final observed time of event activity as the limiting condition for simulated event sequences. If a continuous interval is specified of the form c(0,2, 3), unique networks are created containing all ties that occur between those intervals. For example, if c(1,2,5) is specified, unique networks are constructed using all ties that occur between the 1 and 2 time periods and the 2 and 5 time periods, and the MEMS is computed over those two networks. If unspecified, defaults to all time periods.
#edgelist is an optional parameter only used with rem.dyad models when the time interval is left unspecified or set to "aggregate". In this case, the edgelist should be a 3 column matrix used as the dependent variable for rem.dyad. The first column contains the time variable, the second column the senders, and the third column the receivers.

MEMS_rem<- function(model,
                     micro_process,
                     macro_function,
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
                     control_functions=control_functions) {

message("Global and node-level statistics that are a function of *exogenous* node-level and edge-level characteristics are not intrinsically available for rem.dyad type objects. If the macro statistic of interest is a function of node or edge-level exogenous characteristics, adjust macro_function to assign exogenous node and edge attributes before computing macro statistic of interest. See help file for examples.")


  if(algorithm=="parametric"){
    results<-MEMS_rem_param(model=model,
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
                            mediator=mediator,
                            link_id=link_id,
                            controls=controls,
                            control_functions=control_functions)

  }else{

   stop("Nonparametric estimation using rem.dyad object is not currently implemented due to computational demands. To use nonparametric estimation with REM, try respecifying the model as a GLM.")


  return(results)

  }
}
