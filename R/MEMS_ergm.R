#
# This function calls the appropriate ERGM algorithm for estimating the MEMS
#
#
# model is the ERGM object
#micro_process is the micro process of interest, provided as a character string
#macro_function is the function to be calculated on the network object. NOTE: currently only supports statistics calculated on network and igraph objects.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#full_output tells R whether to return simulated distribution in addition to summary statistic
#algorithm tells netmediate which estimation routine to use.


MEMS_ergm<- function(model,
                            micro_process,
                            macro_function,
                            object_type=object_type,
                            interval=interval,
                            nsim=nsim,
                            silent=silent,
                            full_output=full_output,
                            algorithm=algorithm,
                            mediator=mediator,
                            link_id=link_id,
                            controls=controls,
                            control_functions=control_functions) {



      if(algorithm=="parametric"){
          if(class(model)[1]=="btergm"){
            stop("MEMS using bootstrap TERGM must be estimated nonparametrically. Respecify the estimation algorithm.")
          }
        results<-MEMS_ergm_param(model,
                                 micro_process,
                                 macro_function,
                                 object_type=object_type,
                                 interval=interval,
                                 nsim=nsim,
                                 silent=silent,
                                 full_output=full_output,
                                 mediator=mediator,
                                 link_id=link_id,
                                 controls=controls,
                                 control_functions=control_functions)

      }else{

        if(class(model)[1]!="btergm"){

            if(is.curved(model)){
             warning("Bootstrap pseudolikelihood is currently the only nonparametric method implemented. This is only consistent in large networks. If the network is small or moderately sized, consider re-estimating parametrically.")
            }
          }

        results<-MEMS_ergm_nonparam(model,
                                 micro_process,
                                 macro_function,
                                 object_type=object_type,
                                 interval=interval,
                                 nsim=nsim,
                                 silent=silent,
                                 full_output=full_output,
                                 mediator=mediator,
                                 link_id=link_id,
                                 controls=controls,
                                 control_functions=control_functions)

      }

  return(results)

}
