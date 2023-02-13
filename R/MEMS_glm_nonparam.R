#
# This function estimates the MEMS using GLM or GLMER with parametric algorithm
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
#group_id is an optional vector of group identifiers to use when estimating a glm or glmer on grouped data (i.e., multiple time periods, multiple networks) on a dependent tie variable. When specified, the MEMS command will induce unique networks for each grouping factor. If left unspecified, all groups/time period are pooled. If using glmer, the grouping factor does not have to be provided as part of the model or used as a random effect.
#node numbers is a vector listing the number of nodes in each network when using GLM or GLMER. If estimating MEMS aggregated over all networks, this shoud be the total number of nodes in all networks. Required when using GLM or GLMER, ignored otherwise.


MEMS_glm_nonparam <- function(model=model,
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
                              control_functions=control_functions) {

  interval<-sort(interval)

  link_fun<-family(model)
  if(!link_fun$link%in%c("logit","probit","identity")){
    stop("netmediate currently only supports logit, probit, and linear GLMs.")
  }

  if(link_fun$link=="identity"){
    warning("Linear probability models not recommended as predictions may be nonsensical. Proceeding anyway.")
  }





  message("Computing MEMS over ",interval[1],"-",interval[length(interval)]," interval.")

  #set up model
  if(class(model)[1]%in%c("glmerMod","lmerMod")){
    coef<-lme4::fixef(model)

  }else{
    coef<-model$coefficients
  }

  model_mat<-stats::model.frame(model)
  model_form<-stats::formula(model)
  if(nrow(model_mat)>200000){
    message("More than 200,000 observations in the model matrix. Estimation may take awhile. Parametric estimation is generally less computationally expensive.")
  }

  theta<-matrix(NA,nrow=nsim,ncol=length(coef))
  colnames(theta)<-names(coef)

  #start bootstrap
  message("Getting bootstrap parameter vector.")

  for(j in 1:nsim){

    bs_dat<-model_mat[sample(nrow(model_mat),replace=TRUE),]

    if(class(model)[1]%in%c("glmerMod","lmerMod")){
        bs_model<-lme4::glmer(model_form,
                              data=bs_dat,
                              family=link_fun$family)
        theta[j,]<-lme4::fixef(bs_model)

    }else{
      bs_model<-glm(model_form,
                    data=bs_dat,
                    family=link_fun$family)
      theta[j,]<-coef(bs_model)
    }


    if(silent==FALSE){

      message("Getting bootstrap parameters ", j, " of ", nsim, " complete.")

    }

  }




  if(class(model)[1]%in%c("glmerMod","lmerMod")){
    model_mat<-model_mat[,!colnames(model_mat)%in%names(lme4::ranef(model))]
  }

  if(is.null(group_id)){
    group_id<-rep(1,nrow(model_mat))
  }
  if(length(group_id)!=nrow(model_mat)){
    stop("Number of observations in group_id does not match number of observations in model matrix.")
  }

  unique_groups<-unique(group_id)
  mat_list<-list()
  output_list<-list()
  #create list of values to predict over
  for(i in 1:length(interval)){
    mat_list[[i]]<-model_mat

  }



  #create lists of data for each link_id and control. These lists provide the data for AMME function
  if(!is.null(link_id)){

    link_list_data<-list()
    for(i in 1:nsim){
      link_list_data[[i]]<-matrix(NA,nrow=length(unique(link_id[[1]])),ncol=length(interval))
      rownames(link_list_data[[i]])<-unique(link_id[[1]])
    }

    if(length(controls)>0){
      controls_list<-as.list(controls)
      for(i in 1:length(controls)){
        sim_list<-list()
        for(j in 1:nsim){
          sim_list[[j]]<-matrix(NA,nrow=length(unique(link_id[[i+1]])),ncol=length(interval))
          rownames(sim_list[[j]])<-unique(link_id[[i+1]])
        }
        controls_list[[i]]<-sim_list
        names(controls_list)<-controls

      }
    }else{
      controls_list<-NULL
    }
  }



  for(i in 1:length(unique_groups)){
    output_list[[i]]<-matrix(NA,nrow=nsim,ncol=length(interval))
  }
  aMEMS_tracer<-0
  if(length(group_id)>1){aMEMS_tracer<-1}
  model_mat[,1]<-1 ##assign constant value for intercept to use in predictions


  message("Computing MEMS over ",interval[1],"-",interval[length(interval)]," interval")
  #simulate values for j parameters
  for(j in 1:nrow(theta)){

    #create networks for k intervals

    net_list<-list()
    for(i in 1:length(mat_list)){

      pred_mat<-mat_list[[i]]
      cbcoef<-theta[j,]
      cbcoef[micro_process]<-cbcoef[micro_process]*interval[i]

      #predict ties
      lp <- as.matrix(pred_mat)%*%cbcoef

      if(link_fun$link=="logit"){
        result <- 1/(1 + exp(-lp))
      }

      if(link_fun$link=="probit"){
        result<-VGAM::probitlink(lp,inverse=TRUE)
      }

      if(link_fun$link=="identity"){
        result<-lp
        result[result>1]<-1
        result[result<0]<-0

      }

      net_list[[i]]<-result
    }
    ##now need to calculate statistics on each network in list
    #then need to calculate difference statistics and MEMS or aMEMS
    #then need to store output
    for(entry in 1:length(unique_groups)){

      if(object_type[1]%in%c("network")){

        for(i in 1:length(net_list)){
          entry_network_interval <- network::as.network(sna::rgraph(n=node_numbers[entry],
                                                                    tprob=net_list[[i]][group_id==unique_groups[entry]]))
          b<-macro_function(entry_network_interval)

          if(!is.null(link_id)){
            node_start<- which(is.na(link_list_data[[j]][,i]))[1]
            node_index<-node_start+(length(b)-1)
            link_list_data[[j]][node_start:node_index,i]<-b
          } ##store vector of output when calling AMME

          if(length(b)>1){
            b<-mean(b,na.rm=TRUE)
            aMEMS_tracer<-1} ##handle node and multiple obs characteristics for MEMS output
          output_list[[entry]][j,i]<-b
        }

      }else{
        #for igraph functions
        for(i in 1:length(net_list)){
          entry_network_interval <- network::as.network(sna::rgraph(n=node_numbers[entry],
                                                                    tprob=net_list[[i]][group_id==unique_groups[entry]]))
          entry_network_interval<-intergraph::asIgraph(entry_network_interval)
          b<-macro_function(entry_network_interval)

          if(!is.null(link_id)){
            node_start<- which(is.na(link_list_data[[j]][,i]))[1]
            node_index<-node_start+(length(b)-1)
            link_list_data[[j]][node_start:node_index,i]<-b
          } #store vector of output when calling AMME
          if(length(b)>1){
            b<-mean(b,na.rm=TRUE)
            aMEMS_tracer<-1} ##handle node and multiple obs characteristics
          output_list[[entry]][j,i]<-b

        }

      }


      ##get controls

      if(length(controls)>0){

        for(control in 1:length(controls_list)){
          if(object_type[control+1]%in%c("network")){

            for(i in 1:length(net_list)){
              entry_network_interval <- network::as.network(sna::rgraph(n=node_numbers[entry],
                                                                        tprob=net_list[[i]][group_id==unique_groups[entry]]))
              b<-control_functions[[control]](entry_network_interval)
              node_start<- which(is.na(controls_list[[control]][[j]][,i]))[1]
              node_index<-node_start+(length(b)-1)
              controls_list[[control]][[j]][node_start:node_index,i]<-b
            }

          }else{
            #for igraph functions
            for(i in 1:length(net_list)){
              entry_network_interval <- network::as.network(sna::rgraph(n=node_numbers[entry],
                                                                        tprob=net_list[[i]][group_id==unique_groups[entry]]))
              entry_network_interval<-intergraph::asIgraph(entry_network_interval)
              b<-control_functions[[control]](entry_network_interval)
              node_start<- which(is.na(controls_list[[control]][[j]][,i]))[1]
              node_index<-node_start+(length(b)-1)
              controls_list[[control]][[j]][node_start:node_index,i]<-b

            }

          }

        }
      } #closes if controls statement




    }#close entry loop


    if(silent==FALSE){
      print(paste("Simulation ",j," out of", nsim," complete"))
    }
  } #close j loop

  if(aMEMS_tracer==1 &is.null(link_id)){
    message("More than one macro statistic is being calculated. Reporting the aMEMS.")
  }

  for(i in 1:length(output_list)){

    if(i == 1){
      output_data<-output_list[[i]]
    }else{
      output_data<-rbind(output_data,output_list[[i]])
    }

  }


  diff_data<-matrix(NA,nrow=nrow(output_data),ncol=ncol(output_data)-1)

  for(i in 1:ncol(diff_data)){
    k<-i+1

    diff_data[,i]<-output_data[,k]-output_data[,i]
  }

  summary_dat<-matrix(NA,nrow=2,ncol=5)
  rownames(summary_dat)<-c("(a)MEMS","Prop. Change in M")
  colnames(summary_dat)<-c("Estimate","Std. Dev.","lower 95% CI","Upper 95% CI","MC p-val")

  summary_dat[1,1]<-mean(diff_data,na.rm=TRUE)
  summary_dat[1,2]<-sd(diff_data,na.rm=TRUE)

  if(summary_dat[1,1]<0){
    summary_dat[1,5]<-length(diff_data[which(diff_data>=0)])/(nsim*length(unique_groups)*ncol(diff_data))
  }else{
    summary_dat[1,5]<-length(diff_data[which(diff_data<=0)])/(nsim*length(unique_groups)*ncol(diff_data))

  }
  summary_dat[1,3]<-quantile(diff_data,.025,na.rm=TRUE)
  summary_dat[1,4]<-quantile(diff_data,.975,na.rm=TRUE)



  prop_change<-matrix(NA,nrow=nrow(output_data),ncol=ncol(output_data)-1)

  for(i in 1:ncol(prop_change)){
    k<-i+1

    prop_change[,i]<-(diff_data[,i]/output_data[,k])
  }

  prop_change<-prop_change[!is.infinite(prop_change)]

  summary_dat[2,1]<-mean(prop_change,na.rm=TRUE)

  if(full_output==FALSE){
    return(summary_dat)
  }else{
    out_dat<-list(summary_dat=summary_dat,
                  output_data=output_data,
                  mems_samples=diff_data)
    if(is.null(link_id)){

      return(out_dat) #return data no controls

    }else{

      #return data for AMME call
      out_dat<-list(out_dat_main=link_list_data,
                    out_dat_controls=controls_list)

      return(out_dat)


    }
  }
}
