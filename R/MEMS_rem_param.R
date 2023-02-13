#
# This function estimates the MEMS using REM with parametric algorithm
#
#
# model is the REM object
#micro_process is the micro process of interest, provided as a character string
#macro_function is the function to be calculated on the network object. NOTE: currently only supports statistics calculated on network and igraph objects.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#time interval is an optional parameter to be used with rem.dyad models or GLM/GLMER models for REM data. Time interval tells which time points to compute the MEMS. Two options are available: aggregate or a continuous vector of numeric values. If "aggregate" is specified, a single MEMS is calculated by creating an aggregate cross-sectional representation of the event sequence treating the final observed time of event activity as the limiting condition for simulated event sequences. If a continuous interval is specified of the form c(0,2, 3), unique networks are created containing all ties that occur between those intervals. For example, if c(1,2,5) is specified, unique networks are constructed using all ties that occur between the 1 and 2 time periods and the 2 and 5 time periods, and the MEMS is computed over those two networks. If unspecified, defaults to all time periods.
#full_output tells R whether to return simulated distribution in addition to summary statistic
#edgelist is an optional parameter only used with rem.dyad models when the time interval is left unspecified or set to "aggregate". In this case, the edgelist should be a 3 column matrix used as the dependent variable for rem.dyad. The first column contains the time variable, the second column the senders, and the third column the receivers.


MEMS_rem_param <- function(model=model,
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
                           control_functions=control_functions) {

  if(is.null(time_interval)){
    time_interval<-"aggregate"
  }

  if(is.null(edgelist) &"aggregate"%in%time_interval){
    stop("An edgelist must be provided if time interval is left null or aggregating over all events.")
  }
  time_interval<-unique(time_interval)

  if(length(time_interval)>20){
    message("More than 10 time intervals specified. Long run times are possible.")
  }

  aMEMS_tracer<-0
  coef<-model$coef
  cov_mat<-model$cov

  if(any(is.na(coef))|
     any(is.infinite(coef))){
    stop("Infinite or missing values in parameter estimates. Algorithm cannot continue.")

  }

  if(any(is.na(cov_mat))|
     any(is.infinite(cov_mat))){
    stop("Infinite or missing values in covariance matrix estimates. Algorithm cannot continue.")

  }
  theta<-MASS::mvrnorm(n=nsim,
                       mu=coef,
                       Sigma=cov_mat,
                       empirical = TRUE)



  interval<-sort(interval) #order from lowest to highest

  ##create time interval windows
  if(length(time_interval)==1){
    if(time_interval=="aggregate"){
      time_interval<-c(0,max(edgelist[,1]))
    }else{
      time_interval<-c(0,time_interval)
    }
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



  interval_lengths<-length(time_interval)-1
  output_list<-list()
  for(i in 1:interval_lengths){
    output_list[[i]]<-matrix(NA,nrow=nsim,ncol=length(interval))

  }

  ##start simulation--loop over each time interval
  message("Computing MEMS over ",interval[1],"-",interval[length(interval)]," interval")

  for(entry in 1:length(output_list)){

     future_interval<-time_interval[entry+1]

    #simulate values for j parameters
    for(j in 1:nrow(theta)){

      #create networks for k intervals

      net_list<-list()
      for(i in 1:length(interval)){

        theta2<-theta[j,]
        theta2[micro_process]<-theta2[micro_process]*interval[i]
        if(is.null(covar_list)){ #include covariates if provided
          sim_sequence<-simulate(model,coef=theta2)
        }else{
          sim_sequence<-simulate(model,coef=theta2,covar=covar_list)

        }

          sim_sequence<-sim_sequence[which(sim_sequence[,1]<=time_interval[time_interval==future_interval] &
                                             sim_sequence[,1]>=time_interval[time_interval==time_interval[entry]]),]


        #create network
        el<-as.matrix(sim_sequence[,-c(1)])
        net_list[[i]]<- as.network(matrix(0, nrow=model$n, ncol=model$n))
        net_list[[i]]<-network::add.edges(net_list[[i]],tail=el[,1],head=el[,2])

      }

      ##now need to calculate statistics on each network in list
      #then need to calculate difference statistics and MEMS or aMEMS
      #then need to store output

      if(object_type[1]%in%c("network")){

        for(i in 1:length(interval)){
          b<-macro_function(net_list[[i]])
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

      }else{
        #for igraph functions
        for(i in 1:length(interval)){
          net_list[[i]]<-intergraph::asIgraph(net_list[[i]])
          b<-macro_function(net_list[[i]])
          if(!is.null(link_id)){
            node_start<- which(is.na(link_list_data[[j]][,i]))[1]
            node_index<-node_start+(length(b)-1)
            link_list_data[[j]][node_start:node_index,i]<-b
          } #store vector of output when calling AMME
          if(length(b)>1){
            b<-mean(b,na.rm=TRUE)
            aMEMS_tracer<-1} ##handle node and multiple obs characteristics
          output_list[[entry]][j,i]<-b
          if(length(controls)>0){
            net_list[[i]]<-intergraph::asNetwork(net_list[[i]]) # convert back to network for control operations
          }

        }

      }

      ##get controls

      if(length(controls)>0){

        for(control in 1:length(controls_list)){
          if(object_type[control+1]%in%c("network")){

            for(i in 1:length(interval)){
              b<-control_functions[[control]](net_list[[i]])
              node_start<- which(is.na(controls_list[[control]][[j]][,i]))[1]
              node_index<-node_start+(length(b)-1)
              controls_list[[control]][[j]][node_start:node_index,i]<-b
            }

          }else{
            #for igraph functions
            for(i in 1:length(interval)){
              net_list[[i]]<-intergraph::asIgraph(net_list[[i]])
              b<-control_functions[[control]](net_list[[i]])
              node_start<- which(is.na(controls_list[[control]][[j]][,i]))[1]
              node_index<-node_start+(length(b)-1)
              controls_list[[control]][[j]][node_start:node_index,i]<-b
              net_list[[i]]<-intergraph::asNetwork(net_list[[i]]) #convert back to network object for further loops

            }

          }

        }
      } #closes if controls statement




      if(silent==FALSE){
        print(paste("Simulation ",j," out of", nsim," complete for time_interval ", entry))
      }
    }

  }

  if(aMEMS_tracer==1 & is.null(link_id)){
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
  rownames(summary_dat)<-c("aMEMS","Prop. Change in M")
  colnames(summary_dat)<-c("Estimate","Std. Dev.","lower 95% CI","Upper 95% CI","MC p-val")

  summary_dat[1,1]<-mean(diff_data,na.rm=TRUE)
  summary_dat[1,2]<-sd(diff_data,na.rm=TRUE)

  if(summary_dat[1,1]<0){
    summary_dat[1,5]<-length(diff_data[which(diff_data>=0)])/(interval_lengths*nsim*ncol(diff_data))
  }else{
    summary_dat[1,5]<-length(diff_data[which(diff_data<=0)])/(interval_lengths*nsim*ncol(diff_data))

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

  message("Note that a single aMEMS is returned for all networks. To obtain distinct MEMS estimates for each network, call the MEMS function specifying a single model, rather than a list of models.")
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
  }}
