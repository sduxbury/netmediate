#
# This function estimates the MEMS using ERGM with parametric algorithm
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



MEMS_ergm_param <- function(model,
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
                            control_functions=control_functions) {

  coef<-model$coefficients
  cov_mat<-model$covar

  if(any(is.na(coef))|
     any(is.infinite(coef))){
    stop("Infinite or missing values in parameter estimates. Algorithm cannot continue.")

  }

  if(any(is.na(cov_mat))|
     any(is.infinite(cov_mat))){
    stop("Infinite or missing values in covariance matrix estimates. Algorithm cannot continue.")

  }
  interval<-sort(interval) #order from lowest to highest


  theta<-MASS::mvrnorm(n=nsim,
                       mu=coef,
                       Sigma=cov_mat,
                       empirical = TRUE)

  ergm_mat<-ergMargins::edge.prob2(model)
  mat_list<-list()

  #create list of values to predict over
  for(i in 1:length(interval)){
    mat_list[[i]]<-ergm_mat

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

  message("Computing MEMS over ",interval[1],"-",interval[length(interval)]," interval")
  output_data<-matrix(NA,nrow=nsim,ncol=length(interval))
  aMEMS_tracer<-0

  #simulate values for j parameters
  for(j in 1:nrow(theta)){

    #create networks for k intervals

    net_list<-list()
    for(i in 1:length(mat_list)){

      pred_mat<-mat_list[[i]]
      start.drops<-ncol(pred_mat)-5
      pred_mat<-pred_mat[,-c(1,start.drops:ncol(pred_mat))]
      cbcoef<-theta[j,]
      cbcoef[micro_process]<-cbcoef[micro_process]*interval[i]

      #predict ties
      lp <- as.matrix(pred_mat)%*%cbcoef
      result <- c(1/(1 + exp(-lp)))
      pred_mat<-mat_list[[i]]
      pred_mat$y <- rbinom(nrow(mat_list[[i]]),1,result)      # create predicted ties

      #create network
      el<-pred_mat[,c("i","j","y")]
      el<-el[el$y==1,-c(3)]
      el<-as.matrix(el)

      #create network, assign vertex attributes, by creating empty network and adding new edges
      net_list[[i]]<-model$network
      net_list[[i]][,]<-0
      net_list[[i]]<-network::add.edges(net_list[[i]],tail=el[,1],head=el[,2])

    }

  ##now need to calculate statistics on each network in list
    #then need to calculate difference statistics and MEMS or aMEMS
    #then need to store output

  if(object_type[1]%in%c("network")){

    for(i in 1:length(net_list)){
      b<-macro_function(net_list[[i]])
      if(!is.null(link_id)){
        link_list_data[[j]][,i]<-b
      } ##store vector of output when calling AMME

      if(length(b)>1){
        b<-mean(b,na.rm=TRUE)
        aMEMS_tracer<-1} ##handle node and multiple obs characteristics for MEMS output
      output_data[j,i]<-b
    }

  }else{
            #for igraph functions
      for(i in 1:length(net_list)){
        net_list[[i]]<-intergraph::asIgraph(net_list[[i]])
        b<-macro_function(net_list[[i]])
        if(!is.null(link_id)){
          link_list_data[[j]][,i]<-b
        } #store vector of output when calling AMME
        if(length(b)>1){
          b<-mean(b,na.rm=TRUE)
          aMEMS_tracer<-1} ##handle node and multiple obs characteristics
        output_data[j,i]<-b
        if(length(controls)>0){
          net_list[[i]]<-intergraph::asNetwork(net_list[[i]]) # convert back to network for control operations
        }

      }

    }

    ##get controls

    if(length(controls)>0){

      for(control in 1:length(controls_list)){

       if(object_type[control+1]%in%c("network")){

        for(i in 1:length(net_list)){
          b<-control_functions[[control]](net_list[[i]])
          controls_list[[control]][[j]][,i]<-b
        }

      }else{
        #for igraph functions
        for(i in 1:length(net_list)){
          net_list[[i]]<-intergraph::asIgraph(net_list[[i]])
          b<-control_functions[[control]](net_list[[i]])
          controls_list[[control]][[j]][,i]<-b
          net_list[[i]]<-intergraph::asNetwork(net_list[[i]]) #convert back to network object for further loops

          }

        }

      }
    } #closes if controls statement

    if(silent==FALSE){
      print(paste("Simulation ",j," out of", nsim," complete"))
    }
  }

  if(aMEMS_tracer==1 & is.null(link_id)){
    message("More than one macro statistic is being calculated. Reporting the aMEMS.")
  }







  ###package the output
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
    summary_dat[1,5]<-length(diff_data[which(diff_data>=0)])/(nsim*ncol(diff_data))
  }else{
    summary_dat[1,5]<-length(diff_data[which(diff_data<=0)])/(nsim*ncol(diff_data))

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
