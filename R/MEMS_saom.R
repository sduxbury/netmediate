#
# This function estimates the MEMS using ERGM with parametric algorithm
#
#
# model is the ERGM object
#micro_process is the micro process of interest, provided as a character string. For SAOM, it should be as per the effectName parameter in your SAOM effects object. For example, effects_obj$effectName
#macro_function is the function to be calculated on the network object. NOTE: currently only supports statistics calculated on network and igraph objects.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#full_output tells R whether to return simulated distribution in addition to summary statistic



MEMS_saom <- function(model,
                      micro_process,
                      macro_function,
                      object_type=object_type,
                      interval=interval,
                      nsim=nsim,
                      silent=silent,
                      full_output=full_output,
                      SAOM_data=SAOM_data,
                      SAOM_var=SAOM_var,
                      mediator=mediator,
                      algorithm=algorithm,
                      link_id=link_id,
                      controls=controls,
                      control_functions=control_functions) {

  if(class(model)[1]=="sienaGroup"){
    stop("sienaGroup objects not currently supported. To estimate MEMS or AMME with multiple longitudinal networks, re-estimate each siena group independently and then supply the sienaFit objects as a list.")
  }

  coef<-model$theta
  cov_mat<-model$covtheta
  aMEMS_tracer<-0

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



  rate_theta<-MASS::mvrnorm(n=nsim,
                            mu=model$rate,
                            Sigma= diag(model$vrate^2,ncol=length(model$vrate),nrow=length(model$vrate)),
                            empirical=TRUE)

  theta<-cbind(rate_theta,theta)

  output_data<-matrix(NA,nrow=nsim,ncol=length(interval))
  effects_obj<-rbind(attributes(model$f)$condEffects,model$effects) #create effects object
  theta_index<-match(micro_process,effects_obj$effectName) #get index for parameter to be altered
  node_number<-length(SAOM_data$nodeSets$Actors)


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






  ###initiate simulation
  message("Computing MEMS over ",interval[1],"-",interval[length(interval)]," interval")

  for(i in 1:length(interval)){

    #create manipulated parameter--this is equivalent to setting the variable at a specific value
    theta2<-theta
    theta2[,theta_index]<-theta2[,theta_index]*interval[i]
    sim.model<-RSiena::sienaAlgorithmCreate(cond = FALSE,
                                  useStdInits = FALSE, nsub = 0 ,
                                  simOnly = TRUE,
                                  n3=nsim)

    if(silent==FALSE){
      message("Simulating networks with ", micro_process," held at ", interval[i])
    }
    sim_nets<-RSiena::siena07(sim.model,data=SAOM_data,effects=effects_obj,prevAns=model,
                            returnDeps=TRUE,thetaValues = theta2,silent=TRUE)


    #convert networks from edgelist to network object
    net_list<-as.list(vector(length=length(sim_nets$sims),"numeric"))

    #loop over to convert networks from edgelist and then compute output values
    for(j in 1:length(net_list)){

      #convert to network object from edgelist
      net_list[[j]]<-c(list(t(sim_nets$f$Data1$nets$Network[[1]]$mat1)), #starting observed network
                       sim_nets$sims[[j]][[1]][[1]]) #net_list[[i]] contains edge lists for all unique networks for the number of panels minus 1

      for(entry in 1:length(net_list[[j]])){
        A<- matrix(0, nrow=node_number, ncol=node_number)
        b<-net_list[[j]][[entry]]
        A[ b[,1:2] ]<- b[,3]

        net_list[[j]][[entry]]<-network::as.network(A)

        #assign attributes
        if(!is.null(SAOM_var)){

          for(k in 1:length(SAOM_var)){

            if(class(SAOM_var[[k]])[1]=="varCovar"){
              if(is.null(names(SAOM_var)[[k]])){
                SAOM_names<-paste("varCovar",k)}else{
                SAOM_names<-names(SAOM_var)[[k]]}
              network::set.vertex.attribute(net_list[[j]][[entry]],SAOM_names,SAOM_var[[k]][,entry])
            }
            if(class(SAOM_var[[k]])[1]=="varDyadCovar"){
              if(is.null(names(SAOM_var)[[k]])){
                SAOM_names<-paste("varDyadsCovar",k)}else{
                  SAOM_names<-names(SAOM_var)[[k]]}
              network::set.vertex.attribute(net_list[[j]][[entry]],SAOM_names,SAOM_var[[k]][,entry])
            }
          }

        } #close if(!is.null(SAOM var)) statement


         if(length(SAOM_data$cCovars)>0){
           for(k in 1:length(SAOM_data$cCovars)){
            network::set.vertex.attribute(net_list[[j]][[entry]],names(SAOM_data$cCovars)[k],SAOM_data$cCovars[[k]][,1])
          }
        }

       if(length(SAOM_data$dycCovars)>0){
          for(k in 1:length(SAOM_data$dycCovars)){
            network::set.vertex.attribute(net_list[[j]][[entry]],names(SAOM_data$dycCovars)[k],SAOM_data$dycCovars[[k]][,1])
          }
        }

        ###calculate output values
        if(object_type[1]%in%c("network")){

          b<-macro_function(net_list[[j]][[entry]])

          if(!is.null(link_id)){
            node_start<- which(is.na(link_list_data[[j]][,i]))[1]
            node_index<-node_start+(length(b)-1)
            link_list_data[[j]][node_start:node_index,i]<-b
            } ##store vector of output when calling AMME

           if(length(b)>1){
              b<-mean(b,na.rm=TRUE)
              aMEMS_tracer<-1} ##handle node and multiple obs characteristics for MEMS output
          output_data[j,i]<-b


        }else{
          #for igraph functions
            net_list[[j]][[entry]]<-intergraph::asIgraph(net_list[[j]][[entry]])
            b<-macro_function(net_list[[j]][[entry]])

            if(!is.null(link_id)){
              node_start<- which(is.na(link_list_data[[j]][,i]))[1]
              node_index<-node_start+(length(b)-1)
              link_list_data[[j]][node_start:node_index,i]<-b
            } #store vector of output when calling AMME
            if(length(b)>1){
              b<-mean(b,na.rm=TRUE)
              aMEMS_tracer<-1} ##handle node and multiple obs characteristics
            output_data[j,i]<-b
            if(length(controls)>0){
              net_list[[j]][[entry]]<-intergraph::asNetwork(net_list[[j]][[entry]]) # convert back to network for control operations
            }



        }#close macro statistic ifelse


        ##get controls

        if(length(controls)>0){

          for(control in 1:length(controls_list)){
            if(object_type[control+1]%in%c("network")){

             b<-control_functions[[control]](net_list[[j]][[entry]])
             node_start<- which(is.na(controls_list[[control]][[j]][,i]))[1]
             node_index<-node_start+(length(b)-1)
             controls_list[[control]][[j]][node_start:node_index,i]<-b

            }else{
              #for igraph functions
             net_list[[j]][[entry]]<-intergraph::asIgraph(net_list[[j]][[entry]])
             b<-control_functions[[control]](net_list[[j]][[entry]])
             node_start<- which(is.na(controls_list[[control]][[j]][,i]))[1]
             node_index<-node_start+(length(b)-1)
             controls_list[[control]][[j]][node_start:node_index,i]<-b
             net_list[[j]][[entry]]<-intergraph::asNetwork(net_list[[j]][[entry]]) #convert back to network object for further loops


            }

          }
        }#close if controls statement




      }#close entry loop


    }#closes j loop



  }


  if(aMEMS_tracer==1 &is.null(link_id)){
    message("More than one macro statistic is being calculated. Reporting the aMEMS.")
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
