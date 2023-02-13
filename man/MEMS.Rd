\name{MEMS}
\alias{MEMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate the micro effect on macro structure (MEMS).}
\description{
\code{MEMS} implements parametric and nonparametric estimation routines to estimate the
micro effect on macro structure when using a generative network model (i.e., a model
where the dyad, dyad-time period, or dyad-group is the unit of analysis). The MEMS is defined in postestimation as a function of the possibly endogenous micro process \eqn{m}, which is assumed to be a predictor in the micro model of the form \eqn{f(A|m,Z)=\theta m + \gamma ^TZ}, where \eqn{Z} is a matrix of possibly endogenous controls and \eqn{A} is the network of interest. The MEMS when \eqn{\theta} changes from 0 to 1 is given by \deqn{MEMS=\sum_i \frac{M(\theta, m, \gamma, Z)_i-M(\gamma, Z)_i}{n}}, for \eqn{n} observations. Tuning parameters can be assigned to toggle the strength of \eqn{\theta} in model-implied estimates of \eqn{MEMS}. \code{MEMS} currently accepts \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects and implements both parametric and nonparametric estimation. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} objects.

}
\usage{
MEMS(model,
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
      control_functions=NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{the micro-model to be analyzed. Currently accepts \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} objects. To implement pooled estimation, \code{model} should be provided as a list of \code{ergm} or \code{sienaFit} objects.
}
\item{micro_process}{a character string containing the name of the micro process of interest. The character string should exactly match coefficient names in \code{model} output.}
\item{macro_function}{a \code{function} that calculates the macro statistic of interest. Currently accepts user defined functions as well as functions inherent in the \code{\link[igraph]{igraph}} and \code{statnet} packages for \R.}
\item{object_type}{A character string that tells netmediate the type of object to apply the
\code{macro_function} to. Currently accepts \code{\link[igraph]{igraph}} and \code{\link[network]{network}} objects. If left \code{NULL}, \code{\link[network]{network}} objects are assumed. Can be over-ridden to use other object types with a user-function by defining a function that accepts either a \code{\link[network]{network}} or \code{\link[igraph]{igraph}} object and returns a numeric value or vector of numeric values (see examples).}
\item{interval}{The value of tuning parameters to assign to \eqn{\theta}. Should be provided as a vector of numeric values with 2 entries.}
\item{nsim}{The number of simulations or bootstrap samples to use during estimation.}
\item{algorithm}{The estimation algorithm to be used. Currently accepts \code{"parametric"} and \code{"nonparametric"}. If \code{"parametric"}, estimation is obtained with Monte Carlo sampling. If \code{"nonparametric"}, estimation uses bootstrap resampling.}
\item{silent}{logical parameter. Whether to provide updates on the progress of the simulation or not.}
\item{full_output}{logical parameter. If set to \code{TRUE}, the entire distribution of simulated statistics will be provided as part of the model output.}
\item{SAOM_data}{required when the model is a \code{\link[RSiena]{sienaFit}} object; ignored otherwise. If a \code{\link[RSiena]{sienaFit}} object is provided, \code{SAOM_data} should be the \code{siena} object that contains the data for SAOM estimation. If using pooled estimation on multiple \code{\link[RSiena]{sienaFit}} objects (i.e., providing a list of \code{\link[RSiena]{sienaFit}} objects), then \code{SAOM_data} should be provided as an ordered list with each entry containing the \code{\link[RSiena]{siena}} object corresponding to list of \code{\link[RSiena]{sienaFit}} objects.}
\item{SAOM_var}{optional parameter when the model is a \code{\link[RSiena]{sienaFit}} object. \code{SAOM_var} is a list of of the \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects used to assign time varying node and dyad covariates when calling \code{\link[RSiena]{sienaDataCreate}}. If provided, \code{netmediate} assigns the varying node covariates and dyad covariates to each simulated network. This parameter is required when \code{macro_function} computes a statistic that varies as a function of time varying node or dyad covariates (i.e., network segregation, assorativity). Time invariant characteristics (\code{\link[RSiena]{coCovar}} and \code{\link[RSiena]{coDyadCovar}}) are handled internally by \code{MEMS} and should not be provided. When providing a list of \code{\link[RSiena]{sienaFit}} objects for pooled estimation, \code{SAOM_var} should be provided as a list of lists, where each entry in the list contains a list of \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects associated with corresponding \code{\link[RSiena]{sienaFit}} object. }
\item{time_interval}{an optional parameter to be used with \code{\link[relevent]{rem.dyad}} objects. May be provided as a numeric vector or the character string \code{"aggregate"}. If a numeric vector is provided unique network snapshots at each interval. For example, \code{time_interval=c(0,2,3)} would induce two networks, one for the 0 - 2 time period and one for the 2 - 3 time period. If specified as \code{"aggregate"}, the MEMS is calculated by creating an aggregated cross-sectional representation of the entire event sequence. If left \code{NULL}, defaults to |\code{"aggregate"}.}
\item{covar_list}{an optional list of sender/receiver covariates used in \code{\link[relevent]{rem.dyad}} estimation. Only required for \code{\link[relevent]{rem.dyad}} objects when covariates are included. The list format should correspond to the format required by \code{\link[relevent]{rem.dyad}}}.
\item{edgelist}{an optional three column edgelist providing the sender, receiver, and time of event occurrence when using rem.\code{\link[relevent]{rem.dyad}}. Only required when \code{time_interval} is set to \code{NULL} or \code{"aggregate"}. Ignored for other types of models.}
\item{net_logit_y}{the dependent variable for \code{\link[sna]{netlogit}} objects. Should be provided as a vector. Only required when model is a \code{\link[sna]{netlogit}} object.}
\item{net_logit_x}{the matrix of independent variables for \code{\link[sna]{netlogit}} type objects. Only required when model is a \code{\link[sna]{netlogit}} object.}
\item{group_id}{optional vector of group identifiers to use when estimating a \code{glm} or \code{\link[lme4]{glmer}} on grouped data (i.e., multiple time periods, multiple networks). When specified,  \code{MEMS} will induce unique networks for each grouping factor. If left unspecified, all groups/time periods are pooled. If using \code{\link[lme4]{glmer}}, the grouping factor does not have to be provided as part of the model or used as a random effect.}
\item{node_numbers}{a numeric vector containing the number of nodes in each group_id when using \code{glm} or \code{\link[lme4]{glmer}}. If estimating MEMS aggregated over all networks (i.e., \code{group_id=NULL}), this shoud be the total number of nodes in all networks. Required when using \code{glm} or \code{\link[lme4]{glmer}}, ignored otherwise.}
\item{mediator}{a character string detailing the mediator of interest. Intended for internal use with the \code{\link{AMME}} function; not intended for end users.}
\item{link_id}{a vector or list of vectors corresponding to unique identifiers. Intended for internal use with the \code{\link{AMME}} function; not intended for end users. }
\item{controls}{a vector of character strings listing the controls to be calculated when using \code{\link{AMME}}. Intended for internal use with the \code{\link{AMME}} function; not intended for end users. }
\item{control_functions}{a list of functions to calculate the macro control variables provided in controls. Intended for internal use with the \code{\link{AMME}} function; not intended for end users. }

%-end arguments statement
}
\details{
Estimates the MEMS over the provided intervals. If the macro statistic is calculated on the node or subgraph levels or on multiple network observations, the aMEMS is provided instead. Standard errors and confidence intervals are based on the sampling distribution of simulated values, which are calculated either parametrically or nonparametrically according to \code{algorithm}. Parametric estimation is typically faster, but cannot be used for nonparametric network models (e.g., quadratic assignment procedure).

\code{macro_function} is the workhorse component of MEMS. The function should calculate the macro statistic of interest. \code{netmediate} currently supports functions calculated on \code{\link[igraph]{igraph}} and \code{\link[network]{network}} objects, which should be specified as using the \code{object_type} argument. These may be functions inherent to the \code{\link[statnet]{statnet}} and \code{\link[igraph]{igraph}} software package or they may be functions from other packages that accept \code{\link[network]{network}}/\code{\link[igraph]{igraph}} objects. They may also be user-defined functions that accept \code{\link[network]{network}} or \code{\link[igraph]{igraph}} objects as input and return a numeric value or vector of numeric values as output. It is also possible to over-ride the \code{\link[network]{network}} and \code{\link[igraph]{igraph}} object requirements within a user function. To do so, set the \code{object_type} argument to either \code{\link[network]{network}} or \code{\link[igraph]{igraph}} and then define a user-function that accepts a \code{\link[network]{network}} or \code{\link[igraph]{igraph}} object as its input, converts the object to the desired data structure, calculates the statistic of interest, and finally returns a numeric value or vector of numeric values. See examples below for an illustration.

By default, the MEMS is provided by averaging over the distribution of simulated values. If \code{full_output} is set to \code{TRUE}, the entire distribution of simualted statistics is returned. This may be useful when the median or mode of the simulated distribution is required or if the researcher wants to inspect the distributional shape of simulated values.

\code{MEMS} also supports pooled estimation for multiple \code{\link[ergm]{ergm}} or \code{\link[RSiena]{sienaFit}} objects. To use pooled estimation, the model parameter should be specified as a list of \code{\link[ergm]{ergm}} or \code{\link[RSiena]{sienaFit}} objects. If using \code{\link[RSiena]{sienaFit}}, the \code{SAOM_data} argument will also need to be specified as an ordered list with elements corresponding to entries in the list of \code{\link[RSiena]{sienaFit}} objects. Similarly, the \code{SAOM_var} parameter will need to be specified as a list of lists, where each entry in the list is, itself, a list containing all \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects used to calculate macro statistics of interest. Note that \code{SAOM_var} should not be provided if the macro statistic of interest is not a function of the variables contained in \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}}.

When estimating a relational event model with a \code{\link[relevent]{rem.dyad}} object, \code{time_interval} can be specified to provide exact time intervals over which to induce unique networks. This utility is often useful when combining \code{\link[relevent]{rem.dyad}} estimation with \code{\link{AMME}} when the \code{macro_model} is panel data with coarse timing information. The same behavior can be obtained when estimating a relational event model using \code{glm} or \code{\link[lme4]{glmer}} by assigning the desired time intervals in the model matrix and then providing the vector of time intervals to the \code{group_id} parameter when calling \code{MEMS}.


}
\value{
If \code{full_output=FALSE}, then a table is returned with the MEMS, its standard error, confidence interval, and p-value.

If \code{full_output=TRUE}, then a list is returned with the following three elements.
\item{summary_dat}{is the table of summary output containing the MEMS, its standard error, confidence interval, and p-value.}
\item{output_data}{is a matrix where each row is a simulated draw of the MEMS (or a simulation draw for a specific network in the case of temporal data or pooled estimation) and each column corresponds to a unique value provided in the interval argument.}
\item{mems_samples}{is vector matrix corresponding where each row is a simulated draw of the MEM (or a simulation draw for a specific network in the case of temporal data or pooled estimation) and each column represents the differences in MEMS/aMEMS when subtracting the value of a macro statistic at one interval level from the next highest interval level.}
}
\references{
Duxbury, Scott W. "Micro Effects on Macro Structure in Social Networks." Working paper.

}
\author{
Duxbury, Scott W. Assistant Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}

\seealso{
\code{\link{AMME}}
\code{\link[ergMargins]{ergm.mma}}
\code{\link[mediation]{mediate}}
}
\examples{
\dontshow{
require(ergm)
require(network)
require(sna)


set.seed(21093)
a1<-network::as.network(matrix(c(rbinom(10, 1,.3),
            rbinom(10, 1,.3),
             rbinom(10, 1,.3),
              rbinom(10, 1,.3),
               rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3)),
          nrow=10,ncol=10))

network::set.vertex.attribute(a1,"var.1",rbinom(10,1,.3))
a<-ergm(a1~edges+nodeifactor("var.1")+nodeofactor("var.1"))

MEMS(a,
      micro_process="nodeifactor.var.1.1",
      macro_function=gtrans,
      nsim=50)

}

\donttest{
########################################
# ERGM examples and basic utilities
#######################################


####start with a simple model
library(statnet)

data("faux.mesa.high")

model1<-ergm(faux.mesa.high~edges+
               nodecov("Grade")+
               nodefactor("Race")+
               nodefactor("Sex")+
               nodematch("Race")+
               nodematch("Sex")+
               absdiff("Grade"),
             control = control.ergm(parallel=4))



##calculate the MEMS when the absolute difference in grade is changed from an interval of 0 to 1
  #with default specifications for gtrans
MEMS(model1,
     micro_process="absdiff.Grade",
     macro_function = gtrans,
     object_type = "network",
     nsim=100,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "parametric")

#call an argument from gtrans by specifying it as a function
  #use nonparametric estimation
MEMS(model1,
     micro_process="absdiff.Grade",
     macro_function = function(x){gtrans(x,measure="strongcensus")},
     object_type = "network",
     nsim=100,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "nonparametric")




####calculate the MEMS using igraph
MEMS(model1,
     micro_process="absdiff.Grade",
     macro_function = function(x){igraph::transitivity(x,type="local")},
     object_type = "igraph",
     nsim=100,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "parametric")



##specify a user function that counts the number of communities
community_counts<-function(x){
  walktrap<-igraph::walktrap.community(x) #use walktrap community detection
  return(length(unique(walktrap$membership))) #return the number of communities
}

MEMS(model1,
     micro_process="absdiff.Grade",
     macro_function = community_counts,
     object_type = "igraph",
     nsim=100,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "parametric")



##calculate a function using exogenous node attributes
assortativity_grade<-function(x){
  require(igraph)
  return(assortativity_nominal(x,V(x)$Grade))
}

MEMS(model1,
     micro_process="absdiff.Grade",
     macro_function = assortativity_grade,
     object_type = "igraph",
     nsim=100,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "parametric")

##specify a user function that does not depend on either igraph or statnet
  #assuming a network input object, we have
manual_user_function<-function(x){
  x<-as.sociomatrix(x)
  return(colSums(x))
}

MEMS(model1,
     micro_process="absdiff.Grade",
     macro_function = manual_user_function,
     object_type = "network",
     nsim=100,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "parametric")







####estimation for POOLED ERGM
data("faux.magnolia.high")

model2<-ergm(faux.magnolia.high~edges+
               nodecov("Grade")+
               nodefactor("Race")+
               nodefactor("Sex")+
               nodematch("Race")+
               nodematch("Sex")+
               absdiff("Grade"),
             control = control.ergm(parallel=4))



MEMS(list(model1,model2),
     micro_process="absdiff.Grade",
     macro_function = assortativity_grade,
     object_type = "igraph",
     nsim=50,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "parametric")



#################################
#   Estimation with GLM and GLMER
#################################
library(btergm)

#use models 1 and 2 from examples above
glm_dat<-edgeprob(model1)
glm_dat2<-edgeprob(model2)
glm_dat2<-glm_dat2[,-c(4)]


##create stacked dataset for the purposes of grouped estimation
glm_dat$net_id<-"mesa" #specify ID for each network
glm_dat2$net_id<-"magnolia"
glm_dat<-rbind(glm_dat,glm_dat2)


##estimate as a linear probability model
net_glm<-glm(tie~nodecov.Grade+
               nodefactor.Race.Hisp+
               nodefactor.Race.NatAm+
               nodefactor.Race.Other+
               nodefactor.Sex.M+
               nodematch.Race+
               nodematch.Sex+
               absdiff.Grade,
             data=glm_dat)



MEMS(net_glm,
     micro_process="nodematch.Race", #should be written as in netlogit output
     macro_function = function(x){gtrans(x)},
     object_type = "network",
     nsim=100,
     interval=c(0,.5),
     silent=FALSE,
     full_output = FALSE,
     algorithm = "parametric",
     group_id=glm_dat$net_id, #provide network ID for estimation
     node_numbers =c(network.size(faux.mesa.high), #provide the number of nodes in each network
                      network.size(faux.magnolia.high)))


##estimate as a multilevel model
library(lme4)
net_glmer<-glmer(tie~nodecov.Grade+
               nodefactor.Race.Hisp+
               nodefactor.Race.NatAm+
               nodefactor.Race.Other+
               nodefactor.Sex.M+
               nodematch.Race+
               nodematch.Sex+
               absdiff.Grade+
                 (1|net_id),
             data=glm_dat,
             family=gaussian)



MEMS(net_glmer,
     micro_process="nodematch.Race", #should be written as in netlogit output
     macro_function = function(x){gtrans(x)},
     object_type = "network",
     nsim=50,
     interval=c(0,.5),
     silent=FALSE,
     full_output = FALSE,
     algorithm = "parametric",
     group_id=glm_dat$net_id,
     node_numbers =c(203,974))




##############################################
##nonparametric estimation for bootstrap TERGM
##############################################

library(btergm)
data(alliances)
ally_data<-list(LSP[[1]],
                LSP[[2]],
                LSP[[3]])

bt_model<-btergm(ally_data~edges+
                   gwesp(.7,fixed=T)+
                   mutual,R=200)



MEMS(bt_model,
     micro_process="gwesp.fixed.0.7",
     macro_function = gtrans,
     object_type = "network",
     nsim=50,
     interval=c(0,1),
     silent=FALSE,
     algorithm = "nonparametric")





################################
# Parametric estimation using SAOM
##################################
library(RSiena)
#specify 3 wave network panel data as DV
network_list<-array(c(s501,s502,s503),dim = c(50,50,3))

Network<-sienaDependent(network_list)
Smoking<-varCovar(s50s)
Alcohol<-varCovar(s50a)
SAOM.Data<-sienaDataCreate(Network=Network,Smoking,Alcohol)

#specify
SAOM.terms<-getEffects(SAOM.Data)
SAOM.terms<-includeEffects(SAOM.terms,egoX,altX,sameX,interaction1="Alcohol")
SAOM.terms<-includeEffects(SAOM.terms,egoX,altX,sameX,interaction1="Smoking")
SAOM.terms<-includeEffects(SAOM.terms,transTies,inPop)


create.model<-sienaAlgorithmCreate(projname="netmediate",
                                   nsub=5,
                                   n3=2000)


##estimate the model using siena07
SAOM_model<-siena07(create.model,
                        data=SAOM.Data,
                        effects=SAOM.terms,
                        verbose=TRUE)


SAOM_model




##basic specification for reciprocity effects on outdegree distribution
MEMS(SAOM_model,
     micro_process="reciprocity", #should be written as in SIENA output
     macro_function = function(x){igraph::degree(x,mode="out")},
     object_type = "igraph",
     interval=c(0,.5),
     SAOM_data=SAOM.Data,
     silent=FALSE,
     algorithm = "parametric")



##include user functions on time varying covariates
assortativity_smoking<-function(x){
  return(assortativity_nominal(x,V(x)$Smoking))
}


MEMS(SAOM_model,
     micro_process="reciprocity",
     macro_function = assortativity_smoking,
     object_type = "igraph",
     interval=c(0,.5),
     SAOM_data=SAOM.Data,
     SAOM_var=list(Smoking=Smoking,Alcohol=Alcohol), #Smoking and Alcohol are varCovar objects
     silent=FALSE,
     full_output = FALSE,
     algorithm = "parametric")




###Pooled SAOM
MEMS(list(SAOM_model,SAOM_model),
     micro_process="reciprocity",
     macro_function = gtrans,
     object_type = "network",
     interval=c(0,.5),
     SAOM_data=list(SAOM.Data,SAOM.Data),
     silent=FALSE,
     full_output = FALSE,
     nsim=100,
     algorithm = "parametric")


#Pooled SAOM with user functions and time varying attributes
assortativity_smoking<-function(x){
  return(assortativity_nominal(x,V(x)$Smoking))
}



MEMS(list(SAOM_model,SAOM_model),
     micro_process="reciprocity",
     macro_function = assortativity_smoking,
     object_type = "igraph",
     interval=c(0,.5),
     SAOM_data=list(SAOM.Data,SAOM.Data),
     SAOM_var=list(list(Smoking=Smoking,Alcohol=Alcohol),
                    list(Smoking=Smoking,Alcohol=Alcohol)),
     silent=FALSE,
     full_output = FALSE,
     nsim=100,
     algorithm = "parametric")







#######################################
# Relational event models using relevent
#######################################
set.seed(21093)
library(relevent)
##generate a network with 15 discrete time periods
  #example based on relevent rem.dyad example
library(relevent)
roweff<-rnorm(10) #Build rate matrix
roweff<-roweff-roweff[1] #Adjust for later convenience
coleff<-rnorm(10)
coleff<-coleff-coleff[1]
lambda<-exp(outer(roweff,coleff,"+"))
diag(lambda)<-0
ratesum<-sum(lambda)
esnd<-as.vector(row(lambda)) #List of senders/receivers
erec<-as.vector(col(lambda))
time<-0
edgelist<-vector()
while(time<15){ # Observe the system for 15 time units
  drawsr<-sample(1:100,1,prob=as.vector(lambda)) #Draw from model
  time<-time+rexp(1,ratesum)
  if(time<=15) #Censor at 15
    edgelist<-rbind(edgelist,c(time,esnd[drawsr],erec[drawsr]))
  else
    edgelist<-rbind(edgelist,c(15,NA,NA))
}
effects<-c("CovSnd","FERec")



##estimate model
fit.time<-rem.dyad(edgelist,10,effects=effects,
                   covar=list(CovSnd=roweff),
                   ordinal=FALSE,hessian=TRUE)


###aggregate estimation
MEMS(fit.time,
     micro_process="CovSnd.1", #should be written as in relevent output
     macro_function = function(x){sna::degree(x)},
     object_type = "network",
     nsim=10,
     interval=c(0,.5),
     silent=FALSE,
     covar_list=list(CovSnd=roweff), #covariate effects
     time_interval="aggregate", ##aggregated estimation
     edgelist=edgelist,
     algorithm = "parametric")


##time interval estimation
##estimation with time intervals
MEMS(fit.time,
     micro_process="CovSnd.1",
     macro_function = function(x){igraph::degree(x)},
     object_type = "igraph",
     nsim=10,
     interval=c(0,.1),
     silent=TRUE,
     covar_list=list(CovSnd=roweff),
     time_interval=c(0,5,10,15), #specify three time intervals, 0 - 5, 5 - 10, and 10 - 15
     algorithm = "parametric")







########################################################
# Network regression with quadratic assignment procedure
########################################################
library(sna)
##generate network data
set.seed(21093)
x<-rgraph(20,4)
y.l<-x[1,,]+4*x[2,,]+2*x[3,,]
y.p<-apply(y.l,c(1,2),function(a){1/(1+exp(-a))})
y<-rgraph(20,tprob=y.p)

nl<-netlogit(y,x,reps=100)
summary(nl)



MEMS(nl,
     micro_process="x2", #should be written as in netlogit output
     macro_function = function(x){degree(x)},
     object_type = "igraph",
     nsim=20,
     interval=c(0,1),
     silent=FALSE,
     full_output = FALSE,
     net_logit_y=y,
     net_logit_x=x,
     algorithm = "nonparametric")













}%close don't test

}%close examples
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~mediation }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~networks }% __ONLY ONE__ keyword per line
\keyword{ ~macro }% __ONLY ONE__ keyword per line
\keyword{ ~micro }% __ONLY ONE__ keyword per line
