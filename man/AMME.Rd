\name{AMME}
\alias{AMME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate the average micro mediated effect (AMME).}
\description{
\code{AMME} implements parametric and nonparametric estimation routines to estimate the
average mediated micro effect. It requires two models. The first is a generative network model (i.e., a model where the dyad, dyad-time period, or dyad-group is the unit of analysis) of the form \eqn{f(A|m,Z)=\theta m + \gamma ^TZ}, where \eqn{A} is a cross-sectional or longitudinal network or group of longitudinal or cross-sectional networks, \eqn{m} is the possibly endogenous micro process of interest and \eqn{Z} is a matrix of possibly endogenous micro level control variables.

The second model is a cross-sectional or longitudinal macro model (i.e., a model where the unit of analysis is a node, subgraph, or network or a combination of nodes, subgraphs, and networks measured collected from multiple settings [such as distinct schools or organizations]) of the form \eqn{g(y|M,X)=bM+\beta ^TX}, where \eqn{y} is the outcome variable, \eqn{M} is the mediating macro variable, and \eqn{X} is a matrix of control variables that possibly vary as a function of microprocess \eqn{m}. For a single macro statistic, the AMME when \eqn{m} is increased from 0 to 1 is given by \deqn{AMME=\frac{1}{nK} y_i(M(\theta, m, \gamma, Z),X)-y_i(M(\gamma, Z),X)}, where \eqn{n} is the number of observations and \eqn{K=2} is the number of possible combinations of values of \eqn{m} and \eqn{M}. Tuning parameters can be specified to vary the strength of \eqn{\theta}. AMME currently accepts the following micro models: \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects. The following macro model objects are accepted: \code{lm}, \code{glm}, \code{lmer}, \code{glmer}, \code{gam}, \code{plm}, and \code{lnam} objects. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} micro models. Both parametric and nonparametric estimation are available.

}
\usage{
AMME(micro_model,
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
      node_numbers=NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{micro_model}{the micro-model. Currently accepts \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} objects. To implement pooled estimation, \code{model} should be provided as a list of \code{ergm} or \code{sienaFit} objects.
}
\item{macro_model}{the macro model. Currently accepts \code{lm}, \code{glm}, \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}}, \code{\link[gam]{gam}}, \code{\link[plm]{plm}}, and \code{\link[sna]{lnam}} objects.
}
\item{micro_process}{a character string containing the name of the micro process of interest. The character string should exactly match the relevant coefficient name in \code{micro_model} output.
}
\item{mediator}{a character string containing the name of the mediating variable of interest.  The character string should exactly match the relevant coefficient name in \code{macro_model} output.
}
\item{macro_function}{a \code{function} that calculates \code{mediator} on the simulated networks. Currently accepts user defined functions as well as functions inherent in the \code{\link[igraph]{igraph}} and \code{statnet} packages for \R.
}
\item{link_id}{a required vector of IDs used to link the \code{micro_model} output to the \code{macro_model} input. If calculating a network-level mediator, this should be the network identifier or network-group/network-time period identifier. If calculating a node-level mediator, this should be the node ID or node-time-period/node-group identifier. Observations should correspond exactly to rows in the \code{macro_model} data matrix. If calculating multiple network statistics at different levels of analysis when \code{controls} are included, \code{link_id} may be provided as an ordered list of identifiers. In this case, each entry in the list is a vector of IDs corresponding to the unique entries of the relevant statistics. If provided as a list, the first entry should correspond to \code{macro_function} (i.e., the mediator) and the remaining entries should correspond to \code{control_functions} (i.e., the controls).
 }
 \item{controls}{a vector of character strings listing the control variables in \code{macro_model} that may vary as a function of \code{micro_process}. Each element in \code{controls} should correspond exactly to a coefficient in \code{macro_model} output. If \code{controls} is left \code{NULL},then the AMME is calculated without controlling for confounding network variables.
 }
\item{control_functions}{a list of functions used to calculate \code{controls}. The elements in \code{control_functions} should correspond exactly to the elements in \code{controls} and should be provided in the same order. If \code{micro_process} appears as an independent variable in \code{macro_model}, then this can be specified by specifying the \code{netmediate} helper function \code{identity_function} to \code{control_functions}.
}
\item{object_type}{A character string or vector of character strings that tells \code{netmediate} the type of object to apply the \code{macro_function} and \code{control_functions} to. If \code{controls} are included into the \code{AMME} call, then \code{object_type} should be provided as a vector of character strings where the first element is the \code{object_type} for \code{macro_function} and the remaining elements are the ordered \code{object_type} for \code{control_functions}. Currently accepts \code{\link[igraph]{igraph}} and \code{\link[network]{network}} objects. If left \code{NULL}, \code{\link[network]{network}} objects are assumed. Can be over-ridden to use other object types with a user-function by defining a function that accepts either a \code{\link[network]{network}} or \code{\link[igraph]{igraph}} object and returns a numeric value or vector of numeric values (see examples).
}
\item{interval}{Tuning parameters to vary the strength of \eqn{\theta}. Should be provided as a vector of numeric values with 2 entries.
}
\item{nsim}{The number of simulations or bootstrap samples to use during estimation.
}
\item{algorithm}{The estimation algorithm to be used. Currently accepts \code{"parametric"} and \code{"nonparametric"}. If \code{"parametric"}, estimation is obtained with Monte Carlo sampling. If \code{"nonparametric"}, estimation uses bootstrap resampling.
}
\item{silent}{logical parameter. Whether to provide updates on the progress of the simulation or not.
}
\item{full_output}{logical parameter. If set to \code{TRUE}, the entire distribution of simulated statistics will be provided as part of the model output.
}
\item{SAOM_data}{required when \code{micro_model} is a \code{\link[RSiena]{sienaFit}} object; ignored otherwise. If a \code{\link[RSiena]{sienaFit}} object is provided, \code{SAOM_data} should be the \code{siena} object that contains the data for SAOM estimation. If using pooled estimation on multiple \code{\link[RSiena]{sienaFit}} objects (i.e., providing a list of \code{\link[RSiena]{sienaFit}} objects), then \code{SAOM_data} should be provided as an ordered list with each entry containing the \code{\link[RSiena]{siena}} object corresponding to list of \code{\link[RSiena]{sienaFit}} objects.
}
\item{SAOM_var}{optional parameter when \code{micro_model} is a \code{\link[RSiena]{sienaFit}} object. \code{SAOM_var} is a list of of the \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects used to assign time varying node and dyad covariates when calling \code{\link[RSiena]{sienaDataCreate}}. If provided, \code{netmediate} assigns the varying node covariates and dyad covariates to each simulated network. This parameter is required when \code{macro_function} computes a statistic that varies as a function of time varying node or dyad covariates (i.e., network segregation, assorativity). Time invariant characteristics (\code{\link[RSiena]{coCovar}} and \code{\link[RSiena]{coDyadCovar}}) are handled internally by \code{MEMS} and should not be provided. When providing a list of \code{\link[RSiena]{sienaFit}} objects for pooled estimation, \code{SAOM_var} should be provided as a list of lists, where each entry in the list contains a list of \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects associated with corresponding \code{\link[RSiena]{sienaFit}} object.
}
\item{time_interval}{an optional parameter to be used when \code{micro_model} is a  \code{\link[relevent]{rem.dyad}} object. May be provided as a numeric vector or the character string \code{"aggregate"}. If a numeric vector is provided unique network snapshots at each interval. For example, \code{time_interval=c(0,2,3)} would induce two networks, one for the 0 - 2 time period and one for the 2 - 3 time period. If specified as \code{"aggregate"}, the AMME is calculated by creating an aggregated cross-sectional representation of the entire event sequence. If left \code{NULL}, defaults to |\code{"aggregate"}. Note that \code{time_interval} must correspond to the time periods observed in \code{macro_model}. That is, \code{time_interval} must be set to \code{"aggregate"} when \code{macro_model} is cross-sectional and the entries in \code{time_interval} must correspond to the time periods observed in the repated measurement data when \code{macro_model} is longitudinal.
}
\item{covar_list}{an optional list of sender/receiver covariates used in \code{\link[relevent]{rem.dyad}} estimation. Only required when a \code{\link[relevent]{rem.dyad}} object is the \code{micro_model} and covariates are in the \code{\link[relevent]{rem.dyad}} call. The list format should correspond to the format required by \code{\link[relevent]{rem.dyad}}.
}
\item{edgelist}{an optional three column edgelist providing the sender, receiver, and time of event occurrence when \code{micro_model} is a \code{\link[relevent]{rem.dyad}} object. Only required when \code{time_interval} is set to \code{NULL} or \code{"aggregate"}. Ignored for other types of models.
}
\item{net_logit_y}{the dependent variable when \code{micro_model} is a \code{\link[sna]{netlogit}} object. Should be provided as a vector.
}
\item{net_logit_x}{the matrix of independent variables when \code{micro_model} is a \code{\link[sna]{netlogit}} object
}
\item{group_id}{optional vector of group identifiers to use when \code{micro_model} is a \code{glm} or \code{\link[lme4]{glmer}} on grouped data (i.e., multiple time periods, multiple networks). When specified,  \code{AMME} will induce unique networks for each grouping factor. If left unspecified, all groups/time periods are pooled. If using \code{\link[lme4]{glmer}}, the grouping factor does not have to be provided as part of the model or used as a random effect. If specified, the entries in the \code{macro_model} model matrix are assumed to be sequentially ordered by unit_id-group_id.
}
\item{node_numbers}{a numeric vector containing the number of nodes in each group_id when using \code{glm} or \code{\link[lme4]{glmer}}. If estimating AMME aggregated over all networks (i.e., \code{group_id=NULL}), this shoud be the total number of nodes in all networks. Required when using \code{glm} or \code{\link[lme4]{glmer}}, ignored otherwise.
}


%-end arguments statement
}
\details{
Estimates the AMME over the provided intervals. Standard errors and confidence intervals are based on the sampling distribution of simulated values, which are calculated either parametrically or nonparametrically according to \code{algorithm}. Parametric estimation is typically faster, but cannot be used for nonparametric network models (e.g., quadratic assignment procedure).

\code{macro_function} and \code{control_functions} make up the core utilites of \code{AMME}. \code{macro_function} calculates the mediating variable of interest, while \code{control_functions} calculates all control variables that vary as a function of \code{micro_process} and potentially confound the effect of \code{mediator}. When \code{controls} are left \code{NULL}, then \code{AMME} estimates the AMME without accounting for confounding variables. Specifying \code{controls} and \code{control_functions} ensures that estimates of the AMME account for alternative pathways from \code{micro_process} to the outcome variable in \code{macro_model}. In cases where \code{micro_process} is included as a predictor variable in \code{macro_model}, this can be specified by including the \code{netmediate} helper function \code{identity_function} into \code{control_functions}.

\code{netmediate} currently supports functions calculated on \code{\link[igraph]{igraph}} and \code{\link[network]{network}} objects, which should be specified using the \code{object_type} argument. These may be functions inherent to the \code{\link[statnet]{statnet}} and \code{\link[igraph]{igraph}} software package or they may be functions from other packages that accept \code{\link[network]{network}}/\code{\link[igraph]{igraph}} objects. The functions provided to \code{macro_function} and \code{control_functions} may also be user-defined functions that accept \code{\link[network]{network}} or \code{\link[igraph]{igraph}} objects as inputs and return a numeric value or vector of numeric values as output. It is also possible to over-ride the \code{\link[network]{network}} and \code{\link[igraph]{igraph}} object requirements within a user function. To do so, set the \code{object_type} argument (or relevant element within the \code{object_type} argument when \code{object_type} is a list) to either \code{\link[network]{network}} or \code{\link[igraph]{igraph}} and then define a user-function that accepts a \code{\link[network]{network}} or \code{\link[igraph]{igraph}} object as its input, converts the object to the desired data structure, calculates the statistic of interest, and  returns a numeric value or vector of numeric values. See examples below for an illustration.

By default, the AMME is calculated by averaging over the distribution of simulated values. If \code{full_output} is set to \code{TRUE}, the distribution of simualted statistics is returned. This may be useful when the median or mode of the simulated distribution is required or if the researcher wants to inspect the distributional shape of simulated values.

\code{AMME} also supports pooled estimation for when multiple \code{\link[ergm]{ergm}} or \code{\link[RSiena]{sienaFit}} objects are used as the \code{micro_model}. To use pooled estimation, the model parameter should be specified as a list of \code{\link[ergm]{ergm}} or \code{\link[RSiena]{sienaFit}} objects. If using \code{\link[RSiena]{sienaFit}}, the \code{SAOM_data} argument will also need to be specified as an ordered list with elements corresponding to entries in the list of \code{\link[RSiena]{sienaFit}} objects. Similarly, the \code{SAOM_var} parameter will need to be specified as a list of lists, where each entry in the list is, itself, a list containing all \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects used to calculate macro statistics of interest. Note that \code{SAOM_var} should not be provided if the macro statistic of interest is not a function of the variables contained in \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}}.


}
\value{
If \code{full_output=FALSE}, then a table is returned with the AMME, its standard error, confidence interval, and p-value.

If \code{full_output=TRUE}, then a list is returned with the following three elements.
\item{summary_dat}{is the table of summary output ucontaining the AMME, its standard error, confidence interval, and p-value.}
\item{AMME_obs}{is vector of observations where each entry is the AMME for a single simulation trial.}
\item{prop_explained_obs}{is vector containing the proportion explained values for each simulation trial.}
}

\references{
Duxbury, Scott W. "Micro-macro Indirect Effects in Social Networks." Working paper.
}
\author{
Duxbury, Scott W. Assistant Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}

\seealso{
\code{\link{MEMS}}
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

node_data<-data.frame(outdegree=rowSums(network::as.sociomatrix(a1)),
                      var1=rnorm(10))

b<-lm(var1~outdegree,data=node_data)
AMME(micro_model=a,
     macro_model=b,
     micro_process="nodeifactor.var.1.1",
     mediator="outdegree",
     macro_function=function(x){colSums(network::as.sociomatrix(x))},
     link_id=1:10,
     nsim=20)



}
\donttest{
##############################
#   Basic AMME specifications
#############################


####create ERGM generative model
library(statnet)
data("faux.mesa.high")
ergm_model<-ergm(faux.mesa.high~edges+
                   nodecov("Grade")+
                   nodefactor("Race")+
                   nodefactor("Sex")+
                   nodematch("Race")+
                   nodematch("Sex")+
                   absdiff("Grade"),
                 control = control.ergm(parallel=4))


###create node-level data for second stage analysis with
node_level_data<-data.frame(grade=faux.mesa.high\%v\%"Grade",
                            race=faux.mesa.high\%v\%"Race",
                            sex=faux.mesa.high\%v\%"Sex",
                            degree=degree(faux.mesa.high))

node_level_data$senior<-0
node_level_data$senior[node_level_data$grade==max(node_level_data$grade)]<-1
node_level_data$v_id<-1:network.size(faux.mesa.high) #define ID for each observation

probit_model<-glm(senior~race+sex+degree,
                data=node_level_data,
                family=binomial(link="probit"))

###estimate the indirect effect of grade homophily on senior status acting through degree centrality
  #in a model with no network control variables
AMME(micro_model=ergm_model,
     macro_model=probit_model,
     micro_process="absdiff.Grade",
     mediator="degree",
     macro_function=degree,
     link_id=node_level_data$v_id, #specify vertex IDs
     object_type="network",
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE)

#use nonparametric estimation for a generalized additive model
library(gam)

gam_model<-gam(senior~race+sex+s(degree),
               data=node_level_data)

AMME(micro_model=ergm_model,
     macro_model=gam_model,
     micro_process="absdiff.Grade",
     mediator="s(degree)",
     macro_function=degree,
     link_id=node_level_data$v_id,
     object_type="network",
     interval=c(0,1),
     nsim=50,
     algorithm="nonparametric",
     silent=FALSE)



###estimate AMME with linear network autocorrelation model

lnam_model<-lnam(node_level_data$grade,
                 x=as.matrix(node_level_data[,4:5]),
                 W1=as.sociomatrix(faux.mesa.high))


AMME(micro_model=ergm_model,
     macro_model=lnam_model,
     micro_process="absdiff.Grade",
     mediator="degree",
     macro_function=degree,
     link_id=node_level_data$v_id,
     object_type="network",
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE)




############################
#   Including controls
###########################

##single control
node_level_data<-data.frame(grade=faux.mesa.high\%v\%"Grade",
                            race=faux.mesa.high\%v\%"Race",
                            sex=faux.mesa.high\%v\%"Sex",
                            degree=degree(faux.mesa.high),
                            betweenness=betweenness(faux.mesa.high))

node_level_data$senior<-0
node_level_data$senior[node_level_data$grade==max(node_level_data$grade)]<-1
node_level_data$v_id<-1:network.size(faux.mesa.high) #define ID for each observation

probit_model<-glm(senior~race+sex+degree+betweenness,
                  data=node_level_data,
                  family=binomial(link="probit"))


AMME(micro_model=ergm_model,
     macro_model=probit_model,
     micro_process="absdiff.Grade",
     mediator="degree",
     macro_function=degree,
     link_id=node_level_data$v_id, #specify vertex IDs
     controls="betweenness", #should match model output exactly
     control_functions=betweenness,
     object_type="network",
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE)



##multiple controls
##include an AR 1 parameter to make it a nonlinear network autocorrelation model
node_level_data$AR1<-as.sociomatrix(faux.mesa.high)\%*\%node_level_data$senior
probit_model<-glm(senior~race+sex+degree+betweenness+AR1,
                  data=node_level_data,
                  family=binomial(link="probit"))

#specify user function
ar_function<-function(x){
  return(as.sociomatrix(x)\%*\%node_level_data$senior)
}


AMME(micro_model=ergm_model,
     macro_model=probit_model,
     micro_process="absdiff.Grade",
     mediator="degree",
     macro_function=degree,
     link_id=node_level_data$v_id,
     controls=c("betweenness","AR1"), #should match model output exactly
     control_functions=list(betweenness,ar_function), #provide functions as a list
     object_type="network",
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE)



##using identity_function when micro_process has a direct effect on y
  #to use identity_function, the control and micro_process need to have the same
  #name and the macro control variable has to be numeric

node_level_data$Sex<-as.numeric(as.factor(node_level_data$sex))
logit_model<-glm(senior~race+Sex+degree+betweenness+AR1,
                  data=node_level_data,
                  family=binomial)



AMME(micro_model=ergm_model,
     macro_model=logit_model,
     micro_process="nodefactor.Sex.M",
     mediator="degree",
     macro_function=degree,
     link_id=node_level_data$v_id,
     controls=c("betweenness","AR1","Sex"), #should match model output exactly
     control_functions=list(betweenness,ar_function,identity_function),
     object_type="network",
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE)





################################
#   More complex data structures
###############################


###############################
# AMME with longitudinal data
##############################

#bootstrap TERGM and panel data model
library(btergm)
library(plm)
data(alliances)

ally_data<-list(LSP[[1]],
                LSP[[2]],
                LSP[[3]])

#fit bootstrap TERGM with 200 replications
bt_model<-btergm(ally_data~edges+
                   gwesp(.7,fixed=T)+
                   mutual,R=200)


#create node data
ally_node_data<-data.frame(outdeg=c(rowSums(LSP[[1]]),rowSums(LSP[[2]]),rowSums(LSP[[3]])),
                           indeg=c(colSums(LSP[[1]]),colSums(LSP[[2]]),colSums(LSP[[3]])))

ally_node_data$v_id<-rep(rownames(LSP[[1]]),3) #create node IDS
ally_node_data$t_id<-c(rep(1, nrow(ally_data[[1]])), #create time IDS
                       rep(2, nrow(ally_data[[1]])),
                       rep(3, nrow(ally_data[[1]])))
ally_node_data$link_id<-paste(ally_node_data$v_id,ally_node_data$t_id)#create node-panel identifiers

ally_node_data$v_id<-as.factor(as.character(ally_node_data$v_id))

#estimate a linear model with node fixed effects
lm_model<- lm(outdeg~indeg +v_id,
          data = ally_node_data)



AMME(micro_model=bt_model,
     macro_model=lm_model,
     micro_process="gwesp.fixed.0.7",
     mediator="indeg",
     macro_function=function(x){degree(x,cmode="indegree")},
     link_id=ally_node_data$link_id, #provide node-panel identifiers
     object_type="network",
     interval=c(0,1),
     nsim=11,
     algorithm="nonparametric",
     silent=FALSE)


##include controls at different units of analysis
  #include global transitivity statistic at each network panel
transitivity_list<-c(gtrans(as.network(LSP[[1]])),
                     gtrans(as.network(LSP[[2]])),
                     gtrans(as.network(LSP[[3]])))


ally_node_data$transitivity<-c(rep(transitivity_list[1],nrow(LSP[[1]])),
                               rep(transitivity_list[2],nrow(LSP[[2]])),
                               rep(transitivity_list[3],nrow(LSP[[3]])))



lm_model<- lm(outdeg~indeg+transitivity +v_id,
              data = ally_node_data)



AMME(micro_model=bt_model,
     macro_model=lm_model,
     micro_process="gwesp.fixed.0.7",
     mediator="indeg",
     macro_function=function(x){degree(x,cmode="indegree")},
     link_id=list(ally_node_data$link_id,ally_node_data$t_id),#list of IDs for nodes and time
     controls="transitivity",
     control_functions = gtrans,
     object_type="network",
     interval=c(0,1),
     nsim=11,
     algorithm="nonparametric",
     silent=FALSE)




#SAOM and panel data model with PLM package
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


##estimate the SAOM
SAOM_model<-siena07(create.model,
                        data=SAOM.Data,
                        effects=SAOM.terms,
                        verbose=TRUE)


##create node-level data
node_level_data<-data.frame(smoking=s50s[,1], #smoking behavior for DV
                            alcohol=s50a[,1],
                            v_id=rownames(s501), #unique node IDS
                            wave="Wave 1",       #unique time IDS
                            outdegree=rowSums(s501),
                            indegree=colSums(s501),
                            AR1=s501\%*\%s50s[,1],  #assign network autocorrelation
                            gcc=gtrans(as.network(s501)))

node_level_data<-rbind(node_level_data,data.frame(smoking=s50s[,2],
                                                  alcohol=s50a[,2],
                                                  v_id=rownames(s502),
                                                  wave="Wave 2",
                                                  outdegree=rowSums(s502),
                                                  indegree=colSums(s502),
                                                  AR1=s502\%*\%s50s[,2],
                                                  gcc=gtrans(as.network(s502))))



node_level_data<-rbind(node_level_data,data.frame(smoking=s50s[,3],
                                                  alcohol=s50a[,3],
                                                  v_id=rownames(s503),
                                                  wave="Wave 3",
                                                  outdegree=rowSums(s503),
                                                  indegree=colSums(s503),
                                                  AR1=s503\%*\%s50s[,3],
                                                  gcc=gtrans(as.network(s503))))


##create unique identifiers for node-panel
node_level_data$unique_ids<-paste(node_level_data$v_id,node_level_data$wave)

##estimate one-way fixed effects model with PLM
library(plm)
FE_model<-plm(smoking~alcohol+outdegree+indegree+AR1+gcc,
               data=node_level_data,
               index=c("v_id","wave"))



##create AR function to provide to AMME
ar_function<-function(x){return(as.sociomatrix(x)\%*\%(x\%v\%"Smoking"))}


AMME(micro_model=SAOM_model,
     macro_model=FE_model,
     micro_process="reciprocity",
     mediator="indegree",
     macro_function=function(x){degree(x,cmode="indegree")},
     link_id=list(node_level_data$unique_id,node_level_data$unique_id,
                      node_level_data$unique_id,node_level_data$wave),
     object_type="network",
     controls=c("outdegree","AR1","gcc"),
     control_functions=list(function(x){degree(x,cmode="outdegree")},ar_function,gtrans),
     interval=c(0,.1),
     nsim=500,
     algorithm="parametric",
     silent=FALSE,
     SAOM_data = SAOM.Data,
     SAOM_var=list(Smoking=Smoking,Alcohol=Alcohol)) #provide var_list






################################
# AMME with pooled ERGM and SAOM
################################



#pooled ERGM
  #fit two ERGMs to two networks
data("faux.mesa.high")
model1<-ergm(faux.mesa.high~edges+
               nodecov("Grade")+
               nodefactor("Race")+
               nodefactor("Sex")+
               nodematch("Race")+
               nodematch("Sex")+
               absdiff("Grade"),
             control = control.ergm(parallel=4))

data("faux.magnolia.high")
model2<-ergm(faux.magnolia.high~edges+
               nodecov("Grade")+
               nodefactor("Race")+
               nodefactor("Sex")+
               nodematch("Race")+
               nodematch("Sex")+
               absdiff("Grade"),
             control = control.ergm(parallel=4))


#create node level data
node_level_data<-data.frame(grade=faux.mesa.high\%v\%"Grade",
                            sex=faux.mesa.high\%v\%"Sex",
                            degree=degree(faux.mesa.high),
                            betweenness=betweenness(faux.mesa.high),
                            gcc=gtrans(faux.mesa.high),
                            net_id="Mesa")

node_level_data$senior<-0
node_level_data$senior[node_level_data$grade==max(node_level_data$grade)]<-1
node_level_data$v_id<-1:network.size(faux.mesa.high)


node_level_data2<-data.frame(grade=faux.magnolia.high\%v\%"Grade",
                            sex=faux.magnolia.high\%v\%"Sex",
                            degree=degree(faux.magnolia.high),
                            betweenness=betweenness(faux.magnolia.high),
                            gcc=gtrans(faux.magnolia.high),
                            net_id="Magnolia")

node_level_data2$senior<-0
node_level_data2$senior[node_level_data$grade==max(node_level_data2$grade)]<-1
node_level_data2$v_id<-206:(network.size(faux.magnolia.high)+205)
node_level_data<-rbind(node_level_data,node_level_data2)


#estimate glm macro model with an AR 1 process
probit_model<-glm(senior~sex+degree+betweenness+gcc,
                data=node_level_data,
                family=binomial(link="probit"))



AMME(micro_model=list(model1,model2),
     macro_model=probit_model,
     micro_process="nodematch.Sex",
     mediator="degree",
     macro_function=degree,
     link_id=list(node_level_data$v_id,node_level_data$v_id,node_level_data$net_id),
     object_type="network",
     controls=c("betweenness","gcc"),
     control_functions=list(betweenness,gtrans),
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE)



##pooled SAOM with control functions using time varying covariates

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


##estimate the SAOM
SAOM_model<-siena07(create.model,
                    data=SAOM.Data,
                    effects=SAOM.terms,
                    verbose=TRUE)


##create node-level data
node_level_data<-data.frame(smoking=s50s[,1], #smoking behavior for DV
                            alcohol=s50a[,1],
                            v_id=rownames(s501), #unique node IDS
                            wave="Wave 1",       #unique time IDS
                            outdegree=rowSums(s501),
                            indegree=colSums(s501),
                            AR1=s501\%*\%s50s[,1],  #assign network autocorrelation
                            gcc=gtrans(as.network(s501)))

node_level_data<-rbind(node_level_data,data.frame(smoking=s50s[,2],
                                                  alcohol=s50a[,2],
                                                  v_id=rownames(s502),
                                                  wave="Wave 2",
                                                  outdegree=rowSums(s502),
                                                  indegree=colSums(s502),
                                                  AR1=s502\%*\%s50s[,2],
                                                  gcc=gtrans(as.network(s502))))



node_level_data<-rbind(node_level_data,data.frame(smoking=s50s[,3],
                                                  alcohol=s50a[,3],
                                                  v_id=rownames(s503),
                                                  wave="Wave 3",
                                                  outdegree=rowSums(s503),
                                                  indegree=colSums(s503),
                                                  AR1=s503\%*\%s50s[,3],
                                                  gcc=gtrans(as.network(s503))))


#recycle the same model for illustrative purposes
node_level_data$net_ID<-"Model 1"
node_level_data<-rbind(node_level_data,node_level_data)
node_level_data$net_ID[151:300]<-"Model 2"

##create unique identifiers for node-panel
  #ID for node-panel-model
node_level_data$unique_id<-paste(node_level_data$v_id,node_level_data$wave,node_level_data$net_ID)
  #ID for panel-model
node_level_data$unique_waves<-paste(node_level_data$wave,node_level_data$net_ID)

#estimate a linear network autocorrelation model with node fixed effects
FE_model<-lm(smoking~alcohol+outdegree+indegree+AR1+gcc+v_id,
              data=node_level_data)



##create user function calculate AR1 process on time varying node attributes
ar_function<-function(x){return(as.sociomatrix(x)\%*\%(x\%v\%"Smoking"))}

##estimate AMME
AMME(micro_model=list(SAOM_model,SAOM_model), #provide list of sienaFit objects
     macro_model=FE_model,
     micro_process="reciprocity",
     mediator="indegree",
     macro_function=function(x){degree(x,cmode="indegree")},
     link_id=list(node_level_data$unique_id,node_level_data$unique_id,
                        node_level_data$unique_id,node_level_data$unique_waves),
     object_type="network",
     controls=c("outdegree","AR1","gcc"),
     control_functions=list(function(x){degree(x,cmode="outdegree")},ar_function,gtrans),
     interval=c(0,.1),
     nsim=100,                  #parametric estimation requires more simulations than coefficients
     algorithm="parametric",
     silent=FALSE,
     SAOM_data = list(SAOM.Data,SAOM.Data), #list of siena objects
     SAOM_var=list(list(Smoking=Smoking,Alcohol=Alcohol),#provide var_list
                  list(Smoking=Smoking,Alcohol=Alcohol)))





#################################
# AMME with nested data
################################

####create dyad-level data

library(lme4)
library(btergm)
##use small data to simplify estimation
glm_dat<-edgeprob(model1)
glm_dat$net_id<-"mesa"
glm_dat2<-edgeprob(model2)
glm_dat2$net_id<-"magnolia"
glm_dat<-rbind(glm_dat,glm_dat2[,-c(4)])


##estimate micro model as glm for btoh networks using pooled ERGM data
net_glm<-glm(tie~nodecov.Grade+
                 nodefactor.Race.Hisp+
                 nodefactor.Race.NatAm+
                 nodefactor.Race.Other+
                 nodefactor.Sex.M+
                 nodematch.Race+
                 nodematch.Sex+
                 absdiff.Grade,
               data=glm_dat)


#create macro data
node_level_data<-data.frame(grade=faux.mesa.high\%v\%"Grade",
                            sex=faux.mesa.high\%v\%"Sex",
                            degree=degree(faux.mesa.high),
                            betweenness=betweenness(faux.mesa.high),
                            gcc=gtrans(faux.mesa.high),
                            net_id="Mesa")

node_level_data$senior<-0
node_level_data$senior[node_level_data$grade==max(node_level_data$grade)]<-1
node_level_data$v_id<-1:network.size(faux.mesa.high)


node_level_data2<-data.frame(grade=faux.magnolia.high\%v\%"Grade",
                             sex=faux.magnolia.high\%v\%"Sex",
                             degree=degree(faux.magnolia.high),
                             betweenness=betweenness(faux.magnolia.high),
                             gcc=gtrans(faux.magnolia.high),
                             net_id="Magnolia")

node_level_data2$senior<-0
node_level_data2$senior[node_level_data$grade==max(node_level_data2$grade)]<-1
node_level_data2$v_id<-206:(network.size(faux.magnolia.high)+205)
node_level_data<-rbind(node_level_data,node_level_data2)


#estimate glm macro model
probit_model<-glm(senior~sex+degree+betweenness+gcc,
                  data=node_level_data,
                  family=binomial(link="probit"))




AMME(micro_model=net_glm,
     macro_model=probit_model,
     micro_process="nodematch.Sex",
     mediator="degree",
     macro_function=degree,
     link_id=list(node_level_data$v_id,node_level_data$v_id,node_level_data$net_id),
     object_type="network",
     controls=c("betweenness","gcc"),
     control_functions=list(betweenness,gtrans),
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE,
     group_id=glm_dat$net_id,
     node_numbers = c(network.size(faux.mesa.high),
                      network.size(faux.magnolia.high)))




###using glmer for micro model
net_glmer<-glmer(tie~nodecov.Grade+
                 nodefactor.Race.Hisp+
                 nodefactor.Race.NatAm+
                 nodefactor.Race.Other+
                 nodefactor.Sex.M+
                 nodematch.Race+
                 nodematch.Sex+
                 absdiff.Grade+
                 (1|net_id),
               data=glm_dat)

probit_glmer<-glm(senior~sex+degree+betweenness+gcc,
                data=node_level_data,
                family=binomial(link="probit"))


AMME(micro_model=net_glm,
     macro_model=probit_glmer,
     micro_process="nodematch.Sex",
     mediator="degree",
     macro_function=degree,
     link_id=list(node_level_data$v_id,node_level_data$v_id,node_level_data$net_id),
     object_type="network",
     controls=c("betweenness","gcc"),
     control_functions=list(betweenness,gtrans),
     interval=c(0,1),
     nsim=50,
     algorithm="parametric",
     silent=FALSE,
     group_id=glm_dat$net_id,
     node_numbers = c(network.size(faux.mesa.high),
                      network.size(faux.magnolia.high)))






}

}%close examples

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~mediation }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~networks }% __ONLY ONE__ keyword per line
\keyword{ ~macro }% __ONLY ONE__ keyword per line
\keyword{ ~micro }% __ONLY ONE__ keyword per line
