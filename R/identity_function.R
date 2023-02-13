#'identity_function is a helper function used internally to estimate the direct effect of m
#
#'takes networks object as its input
#
#'@x is the network object
#'@env calls the parent environment to access objects internal to AMME function

identity_function<-function(x){

  if(class(x)[1]=="igraph"){
    x<-intergraph::asNetwork(x)
  }
  macro_name<-sys.frame(sys.parent())$controls[sys.frame(sys.parent())$control]
  macro_process<-as.numeric(as.factor(x%v%macro_name))
  return(macro_process*sys.frame(sys.parent())$interval[sys.frame(sys.parent())$i])

}
