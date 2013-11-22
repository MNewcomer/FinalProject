#the functions will be designed to take the following inputs:
#the density, in the form of an R expression via expression or quote with variable x

#Matt Boyas
#Code Version 11/22/2013
#Known Bugs: None

density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))
#the above will be defined by the user
#work with standard normal to begin with

prep <- function(f){
logQuote <- bquote(log(.(density)))
deriv<- D(logQuote, "x")
return(list(logQuote, deriv))
}
#compute the log and the derivative of the log in a preparation function so we don't have to evaluate these expressions every time we want h or h'

prep <- prep(density)

h <- function(x){
  return(eval(prep[[1]]))
}

hPrime <- function(x){
  return(eval(prep[[2]]))
}
#h and hPrime take x as the sole input and can be used in the computation of the equations in the paper algorithm