#All
#Code Version 11/24/2013

##h and h'

density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))

prep <- function(f){
  logQuote <- bquote(log(.(density)))
  deriv<- D(logQuote, "x")
  return(list(logQuote, deriv))
}

prep <- prep(density)

h <- function(x){
  return(eval(prep[[1]]))
}

hPrime <- function(x){
  return(eval(prep[[2]]))
}

##sampling one value from the standard normal.

#Set starting values as -1 and 1. 
xabs=c(-1,1)
k=length(xabs)

#Vectorized version to get z's.
z <- c() 
zrange <- function(xabs){
	j=1:(k-1)
	z[j] <- (h(xabs[j+1])-h(xabs[j])-xabs[j+1]*hPrime(xabs[j+1])+xabs[j]*hPrime(xabs[j]))/(hPrime(xabs[j])-hPrime(xabs[j+1]))
	return(z)
}

#Get z, which is 0 in this case.
zvalue <- zrange(xabs)

#Sample from upper bound.(Didn't figure out how to do every piece simultaneously)
uk_x1<- function(x){
	exp(h(xabs[1])+(x-xabs[1])*hPrime(xabs[1]))
}
prob1=integrate(uk_x1,-Inf,zvalue)$value

uk_x2<- function(x){
	exp(h(xabs[2])+(x-xabs[2])*hPrime(xabs[2]))
}
prob2=integrate(uk_x2, zvalue,Inf)$value

prob_sum=prob1+prob2

piece_num=sample(1:k,1,prob=c(prob1/prob_sum,prob2/prob_sum))
lambda=c(hPrime(xabs[piece_num]))
x=rexp(1,abs(lambda))

#Generate w from uniform distribution.
w=runif(1)

#lower bound and upper bound.
lk_x <- function(x){
	j=1:(k-1)
	exp(((xabs[j+1]-x)*h(xabs[j])+(x-x[j])*h(xabs[j+1]))/(xabs[j+1]-xabs[j]))
}
uk_x<- function(x){
	exp(h(xabs[piece_num])+(x-xabs[piece_num])*hPrime(xabs[piece_num]))
}

#Compare. The sampled points either go to the sample or the xabs.

sample=c()

if(x>=xabs[1] & x<=xabs[k]){
	if(u<lk_x(x)/uk_x(x))
	sample[1]=x
	else if(u<exp(h(x))/uk_x(x))
	sample[1]=x
	else
	xabs[k+1]=x
}

if(x<=xabs[1] & x>=xabs[k]){
	if(u<exp(h(x))/uk_x(x))
	sample[1]=x
	else
	xabs[k+1]=x
}



