#Matt Boyas, Michelle Newcomer, Ying Chao Shi
#Code Version 12/01/2013
#Known Bugs: None

set.seed(0)

#------Begin User input-----------------------------
density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))
xabs<-c(-1,1) #user specified starting points
accept<-5000 #total number of points required to accept
#-----End user Input--------------------------------

#------Prep Function to define h(x) and h'(x)-----------------------------

prep <- function(f){
  logQuote <- bquote(log(.(density)))
  deriv<- D(logQuote, "x")
  return(list(logQuote, deriv))
}
#compute the log and the derivative of the log in a preparation function

prep <- prep(density)

h <- function(x){
  return(eval(prep[[1]]))
}

hPrime <- function(x){
  return(eval(prep[[2]]))
}

#h and hPrime take x as the sole input and can be used in the computation of the equations in the paper algorithm
#------End Prep Function to define h(x) and h'(x)-----------------------------


#------Begin Functions for z and lower bound-----------------------------
#Vectorized equation 1

z <- c()
zrange <- function(xabs){
  k <- length(xabs)
  j <- 1:(k-1)
  z[j] <- (h(xabs[j+1])-h(xabs[j])-xabs[j+1]*hPrime(xabs[j+1])+xabs[j]*hPrime(xabs[j]))/(hPrime(xabs[j])-hPrime(xabs[j+1]))
  return(z)
}

#lower bound
lk_x <- function(x){
  j=1:(k-1)
  exp(((xabs[j+1]-x)*h(xabs[j])+(x-x[j])*h(xabs[j+1]))/(xabs[j+1]-xabs[j]))
}
#------End Functions for z and lower bound-----------------------------


#------Begin Sampling-----------------------------
sample <- c()

while(length(sample)<accept){
	k=length(xabs) #the number of absiccae

	#The end points of the upper bound.
	zvalue <- zrange(xabs)
	interval<-c(-Inf,zvalue,Inf) #vector of the intervals for x
	segment<-c(interval,xabs) #vector of the all the points
	segment<-sort(segment)#sorted

	#The probabilities under each piece.
	prob <- c()
	for (i in 1:k){
		uk_x<-function(x){
 		exp(h(xabs[i])+(x-xabs[i])*hPrime(xabs[i]))
		}
  		prob[i]<-integrate(uk_x,interval[i],interval[i+1])$value
	}

	#Scale these probabilities so that they add up to 1 (because the upper bound is a pdf).
	prob_each_piece=prob/sum(prob)

	###Sample the x* using the inverse cdf method.
	##Generate a number from the uniform distribution, and this is our cdf.
	u=runif(1)
	##Find the corresponding x*.
	#First, check which piece we fall into. If we are in the first piece, integrate (the second step) directly. If we are in the nth piece (n>1), we have to subtract the probabilites of the previous n-1 pieces.
	prob_cumulative=cumsum(prob_each_piece)
	n=sum(u>prob_cumulative)+1 #The nth piece is chosen.
	if(n==1){
		remain_prob=u
	} else {
		remain_prob=u-prob_cumulative[n-1]
	}
	#Second, integrate from the nth z to x*, set this integral to the remaining probablity, and solve x*. (Don't spend time reading this mess. We can go over it together when we meet :)
	x_star=log(sum(prob)*hPrime(xabs[n])*remain_prob/exp(h(xabs[n])-hPrime(xabs[n])*xabs[n])+exp(hPrime(xabs[n])*segment[2*n-1]))/hPrime(xabs[n])

	###Now that we have the x*, just calculate the upper and lower bound.
	upper=exp(h(xabs[n])+(x_star-xabs[n])*hPrime(xabs[n]))
	#Remember we discussed that x* can fall into the two "tails". The paper saids if it does, then lower bound is -Inf.
	lower_segment=sort(c(-Inf,xabs,Inf))
	if(x_star<=xabs[1] || x_star>=xabs[k]){
		lower=-Inf
	} else {
		m=sum(u>xabs) 
		lower=exp(((xabs[m+1]-x_star)*h(xabs[m])+(x_star-xabs[m])*h(xabs[m+1]))/(xabs[m+1]-xabs[m]))
	}

	###Either put x* into the sample or add it as another tangent point. 
	#Generate w from uniform distribution.
	w=runif(1)
	if(w<=lower/upper){
		sample <- c(sample, x_star)
  	} else {
		if(w<=exp(h(x_star))/upper){
		sample <- c(sample, x_star)
		} 
     		xabs[k+1] <- x_star
      	xabs <- sort(xabs)
	}
}
#------End Sampling-----------------------------

#------Begin Result-----------------------------
sample
xabs

##ploting for testing
par(mfrow=c(2,2))
#True standard normal
x=seq(-3,3,by=0.01)
plot(x,dnorm(x),type="l")
#Histogram of a sample from rnorm()
hist(rnorm(5000),freq=F)
#Histogram of a sample from our function
hist(sample,freq=F)
#------End Result-----------------------------




