#Matt Boyas, Michelle Newcomer, Ying Chao Shi
#Code Version 12/08/2013
#Known Bugs: None

set.seed(0)
#only set the seed if you are trying to replicate algorithmic output

#------Begin User input-----------------------------

#N(0,1)
density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))
xabs <- c(-1,1) #user-specified starting abcissae
accept <- 1000 #total number of points required to accept
ARS_values <- ARS(density, xabs, accept) #run this line to run the entire sampling scheme

#N(2,2)
density1 <- quote((1/sqrt(2*pi*4))*exp(-(x-2)^2/(2*4)))
xabs1 <- c(1,3)
accept <- 1000 
ARS_values1 <- ARS(density1,xabs1,accept)

#exp(0.5)
density2 <- quote(0.5*exp(-0.5*x))
xabs2 <- c(1,3)
accept <- 1000
ARS_values2 <- ARS(density2,xabs2,accept)

#Gamma(5,1)
density3 <- quote(1/24*x^4*exp(-x))
xabs3 <- c(2,6)
accept <- 1000
endpoints3 <- c(0,Inf)
ARS_values3 <- ARS(density3,xabs3,accept,endpoints3)

#Beta(2,2)
density4 <- quote(6*x*(1-x))
xabs4 <- c(0.3,0.7)
accept <- 1000 
endpoints4 <- c(0,1)
ARS_values4 <- ARS(density4,xabs4,accept,endpoints4)

#Pareto(2,3)
density5 <- quote((3*2^3)/x^(3+1))
xabs5 <- c(2, 5)
accept <- 1000
endpoints5 <- c(1,Inf)
ARS_values5 <- ARS(density5,xabs5,accept,endpoints5)

#-----End user Input--------------------------------



#------Begin ARS Scheme-----------------------------

ARS <- function(density, xabs, accept, endpoints=c(-Inf, Inf)){
  
  #--------Check to make sure abcissae are within the density endpoints------
  if (xabs[1]<endpoints[1] | xabs[2]>endpoints[2]){
    return("Initial abcissae are not in the density's specified domain.  Terminating ARS algorithm.")
  }
  #--------End check to make sure abcissae are within the density endpoints------
  
  
  
  #------Prep Function to define h(x) and h'(x)-----------------------------
  
  prep <- function(f){
    log_quote <- bquote(log(.(density)))
    deriv <- D(log_quote, "x")
    return(list(log_quote, deriv))
  }
  #compute the log and the derivative of the log in a preparation function
  
  prep <- prep(density)
  
  h <- function(x){
    return(eval(prep[[1]]))
  }
  
  h_prime <- function(x){
    return(eval(prep[[2]]))
  }
  
  #h and h_prime take x as the sole input and can be used in the computation of the equations in the paper algorithm
  #------End Prep Function to define h(x) and h'(x)-----------------------------
  

  #------Begin Functions for z and lower bound-----------------------------
  #Vectorized equation 1
  
  z <- c()
  zrange <- function(xabs){
    k <- length(xabs)
    j <- 1:(k-1)
    z[j] <- (h(xabs[j+1])-h(xabs[j])-xabs[j+1]*h_prime(xabs[j+1])+xabs[j]*h_prime(xabs[j]))/(h_prime(xabs[j])-h_prime(xabs[j+1]))
    return(z)
  }
  
  #lower bound
  lk_x <- function(x){
    j=1:(k-1)
    exp(((xabs[j+1]-x)*h(xabs[j])+(x-x[j])*h(xabs[j+1]))/(xabs[j+1]-xabs[j]))
  }
  #------End Functions for z and lower bound-----------------------------
  
  #------Begin Sampling-----------------------------
  #For exponential distributions, the upper bound is the same line as the original function. Thus, the inverse cdf method can be applied directly.
  if(h_prime(xabs[1])==h_prime(xabs[2])){
	u=runif(accept)
      sample=log(1-u)/h_prime(xabs[1])
  }else{

 	 sample <- c()
 	 while(length(sample)<accept){
 	   k=length(xabs) #the number of absiccae
 	   
 	   #The end points of the upper bound.
 	   zvalue <- zrange(xabs)
 	   interval <- c(endpoints[1], zvalue, endpoints[2]) #vector of the intervals for x
 	   segment <- c(interval,xabs) #vector of the all the points
 	   segment <- sort(segment)#sorted
 	   
 	   #The probabilities under each piece.
 	   prob <- c()
 	   for (i in 1:k){
 	     uk_x <- function(x){
 	       exp(h(xabs[i])+(x-xabs[i])*h_prime(xabs[i]))
 	     }

 	     prob[i] <- integrate(uk_x,interval[i],interval[i+1])$value
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
   	 #Second, integrate from the nth z to x*, set this integral to the remaining probablity, and solve x*. 
   	 x_star=log(sum(prob)*h_prime(xabs[n])*remain_prob/exp(h(xabs[n])-h_prime(xabs[n])*xabs[n])+exp(h_prime(xabs[n])*segment[2*n-1]))/h_prime(xabs[n])
   	 
   	 ###Now that we have the x*, just calculate the upper and lower bound.
   	 upper=exp(h(xabs[n])+(x_star-xabs[n])*h_prime(xabs[n]))
   	 #Remember we discussed that x* can fall into the two "tails". The paper saids if it does, then lower bound is -Inf.
   	 lower_segment=sort(c(-Inf,xabs,Inf))
   	 if(x_star<=xabs[1] || x_star>=xabs[k]){
   	   lower=-Inf
   	 } else {
		m=sum(x_star>xabs) 
	      lower=exp(((xabs[m+1]-x_star)*h(xabs[m])+(x_star-xabs[m])*h(xabs[m+1]))/(xabs[m+1]-xabs[m]))
	    }
      
      #Check for log concavity before moving forward
      if (upper<exp(h(x_star)) | lower>exp(x_star)){
 	     return("The density function is not logarithmically concave.  Terminating ARS algorithm.")}
      
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
  }
  	#------End Sampling-----------------------------
return(sample)
}
#------End ARS Scheme-----------------------------

#-----Begin Testing Scheme--------------------------
ARS_testing <- function(alpha=0.001){
 
  #Standard Normal Tests
  density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))
  xabs <- c(-1,1) #user-specified starting abcissae
  accept <- 1000 #total number of points required to accept
  sample <- ARS(density, xabs, accept)
  
  print("Beginning tests with Standard Normal")
  
  ##user visual plot for testing
  par(mfrow=c(2,1))
  #True density
  x <- seq(min(sample),max(sample),by=0.01)
  plot(x,eval(density),type="l",xlab="",ylab="",main="True Density")
  #Histogram of a sample from our function
  hist(sample, freq=F, main="Histogram of ARS sample", xlab="",ylab="")
  input <- readline("Press <y> followed by <return> if these plots appear to match.\nPress any other key followed by <return> to stop the testing script.")
  if (length(grep("y", input, ignore.case=TRUE))==0) {
    print("User determined that plots do not match. Terminating testing algorithm.")
    stop(call. = FALSE)
  }
  
  print("User determined that plots do match")
  print("Beginning Kolmogorov-Smirnov test")

  ###Kolmogorov-Smirnov test for equality of two cdf distributions
  KS <- ks.test(sample, "pnorm", 3, 2) # two-sided, exact
  par(mfrow=c(1,1))
  plot(ecdf(sample),lty=1, col="red", lwd=0.5, main="CDF Comparison Plot")
  lines(x, pnorm(x, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),lty=2,col="blue", lwd=1.25)
  legend("bottomright", legend=c("Empirical CDF of ARS Sample", "Standard Normal CDF"), cex=0.5, col=c("red","blue"), lty=c(1, 2), lwd=c(0.5, 1.25))
  ###Report the pvalue
  if(KS[[2]]<alpha){
    print(paste("P-value = ", KS[[2]]," < ", alpha, ". K-S test passed.", sep=""))
  } else{print(paste("P-value = ", KS[[2]]," !< ", alpha, ". P-value does not indicate statistical equality between the distributions. K-S test failed.", sep=""))}

  ###Pass a non-log concave density Pareto(2,3)
  print("Beginning test for non-log concave density.")
  density <- quote((7*1^7)/x^(7+1))
  xabs <- c(2, 5)
  endpoints <- c(1,Inf)
  sample <- ARS(density, xabs, accept, endpoints)
  if (sample!="The density function is not logarithmically concave.  Terminating ARS algorithm."){
    print("ARS algorithm did not detect non-log concave Pareto(2,3) density.  Test failed.  Terminating testing algorithm.")
    stop(call. = FALSE)
  }
  print("ARS algorithm detected non-log concave Pareto(2,3) density.  Non-log concave test passed."  )
  
  ##Pass a Beta(2,2) density a first starting abcissae<0
  print("Beginning test for an initial abcessa outside of defined density domain.")
  density <- quote(6*x*(1-x))
  xabs <- c(-2,0.7)
  accept <- 1000 
  endpoints <- c(0,1)
  sample <- ARS(density,xabs,accept,endpoints)
  if (sample!="Initial abcissae are not in the density's specified domain.  Terminating ARS algorithm."){
    print("ARS algorithm did not detect that an initial abcissa is outside of the defined density domain.  Test failed.  Terminating testing algorithm.")
    stop(call. = FALSE)
  }
  print("ARS algorithm detected initial abcessa outside of domain.  Initial abcessa domain test passed."  )

}
#-----End Testing Scheme-----------------------------


#------Run Testing Scheme------------------------
ARS_testing()
#------Complete Run of Testing Scheme------------------------
