#Matt Boyas, Michelle Newcomer, Ying Chao Shi
#Code Version 12/12/2013
#Known Bugs: None



#------Begin ARS Scheme-----------------------------

ars <- function(density, xabs, accept, endpoints=c(-Inf, Inf)){
  
  xabs <- sort(xabs)
  #sort initial abcessae input to make sure the user inputs them correctly
  
  #--------Check to make sure abcissae are within the density endpoints------
  if (xabs[1]<endpoints[1] | xabs[2]>endpoints[2]){
    stop("Initial abcissae are not in the density's specified domain. Terminating ARS algorithm.")
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
  #z  
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
  } else {

	  #Check to make sure the first derivatives of the two abcissae have opposite signs.
	  if(sign(h_prime(xabs[1]))*sign(h_prime(xabs[2]))>0){
		stop("Initial abcissae are not valid - their first derivatives have the same sign. Terminating ARS algorithm.")
	  }

	  #Sample the x* following the steps in the paper.
	  sample <- c()
	  while(length(sample)<accept){
	  	k=length(xabs) #the number of absiccae
      
		zvalue <- zrange(xabs)  #z
	  	interval <- c(endpoints[1], zvalue, endpoints[2]) #z and endpoints
	  	segment <- c(interval,xabs) #z, endpoints and tangent points
	  	segment <- sort(segment) #sorted
       
	  	#The probabilities under each piece.
	  	prob <- c()
	  	for (i in 1:k){
			uk_x <- function(x){
	      	exp(h(xabs[i])+(x-xabs[i])*h_prime(xabs[i]))
	      	}
			prob[i] <- integrate(uk_x,interval[i],interval[i+1])$value
	  	}
	         
	  	#Scale these probabilities so that they add up to 1.
	  	prob_each_piece=prob/sum(prob)
	           
	      #Generate a number from the uniform distribution to be our cdf.
	  	u=runif(1)
	
	  	#Check which piece we fall into. 
	  	prob_cumulative=cumsum(prob_each_piece)
	  	n=sum(u>prob_cumulative)+1 #The nth piece is chosen.
	
		#Subtract the probabilites of the previous n-1 pieces. Skip if the first piece was chosen. 
	  	if(n==1){
	  		remain_prob=u
	  	} else {
	  		remain_prob=u-prob_cumulative[n-1]
	  	}
	
		#Integrate from the nth z to x*, set this integral to be the remaining probablity, and solve for x*.
	  	x_star=log(sum(prob)*h_prime(xabs[n])*remain_prob/exp(h(xabs[n])-h_prime(xabs[n])*xabs[n])+exp(h_prime(xabs[n])*segment[2*n-1]))/h_prime(xabs[n])
	           
	  	#Calculate the upper and lower bound for this x*. Lower bound is -Inf if x* falls into the two "tails".
	  	upper=exp(h(xabs[n])+(x_star-xabs[n])*h_prime(xabs[n]))
		lower_segment=sort(c(-Inf,xabs,Inf))
		if(x_star<=xabs[1] || x_star>=xabs[k]){
	  		lower=-Inf
	 	} else {
	 		m=sum(x_star>xabs)
	      	lower=exp(((xabs[m+1]-x_star)*h(xabs[m])+(x_star-xabs[m])*h(xabs[m+1]))/(xabs[m+1]-xabs[m]))
	  	}
      
  		#Check for log concavity before moving forward
	  	if (upper<exp(h(x_star)) | lower>exp(x_star)){
	  		stop("The density function is not logarithmically concave. Terminating ARS algorithm.")
	  	}
	      
	  	#Either put x* into the sample or add it as another tangent point.
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
test <- function(density1=quote((1/sqrt(2*pi))*exp((-x^2)/2)), xabs1=c(-1,1), accept1=1000, endpoints1=c(-Inf, Inf), density2=quote(1/24*x^4*exp(-x)), xabs2=c(2, 6), accept2=1000, endpoints2=c(0,Inf), true_dens=quote(rgamma(1000, 5, 1)), alpha=0.01, density3=quote((7*1^7)/x^(7+1)), xabs3=c(2, 5), accept3=1000, endpoints3=c(1, Inf), density4=quote(6*x*(1-x)), xabs4=c(-2,0.7), accept4=1000, endpoints4=c(0, 1)){
  #Distribution Tests
  sample <- invisible(ars(density1, xabs1, accept1, endpoints1))
  
  print("Beginning distribution similarity test")
  
  ##user visual plot for testing
  x <- seq(min(sample),max(sample),by=0.01)
  hist(sample, freq=F, main="Histogram of ARS sample with true density superimposed", xlab="",ylab="")
  lines(x,eval(density1))
  input <- readline("Press <y> followed by <return> if these plots appear to match.\nPress any other key followed by <return> to stop the testing script.")
  if (length(grep("y", input, ignore.case=TRUE))==0) {
    warning("User determined that plots do not match. Test failed.")
  }
  
  print("User determined that plots do match")
  print("Beginning Kolmogorov-Smirnov test")

  ###Kolmogorov-Smirnov test for equality of two cdf distributions
  sample <- invisible(ars(density2, xabs2, accept2, endpoints2))
  trueSample <- eval(true_dens)
  KS <- ks.test(sample, trueSample, alternative="two.sided", exact=TRUE) # two-sided, exact
  par(mfrow=c(1,1))
  plot(ecdf(sample), lty=1, col="red", lwd=0.5, main="CDF Comparison Plot")
  lines(ecdf(trueSample), lty=2, col="blue", lwd=1.25)
  legend("bottomright", legend=c("Empirical CDF of ARS Sample", "Empirical CDF of True Sample"), cex=0.5, col=c("red","blue"), lty=c(1, 2), lwd=c(0.5, 1.25))
  ###Report the pvalue
  if(KS[[2]]<alpha){
    warning(paste("P-value = ", KS[[2]]," < ", alpha, ". K-S test suggests that the two distributions are not equivalent. K-S test failed.", sep=""))
  } else{print(paste("P-value = ", KS[[2]]," >= ", alpha, ". P-value does not reject the null; evidence for equivalence of the two distributions. K-S test passed.", sep=""))}

  ###Pass a non-log concave density
  print("Beginning test for non-log concave density.")
  err <- try(ars(density3, xabs3, accept3, endpoints3), silent=TRUE)
  if (is(err, "try-error")==FALSE){
    warning("ARS algorithm did not detect non-log concave density. Test failed.")
  } else print("ARS algorithm detected non-log concave density. Non-log concave test passed." )
  
  ##Pass a density with starting abcissae<0
  print("Beginning test for an initial abcessa outside of defined density domain.")
  err <- try(ars(density4, xabs4, accept4, endpoints4), silent=TRUE)
  if (is(err, "try-error")==FALSE){
    warning("ARS algorithm did not detect that an initial abcissa is outside of the defined density domain. Test failed.")
  } else {print("ARS algorithm detected initial abcessa outside of domain. Initial abcessa domain test passed." )}
  
  print("Testing algorithm complete.")
}
#-----End Testing Scheme-----------------------------

#------Begin Example User input-----------------------------
#N(0,1)
density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))
xabs <- c(-1,1) #user-specified starting abcissae
accept <- 1000 #total number of points required to accept
ARS_values <- ars(density, xabs, accept) #run this line to run the entire sampling scheme

#N(2,2)
density1 <- quote((1/sqrt(2*pi*4))*exp(-(x-2)^2/(2*4)))
xabs1 <- c(1,3)
accept <- 1000
ARS_values1 <- ars(density1,xabs1,accept)

#exp(0.5)
density2 <- quote(0.5*exp(-0.5*x))
xabs2 <- c(1,3)
accept <- 1000
ARS_values2 <- ars(density2,xabs2,accept)

#Gamma(5,1)
density3 <- quote(1/24*x^4*exp(-x))
xabs3 <- c(2,6)
accept <- 1000
endpoints3 <- c(0,Inf)
ARS_values3 <- ars(density3,xabs3,accept,endpoints3)

#Beta(2,2)
density4 <- quote(6*x*(1-x))
xabs4 <- c(0.3,0.7)
accept <- 1000
endpoints4 <- c(0,1)
ARS_values4 <- ars(density4,xabs4,accept,endpoints4)

#Pareto(2,3)
density5 <- quote((3*2^3)/x^(3+1))
xabs5 <- c(2, 5)
accept <- 1000
endpoints5 <- c(1,Inf)
ARS_values5 <- ars(density5,xabs5,accept,endpoints5)

#-----End example user Input--------------------------------


#------Run Testing Scheme------------------------
test()
#------Complete Run of Testing Scheme------------------------