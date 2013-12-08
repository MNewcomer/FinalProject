#Matt Boyas, Michelle Newcomer, Ying Chao Shi
#Code Version 12/08/2013
#Known Bugs: None

set.seed(0)
#only set the seed if you are trying to replicate algorithmic output

#------Begin User input-----------------------------
density <- quote((1/sqrt(2*pi))*exp((-x^2)/2))
xabs <- c(-1,1) #user-specified starting abcissae
accept <- 5000 #total number of points required to accept

#Make sure to initialize the function "ARS" and "ARS_testing" before running line 12 
ARS_values <- ARS(density, xabs, accept) #run this line to run the entire sampling scheme
ARS_testing(ARS_values)
#-----End user Input--------------------------------



#------Begin ARS Scheme-----------------------------

ARS <- function(density, xabs, accept, endpoints=c(-Inf, Inf)){
  
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
    #Second, integrate from the nth z to x*, set this integral to the remaining probablity, and solve x*. (Don't spend time reading this mess. We can go over it together when we meet :)
    x_star=log(sum(prob)*h_prime(xabs[n])*remain_prob/exp(h(xabs[n])-h_prime(xabs[n])*xabs[n])+exp(h_prime(xabs[n])*segment[2*n-1]))/h_prime(xabs[n])
    
    ###Now that we have the x*, just calculate the upper and lower bound.
    upper=exp(h(xabs[n])+(x_star-xabs[n])*h_prime(xabs[n]))
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
    
    ###Test for concavity before moving on to the next sample round
    #Send xabs to the the h_prime function
    concave <- h_prime(xabs)
    runs <- (length(rle(sign(concave))$lengths)) # for concavity there should only be two runs and one switch in the runs
    if(runs>2){
      print("The log density function is not concave.")
      break}
    
  }
  #------End Sampling-----------------------------
return(sample)
}
#------End ARS Scheme-----------------------------

#-----Begin Testing Scheme--------------------------
ARS_testing <- function(sample){
  #------Begin Result-----------------------------
  
  ##user visual plot for testing
  
  par(mfrow=c(2,1))
  #True standard normal
  x=seq(-4,4,by=0.01)
  #Histogram of a sample from rnorm() with N(0,1) pdf lines
  hist(rnorm(5000), freq=F, xlim=c(-4, 4), ylim=c(0, 0.4), main="Histogram of rnorm(5000) \n with N(0,1) PDF superimposed", xlab="")
  lines(x, dnorm(x))
  #Histogram of a sample from our function
  hist(sample, freq=F, xlim=c(-4, 4), ylim=c(0, 0.4), main="Histogram of ARS sample", xlab="")
  input <- readline("Press <y> followed by <return> if these plots appear to match.\nPress any other key followed by <return> to stop the testing script.")
  if (length(grep("y", input, ignore.case=TRUE))==0) {
    print("User determined that plots do not match.  Terminating testing algorithm.")
    stop(call. = FALSE)
  }
  else{
    print("User determined that plots do match")
    print("Beginning Kolmogorov-Smirnov test on Standard Normal")
  
    ###Kolmogorov-Smirnov test for equality of two cdf distributions
    KS <- ks.test(sample, "pnorm", 3, 2) # two-sided, exact
    par(mfrow=c(1,1))
    plot(ecdf(sample),lty=1, col="red", lwd=0.5, main="CDF Comparison Plot")
    lines(x, pnorm(x, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),lty=2,col="blue", lwd=1.25)
    legend("bottomright", legend=c("Empirical CDF of ARS Sample", "Standard Normal CDF"), cex=0.3, col=c("red","blue"), lty=c(1, 2), lwd=c(0.5, 1.25))
    ###Report the pvalue
    if(KS[[2]]<0.001){print("P-value < 0.001. K-S test passed.")}
    else{print("P-value does not indicate statistical equality between the distributions")}
    print(KS)
  #------End Result-----------------------------
}
}
#-----End Testing Scheme-----------------------------



