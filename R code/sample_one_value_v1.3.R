#Matt Boyas, Michelle Newcomer, Ying Chao Shi
#Code Version 11/29/2013
#Known Bugs: None


set.seed(0)

#------Begin User input-----------------------------
density <- quote((1/sqrt(2*pi))*exp((-x^2)/2)) 
xabs<-c(-1,1) #user specified boundaries of domain
k=length(xabs) #starting value for the number of absiccae
accept<-300 #total number of points required to accept
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



#Vectorized equation 1

z <- c() 
zrange <- function(xabs){
  k <- length(xabs)
  j  <- 1:(k-1)
  z[j] <- (h(xabs[j+1])-h(xabs[j])-xabs[j+1]*hPrime(xabs[j+1])+xabs[j]*hPrime(xabs[j]))/(hPrime(xabs[j])-hPrime(xabs[j+1]))
  return(z)
}

#lower bound 
lk_x <- function(x){
  j=1:(k-1)
  exp(((xabs[j+1]-x)*h(xabs[j])+(x-x[j])*h(xabs[j+1]))/(xabs[j+1]-xabs[j]))
}

update<-c() #this vector should not change
sample <- c()

while (length(sample)<accept){
#Restart here when k=3,4,5
#initializes vectors and matrices
useg <- matrix(nrow=1000, ncol=k)
prob <- c()
x_star <- c()
upper <- c()
lower <- matrix(nrow=(k-1), ncol=(k-1))
accept_indicator <- 0
k <- length(xabs)

#Get z, which is 0 in this case.
zvalue <- zrange(xabs)

interval<-c(-Inf,zvalue,Inf) #vector of the intervals for x
segment<-c(interval,xabs) #vector of the all the points
segment<-sort(segment)#sorted to make the definition of useg easier



for (i in 1:k){
  useg[,i] <- seq(segment[i+1],segment[i+2],length=1000) #find the x values of the line segments for equation 2                 
  uk_x<-function(x){
    exp(h(xabs[i])+(x-xabs[i])*hPrime(xabs[i]))
  }  
  prob[i]<-integrate(uk_x,interval[i],interval[i+1])$value
}
  
while(accept_indicator==0){
  piece_num=sample(1:k,1,prob=prob/sum(prob))
  x_star<-sample(useg[,piece_num],1)
  upper<-uk_x(x_star) #this is already exp
  lower<-lk_x(x_star)#this is already exp
  
  #Generate w from uniform distribution.
  w=runif(1)
  
  #exp(a+b)=exp(a)*exp(b)
 # if(x_star>=xabs & x_star<=xabs[k]){ what is this IF step for?
  if(w<=(1/upper)*lower){
      accept_indicator <- 1
      sample <- c(sample, x_star)
      #xabs[k+1] <- x_star Are we supposed to update it in this situation?
      #xabs<-sort(xabs)
    }
  else if(w<=exp(h(x_star))*(1/upper)){
      accept_indicator <- 1
      sample <- c(sample, x_star)
      xabs[k+1] <- x_star     
      xabs <- sort(xabs)
    }
  #  else
     # sample=NA
 # }
}
}








