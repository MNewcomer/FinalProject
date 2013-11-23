#the functions will be designed to take the following inputs:
#the density, in the form of an R expression via expression or quote with variable x

#Matt Boyas, Michelle Newcomer
#Code Version 11/23/2013
#Known Bugs: None



#------Begin User input-----------------------------
density <- quote((1/sqrt(2*pi))*exp((-x^2)/2)) #Or could we have the user quote something like dnorm(mu, sigma)?
k<-6 #number of points to make the upper hull. Minimum is 2
d<-c(-3,3) #user specified boundaries of domain
#work with standard normal to begin with
#-----End user Input--------------------------------




#------Prep Function to define h(x) and h'(x)-----------------------------

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
#------End Prep Function to define h(x) and h'(x)-----------------------------





#------Define Functions required for initialization-----------------------------


#Calculate zj to determine the tip of the upper hull
z <- c() #initialize the vector that contains the z values

zrange <- function(k,xabs){
    for (j in 1:(k-1)){
      z[j] <- (h(xabs[j+1])-h(xabs[j])-xabs[j+1]*hPrime(xabs[j+1])+xabs[j]*hPrime(xabs[j]))/(hPrime(xabs[j])-hPrime(xabs[j+1])) #this is equation 1
    }
  return(z)
}

#calculate the equations of the tangent lines
xabs <- c() #initializes the vector that contains the abscissae

uk <- function(d,k,index){

  xabs <-  (d[2]-d[1])/(k-1)*(index-1)+d[1] #find the abscissae
  zvalue <- zrange(k,xabs) #find the z values in equation 1

  useg <- c(d[1],zvalue,d[2]) #find the x values of the line segments for equation 2
  xseq <- matrix(nrow=100, ncol=k) #initializes the matrix x values used for equation 2
  uk_x <- matrix(nrow=nrow(xseq),ncol=k) #initializes the matrix for each tangent line uk(x)
  
  for (j in 1:k){
    xseq[,j] <- seq(useg[j],useg[j+1],length=100) #x segment used to calculate equation 2
    uk_x[,j] <-h(xabs[j])+(xseq[,j]-xabs[j])*hPrime(xabs[j]) #equation 2
  }
  
  return(list(xabs,zvalue,useg,xseq,uk_x))
}
#--------------End of functions required for the initialization step----------------------





#------Begin initialization step. First evaluate if the derivative is (+) or (-) --------------------
initialize <- function(d,k){
  if (hPrime(d[1])>0 & hPrime(d[2])<0){
    index<-seq(1,k,by=1)  #creates a vector for the number of absiccae
    u <- uk(d,k,index) #references the function that calculates u in equation 2

  }
  else {
    break
  }
  return(u)
}
#------End initialization step evaluates if the derivative is (+) or (-) --------------------




#------Run the ARS scheme--------------------
upper <- initialize(d,k) #runs the scheme from the user inputs d and k




#------------------Plot results---------------

x <- seq(d[1],d[2],by=0.01)
plot(x,h(x),type="n",ylim=c(h(d[1]),max(upper[[5]])))
lines(x,h(x))
for (i in 1:k){

  lines(upper[[4]][,i],upper[[5]][,i],lty=2) #plots the tangent lines
  points(upper[[1]],h(upper[[1]])) #plots the absiccae
  #points(upper[[2]][j-1],tail(upper[[5]][,j-1],1))

}


