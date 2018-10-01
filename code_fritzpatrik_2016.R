#modeling the fritz paper
library(deSolve)
SIR<-function()
Lorenz<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dS <- r*S - *S*N - b*S*I
      dE <- b*S*I-(g+)
      dZ <- -X*Y + c*Y - Z
      
      # return the rate of change
      list(c(dX, dY, dZ))
    })   # end with(as.list ...
}

#using gillespie algorithm for a stochastic treatment of SIR model
library(GillespieSSA)
## Kermack-McKendrick SIR model
## Not run:
par(mfrow=c(1,1))
parms <- c(beta=0.001, gamma=0.1)
x0  <- c(S=499,I=1,R=0)
a   <- c("beta*S*I","gamma*I")
nu  <- matrix(c(-1,0,+1,-1,0,+1),nrow=3,byrow=TRUE)
out <- ssa(x0,a,nu,parms,tf=100,simName="SIR model")
ssa.plot(out)

#simple rabies model in anderson 1981 - equations from childs

#1. Defining model parameters and their values

parameters <- c(N=100, beta = 79, gamma = 0.5, sigma = 1/13, alpha = 1/73, r = 0.5, b = 0.5 )

#2. Defining model state variables and their initial conditions
state <- c(X = 1, I = 1 , Y=1)

#3. Implementing the model equations that calculate the rate of change of the
#state variables (here using the function and equations from Anderson 1981)
SIR<-function(t, state, parameters) {
     with(as.list(c(state, parameters)),{
       # rate of change
           dX <- r*X - gamma*X*N - beta*X*Y
           dI <- beta*X*Y - (sigma + b + gamma*N)*I
           dY <- sigma*I - (alpha + b + gamma*N)*Y
      
           # return the rate of change
           list(c(dX, dI, dY))
         })   # end with(as.list ...
}

#Model application 
times <- seq(0, 0.01, by = 0.00001)

#Model integration
#The  model  is  solved  using deSolve function ode,  which  is  the  default  integration  routine.
# Function ode takes as input, a.o.  the state variable vector (y), the times at which output is 
# required (times), the model function that returns the rate of change (func) and the parameter
# vector (parms)

library(deSolve)
out <- ode(y = state, times = times, func = SIR, parms = parameters)
head(out)

#plotting results
par(oma = c(0, 0, 3, 0))
head(out)
plot(out[, "X"], xlab = "time", ylab = "-", pch=".", ylim=range(c(0,2)))
points(out[, "I"], col="red", pch=".")
points(out[, "Y"], col="blue", pch=".")

mtext(outer = TRUE, side = 3, "SIR model", cex = 1.5)

