# Filename: NSimulation.R
# All rights reserved.

# Calculate n. of Simulation 
#
# acc A character indicating the level of accuracy of the estimate of interest
# var A character indicating the variance for the parameter of interest
# alpha A character indicating significance level (default is 0.05)
# digits Rounding of Numbers (default is 0)
#  Returns a number containing the minimum number of simulation.

Nsimulation <- function(acc,var,alpha=0.05, digits=0) {
  z<-qnorm(1-alpha/2)
  B<- round((z^2*(var))/(acc)^2,digits)
  return(B)
}
