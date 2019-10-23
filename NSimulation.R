# Filename: NSimulation.R
# Authors: Francesca Graziano & Paola Rebora
# All rights reserved

# Calculate number of simulation needed
#
#  @param acc A character indicating the level of accuracy of the estimate of interest
#  @param var A character indicating the variance for the parameter of interest
#  @param alpha A character indicating significance level (default is 0.05)
#  @param digits Rounding of Numbers (default is 0)

#  Returns a number containing the minimum number of simulation.

Nsimulation <- function(acc,var,alpha=0.05, digits=0) {
  z<-qnorm(1-alpha/2)
  B<- round((z^2*(var))/(acc)^2,digits)
  return(B)
}
