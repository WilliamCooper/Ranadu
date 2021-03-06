#' @title MachNumber
#' @description Calculates the Mach number 
#' @details returns the Mach number calculated from ambient and dynamic pressure and, 
#' optionally, corrected for water vapor pressure
#' @aliases machnumber
#' @author William Cooper
#' @export MachNumber
#' @param P A numeric (optionally vector) giving the ambient pressure in hPa 
#' @param Q A numeric (optionally vector) giving the dynamic pressure in hPa 
#' @param E A numeric (optionally vector) giving the water vapor pressure in hPa 
#' @return A numeric (possible vector) giving the MACH Number
#' @examples 
#' MACH <- MachNumber(500., 50.)
#' MACH <- MachNumber(700., 60., 1.1)
#' MACH <- MachNumber (RAFdata$PSXC, RAFdata$QCXC, RAFdata$EWX)
#' 
MachNumber <- function (P, Q, E=0.) {
# function to calculate the Mach number, optionally with
# humidity correction. Call with pressure, dynamic pressure,
# and vapor pressure, all in the same units.
# Any units can be used as long as the same are used for all
# three arguments. Only ratios (Q/P, E/P) enter the formula
  CP <- SpecificHeats (E/P)
  return (as.vector(((2.*CP[,2]/CP[,3])*(((P+Q)/P)**(CP[,3]/CP[,1])-1.))**0.5))
}

