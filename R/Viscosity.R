#' @title Viscosity
#' @description Provides the dynamic and kinematic viscosity of dry air in SI units.
#' @details The Sutherland formula is used, with appropriate constants for dry air.
#' The result for dynamic viscosity is assumed independent of pressure. To
#' calculate kinematic viscosity, the air density is calculated for dry air
#' from the provided temperature and pressu8re.
#' @aliases Viscosity
#' @author William Cooper
#' @export Viscosity
#' @param T The air temperature in deg.C. Default: 0
#' @param p The air pressure in hPa. Default 1013.25.
#' @return A two-component list containing the dynamic and kinematic viscosity.
#' @examples 
#' Viscosity (15., 500)
Viscosity <- function (T=0, p=1013.25) {
  Cs <- 113 # Sutherland's constant for air (a compromise between N2 and O2) [K]
  Tr <- 273.15 # Reference temperature [K]
  nu_0 <- 0.00001716 # Reference viscosity [Pa s]
  Tk <- T + 273.15
  # Sutherland's formula:
  nu <- nu_0 * (Tk / Tr)^1.5 * (Cs + Tr) / (Cs + Tk)
  rho_a <- p * 100 / (StandardConstant('Rd') * (Tk)) 
  eta <- nu / rho_a
  return(c(nu, eta))
}