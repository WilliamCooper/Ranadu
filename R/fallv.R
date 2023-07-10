#' @title fallv
#' @description Provides the fall speeds of spherical hydrometeors
#' @aliases fallv
#' @author William Cooper
#' @export fallv
#' @param d The droplet diameter [m]; default=20e-6
#' @param rho Density of the hydrometeor. Default=1.e3
#' @param T temperature [deg.C]
#' @param p pressure [hPa]
#' @return The terminal fall speed [m/s]
#' @examples 
#' fallv (100e-6)
fallv <- function(d, rho=1.e3, T=0, p=1013.25) {
  # units: MKS exc. p is hPa, T is deg.C
  # Careful: d is expected in m
  # result is m/s
  # Results are not valid over 7 mm; will return 7 mm value
  visc <- Viscosity(T, p)
  Tk <- T + 273.15
  rho_a <- p * 100 / (StandardConstant('Rd') * (Tk)) 
  dv <- visc[1]
  sigma0 <- 75.7
  B <- c(-0.318657e1, 0.992696, -0.153193e-2,-0.987059e-3, -0.578878e-3,
         0.855176e-4, -0.327815e-5)
  C <- c(-0.500015e1, 0.523778e1, -0.204914e1, 0.475294, -0.542819e-1,
         0.238449e-2)
  A <- d / 2
  sigma <- (sigma0 - 0.154 * T) * 1.e-3
  # c.....There are three regimes, <19 um diameter, 19--1000 um diameter, and >1 mm diameter
  # c.....Separation is due to deformation of drop from spherical
  # c.....   at sizes larger than 1 mm and slip-flow corrections for small drops.
  # c.
  # c....... first, protect against beyond-range values
  A[A > 0.35e-2] <- 0.35e-2
  vt <- vector('numeric')
  for (a in A) {
    if (a < 9.5e-6) {
        # Stokes solution with slip correction
        elr <- 6.62e-8
        pr <- 101325
        Tr <- 293.15
        etar <- 0.00001818 #kg m-1 s-1
        el <- elr * dv / etar * pr / (p * 100) * ((T+273.15)/Tr)^0.5
        Csc <- 1+2.51*el/(2*a)
        rhoa <- p * 100 / (StandardConstant('Rd') * 273.15)
        fv <- (a*2)^2*(1000-rhoa)*Gravity(45)/(18*dv)*Csc
    }
    else if (a < 535e-6) {
      # c.....Regime I:  diameter smaller than 1 mm
      cdn2 <- 32. * a^3 * (rho - rho_a) * rho_a * Gravity(45) / (3. * dv^2)
      x <- log(cdn2)^(1:6)
      y <- B[1] + sum(B[2:7] * x)
      re <- exp(y)
      fv <- dv * re / (2. * rho_a * a)
    } else {
      # c.....Regime II:  diameter larger than 1000 um
      bo <- Gravity(45) * (rho - rho_a) * a^2 / sigma
      pp <- sigma^3 * rho_a^2 / (dv^4 * Gravity(45) * (rho-rho_a))
      x <- log(16. / 3. * bo * pp^0.166666) ^(1:5)
      y <- C[1] + sum(C[2:6] * x)
      re <- pp^0.166666 * exp(y)
      fv <- dv * re / (2. * rho_a * a)
    }
    vt <- c(vt, fv)
  }
  return(vt)
}
# diam <- 1.e-6 * 10^seq(0, 4.5, by=0.1)
# diam <- diam[diam < 6e-3]
# fvs <- rep(0, length(diam))
# for (i in 1:length(diam)) {fvs[i] <- fallv(diam[i])}
# plot(diam*1.e6, fvs, type='l', col='blue', lwd=2, log='xy', xlim=c(1, 1000), ylim=c(1.e-4, 10))


# function fallv(dd,rho,t,p)
# c.....Input:  dd=diameter, rho=particle density, t,p=air temp and pressure
# c.....All units are cgs except pressure, which enters in mb.
# C.....In particular, output fallv is cm/s
# save b,c,visc0,sigma0
# dimension b(7),c(6)
# data visc0/1.718E-4/,sigma0/75.7/
#   Data b/-0.318657e1,0.992696,-0.153193e-2,-0.987059e-3,
# $   -0.578878e-3,0.855176e-4,-0.327815e-5/
#   data c/-0.500015e1,0.523778e1,-0.204914e1,0.475294,
# $   -0.542819e-1,0.238449e-2/
#   a=dd/2.
# ta=t+273.15
# dens=p*1.e3*28.9644/(8.314e7*ta)
# visc=visc0+0.0049e-4*t
# if(t.lt.0.) visc=visc-1.2e-9*t**2
# sigma=sigma0-0.154*t
# c.....There are two regimes, separated at 1 mm diameter
# c.....Separation is due to deformation of drop from spherical
# c.....   At sizes larger than 1 mm
# c.
# c....... first, protect against beyond-range values
# if(a.gt.0.25) a = 0.25
# If(a.gt.500.e-4) go to 500
# c.....Regime I:  diameter smaller than 1 mm
# cdn2=32.*a**3*(rho-dens)*dens*980.665/(3.*visc**2)
# x=alog(cdn2)
# y=b(1)
# do 10 i=1,6
# y=y+b(i+1)*x**i
# 10    continue
# re=exp(y)
# fallv=visc*re/(2.*dens*a)
# return
# c  
# c.....Regime II:  diameter larger than 1000 um
# 500   continue
# bo=980.665*(rho-dens)*a**2/sigma
# pp=sigma**3*dens**2/(visc**4*980.665*(rho-dens))
# x=alog(16./3.*bo*pp**0.166666)
# y=c(1)
# do 510 i=1,5
# y=y+c(i+1)*x**i
# 510   continue
# re=pp**0.166666*exp(y)
# fallv=visc*re/(2.*dens*a)
# return
# end
