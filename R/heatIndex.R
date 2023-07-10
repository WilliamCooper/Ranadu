#' @title heat index
#' @description Calculates the heat index given the temperature and dewpoint.
#' @details This function duplicates the NWS approximating representation
#' described here: https://www.weather.gov/media/ffc/ta_htindx.PDF .
#' https://rpubs.com/cooperwilliama/797335 
#' If the known quantities are temperature and vapor pressure
#' you can use DPfromE() to find the dewpoint. The
#' input parameters may be vectors, all of the same size, and in that
#' case the result will also be a vector. In addition, the function
#' calculates an apparent temperature described in XXX that
#' can be used as an alternative to the heat index.
#' @aliases heatIndex, HeatIndex, heatindex
#' @author William Cooper
#' @export heatIndex
#' @param AT A numeric representing air temperature in deg. C 
#' @param DPT A numeric representing dewpoint temperature in deg. C. 
#' @param RTA A logical variable specifying if the apparent temperature
#' TA should be returned in addition to the heat index. Default is FALSE.
#' @importFrom stats nlm
#' @return Either a numeric representing the heat index in deg. C (if
#' RTA is FALSE) or a named list containing the heat index and the apparent
#' temperature in deg. C (if RTA is TRUE). Note that
#' these units are not the conventional ones for the heat index.
#' @examples # expect xx
#' HI <- Ranadu::heatIndex(20, 10) 
heatIndex <- function(AT, DPT, RTA=FALSE) {
    TF <- AT * 9 / 5 + 32
    E <- MurphyKoop(DPT)
    ES <- MurphyKoop(AT)
    RH <- 100 * E / ES 
    TF2 <- TF^2
    TFRH <- TF*RH
    RH2 <- RH^2
    # print (sprintf('DPT=%.1f, AT=%.1f, TF=%.1f, RH=%.1f', DPT, AT, TF, RH))
    HI = 0.5 * (TF + 61.0 + (TF - 68.0) * 1.2 + (RH * 0.094))
    CF <- c(-42.379, 2.04901523, 10.14333127, -0.22475541, -6.83783e-3,
           -5.481717e-2, 1.22874e-3, 8.5282e-4, -1.99e-6)
    if (HI > 80) {
        HI2 <-  -42.379 + 2.04901523 * TF + 10.14333127 * RH - 0.22475541 * TF * RH - (6.83783e-3) * TF * TF -
               5.481717e-2 * RH * RH + (1.22874e-3) * TF * TF * RH +
              8.5282e-4 * TF * RH * RH - 1.99e-6 * TF * TF * RH * RH
        HI <- CF[1] + CF[2] * TF + CF[3] * RH + CF[4] * TFRH + CF[5] * TF2 + CF[6] * RH2 +
               CF[7] * TF2 * RH + CF[8] * TF * RH2 + CF[9] * TF2 * RH2
        # print (sprintf ('HI=%.1f, HI2=%.1f', HI, HI2))
        if ((RH < 13) && (TF > 80) && (TF < 112)) {
            HI <- HI - (13 - RH) / 4 * sqrt((17 - abs(TF-95.)) / 17)
        }
        if ((RH > 85) && (TF > 80) && (TF < 87)) {
            HI <- HI + (RH - 85) / 10 * (87-TF) / 5
        }
    }
    if (HI < 40) {HI <- TF}
    HI <- (HI - 32) * 5 / 9  # Convert to degC
    if (RTA) {
      TWB <- wetbulbT(1000, AT, DPT)
      at <- function(A) {
        DPR <- DPfromE(MurphyKoop(A) * 0.20)
        if (DPR > A) {DPR <- A}
        return(abs(wetbulbT(1000, A, DPR) - TWB))
      }
      T1 <- (AT+HI)/1.5  # The first guess
      X <- nlm(at, T1, steptol=0.03) # X$estimate is the solution
      TA <- X$estimate
      return(list(HI=HI, TA=TA))
    } else {
        return(HI)
    }
}


