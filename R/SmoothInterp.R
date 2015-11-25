#' @title SmoothInterp
#' @description Smooths a time series and interpolates for missing values
#' @details Uses Savitzgy-Golay polynomials to smooth the series, after
#' interpolating to fill missing values and then, if they still exist,
#' setting remaining missing values to zero.
#' @aliases SmoothInterp
#' @author William Cooper
#' @export SmoothInterp
#' @param .timeSeries A numeric vector containing a series to be smoothed.
#' @param .maxGap The largest gap across which to interpolate
#' @param .Length The length of the segment for Savitzky-Golay filtering.
#' (Must be odd; will be set odd if supplied as even.) Default: 61
#' @param .order The order of the polynomials to be used for filtering.
#' @return The smoothed and filtered series as a vector.
#' @examples 
#' \dontrun{Wsmoothed <- SmoothInterp (Data$WIC)}
SmoothInterp <- function (.timeSeries, .maxGap=1000, .Length=61, .order=3) {
  ## skip if there are fewer than 100 measurements
  if (length (.timeSeries[!is.na(.timeSeries)]) < 100) {return (.timeSeries)}
  d <- zoo::na.approx (as.vector(.timeSeries), maxgap=.maxGap, na.rm = FALSE)
  if (!(.Length %% 2)) {.Length <- .Length + 1}
  d[is.na(d)] <- 0
  return (signal::filter(signal::sgolay(.order, .Length), d))
}