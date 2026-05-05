#' @title SmoothInterp
#' @description Smooths a time series and interpolates for missing values
#' @details Uses Savitzgy-Golay polynomials to smooth the series, after
#' interpolating to fill missing values and then, if they still exist,
#' setting remaining missing values to zero.
#' @aliases smoothInterp
#' @author William Cooper
#' @export SmoothInterp
#' @importFrom signal filter sgolay
#' @param .timeSeries A numeric vector containing a series to be smoothed.
#' @param .maxGap The largest gap across which to interpolate
#' @param .Length The length of the segment for Savitzky-Golay filtering.
#' (Must be odd; will be set odd if supplied as even.) Default: 61
#' If .Length <= 1 there will be no smoothing, only interpolation.
#' @param .order The order of the polynomials to be used for filtering.
#' @param .minMeas The minimum series length (not including NAs) to
#' perform interpolation. Default: 100. This value previously was coded
#' into the routine, but caused unreported problems for short series, so
#' this parameter was added with a default equal to the previous setting.
#' @return The smoothed and filtered series as a vector.
#' @examples 
#' DPsmoothed <- SmoothInterp (RAFdata$DPXC)

SmoothInterp <- function (.timeSeries, .maxGap=1000, .Length=61, 
                          .order=3, .minMeas = 100) {
  ## skip if there are fewer than .minNeas measurements
  if (length (.timeSeries[!is.na(.timeSeries)]) < .minMeas) {return (.timeSeries)}
  ## tried na.spline, but it often introduces too much variance, so returned to na.approx:
  d <- zoo::na.approx (as.vector(.timeSeries), maxgap=.maxGap, na.rm = FALSE, rule=2)
  if (!(.Length %% 2)) {.Length <- .Length + 1}
  # Special handling of first or last item:
  # if (is.na(d[1])) {d[1] <- d[2]}
  # if (is.na(d[length(d)])) {d[length(d)] <- d[length(d)-1]}
  d[is.na(d)] <- 0
  if (.Length > 2) {
    return (as.vector (signal::filter(signal::sgolay(.order, .Length), d)))
  } else {
    return (d)
  }
}
