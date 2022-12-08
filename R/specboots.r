#' Bootstrapping Spectral Density
#'
#' This function bootstraps kernel spectral density estimates by resampling the periodogram of the original data object.
#' The result is a plot of the bootstrap confidence interval and a dataframe that contains the frequencies and
#' the lower and upper bounds of the estimated spectral density.
#'
#' @import astsa
#' @import stats
#' @param x a time series
#' @param nboot the number of bootstrap series to compute
#' @param spans specify smoothing
#' @param kernel specify kernel
#' @param detrend if TRUE, data is detrended first (TRUE as default)
#' @param demean if TRUE, data is demeaned first (FALSE as default)
#' @param base.stat a function applied to return the baseline statistic of interest (median as default)
#' @param alpha a confidence interval (0.95 as default)
#' @returns
#' A dataframe of frequencies, lower and upper spectral density bounds.
#' \itemize{
#' \item{freq: frequencies at which spectral density are computed}
#' \item{lwr: the spectral density lower bound at all frequencies}
#' \item{upr: the spectral density upper bound at all frequencies}
#' }
#' @examples
#' data("AirPassengers")
#' tsplot(AirPassengers, col=4, lwd=2)
#' kd = kernel("daniell", c(5,5))
#' spec.boots(AirPassengers, 1000, kd, periodogram_plot = FALSE)
#' @references
#' \itemize{
#' \item{Franke, J., and W. Hardle. “On Bootstrapping Kernel Spectral Estimates.” The Annals of Statistics, vol. 20, no. 1, 1992, pp. 121–45. JSTOR, http://www.jstor.org/stable/2242153.}
#' \item{Zoubir, Abdelhak M.. “Bootstrapping spectra: Methods, comparisons and application to knock data”. Signal Process. 90, 5. May, 2010, 1424–1435. https://doi.org/10.1016/j.sigpro.2009.11.030.}
#' }
#' @export
#'
spec.boots <- function(
    x,
    nboot,
    spans = NULL,
    kernel = NULL,
    detrend=TRUE,
    demean = FALSE,
    base.stat = median,
    alpha = 0.95)
{
  # check if x is a time series
  if (!is.ts(x))
  {
    stop("x must be a time series")
  }
  # check if the spans and kernel are specified and valid
  if (is.null(spans) && is.null(kernel))
  {
    stop("must specify spans or a valid kernel")
  }

  if (!is.null(spans))
    kernel <- {
      if (is.tskernel(spans))
        spans
      else kernel("modified.daniell", spans%/%2)
    }

  if (length(kernel[[1]]) > length(x))
  {
    stop("'x' is shorter than the kernel spans")
  }

  # check if alpha is valid
  if (alpha >= 1 || alpha <= 0)
  {
    stop("alpha must be between 0 and 1")
  }

  # extract the frequencies
  n = length(x)
  xfreq = frequency(x)
  nspec = floor(n/2)
  freq = seq(from = xfreq/n, by = xfreq/n, length = nspec)

  # step 1: centering
  if (detrend)
  {
    x = detrend(x)
  }
  else if (demean)
  {
    x = x - mean(x)
  }

  # step 2: initial estimate
  I = Mod(fft(x)/sqrt(n))^2
  I = I[1:floor(n/2)]
  f.hat = kernapply(I, kernel, circular=TRUE)

  # step 3: compute and rescale residuals
  eps = I / f.hat
  scaled_eps = eps / mean(eps)

  # step 4-6: resample residuals and bootstrap spectral denisty
  f.star.dist = c(c())
  for (i in 1:nboot)
  {
    scaled_eps.star = sample(scaled_eps, replace=TRUE)  # resample the scaled residuals
    I.star = f.hat * scaled_eps.star # resample periodogram
    f.star = kernapply(I.star, kernel, circular=TRUE) # bootstrap kernel spectral density
    f.star.dist[i] = list(f.star)
  }

  # step 7: confidence interval estimate
  lwr = c()
  upr = c()
  base = c()
  for (i in 1:nspec)
  {
    f.star.wi = c()
    for (j in 1:nboot)
    {
      f.star.wi[j] = f.star.dist[[j]][i]
    }
    # the lower percentile spectral density at all frequencies
    lwr[i] = quantile(f.star.wi, probs=c((1-alpha)/2))
    # the higher percentile spectral density at all frequencies
    upr[i] = quantile(f.star.wi, probs=c((1+alpha)/2))
    # median as the baseline
    base[i] = base.stat(f.star.wi)
  }

  tsplot(freq, upr, col=4, lty=2, lwd=2,
         ylim=c(min(lwr), max(upr)),
         main="Spectral Density Bootstrap",
         ylab="Spectral Density", xlab="Frequency")
  lines(freq, lwr, col=4, lty=2, lwd=2)
  abline(h = base.stat(base), col=2)

  interval = data.frame(freq=freq, lwr=lwr, upr=upr)
  return(interval)
}
