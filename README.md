# spec.boots: Bootstrapping Spectral Density 

# Introduction 

In time series spectral analysis, we often want to identify the dominant frequencies (or periods) of the observed series. A periodogram gives us a rough sample estimate of the population spectral density but never gets close to the spectrum as the data size gets large. Nevertheless, we can still construct a nonparametric confidence interval for the spectral density based on the periodogram with a few data assumptions. In reality, data does not always follow these assumptions. That is where bootstrapping spectral density becomes handy. This package implements a bootstrap approach to estimate the spectral density by resampling the periodogram of the original data. 

# Installation 

``` r
# install.packages("devtools")
devtools::install_github("uyenle-gh/spec.boots")
```

# Main Feature 

The `spec.boots` function in the package allows us to create a confidence interval for the spectral density of a time series object. 

``` r 
library(spec.boots)
data("AirPassengers")
tsplot(AirPassengers, col=4, lwd=2)
spec.boots(AirPassengers, 1000, c(5,5))
```

# References

Franke, J., and W. Hardle. “On Bootstrapping Kernel Spectral Estimates.” The Annals of Statistics, vol. 20, no. 1, 1992, pp. 121–45. JSTOR, http://www.jstor.org/stable/2242153.

Zoubir, Abdelhak M.. “Bootstrapping spectra: Methods, comparisons and application to knock data”. Signal Process. 90, 5. May, 2010, 1424–1435. https://doi.org/10.1016/j.sigpro.2009.11.030.