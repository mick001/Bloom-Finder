################################################################################
#' This function allows to use different smoothing functions on a numeric vector
#' 
#' Implemented smoothing functions:
#' - Simple moving average
#' - Kernel smoother
#' - Spline smoother
#' - Loess smoother
#' 
#' @param x a numeric vector
#' @param filter_name filter function to be used
#' @param kernel Kernel type
#' @param bandwidth Gaussian kernel parameter
#' @param df df spline parameter
#' @param span loess span parameter
#' @return a numeric vector
#' @example filter_function(x, name="mav", n=5)
#' 
filter_function <- function(x,
                            filter_name = c("mav", "kernelSm", "spline", "loess"),
                            n = 50, # MAV param
                            kernel = "normal", # Kernel type
                            bandwidth = 50, # Gauss kernel param
                            df = 35, # spline param
                            span = 0.3, # loess param
                            ...)
{
    # Match argument
    filter_name <- match.arg(filter_name)
    # Choose appropriate filter
    switch(filter_name,
           mav = TTR::runMean(x = x,
                              n = n,
                              ...),
           kernelSm = ksmooth(x = 1:length(x),
                              y = x,
                              kernel = kernel,
                              bandwidth = bandwidth,
                              ...)$y,
           spline = smooth.spline(x = 1:length(x),
                                  y = x,
                                  df = df,
                                  ...)$y,
           loess = loess(y ~ x1,
                         data = data.frame(x1 = 1:length(x), y = x),
                         span = span,
                         ...)$fitted )
}
