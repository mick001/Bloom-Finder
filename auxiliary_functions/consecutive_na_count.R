################################################################################
#' This function calculates the number of consecutive NAs in a vector and
#' returns a boolean value
#' 
#' If number of consecutive NAs is > n, then F is returned
#' If number of consecutive NAs is <= n, then T is returned
#'
#' @param x vector to be checked for consecutive missing data
#' @param n number of maximum consecutive NAs allowed
#' @return boolean, as stated above
#' @example pixel_consecutive_NA(1:10), pixel_consecutive_NA(c(1:10, NA, NA, NA))
#' 
pixel_consecutive_NA <- function(x, n = 2)
{
    rle_out <- rle(is.na(x))
    rle_out_choose <- rle_out$values & (rle_out$lengths > n)
    
    if(any(rle_out_choose))
    {
        # Remove pixel!
        return(FALSE)
    }else
    {
        # DO NOT remove pixel!
        return(TRUE)
    }
}
