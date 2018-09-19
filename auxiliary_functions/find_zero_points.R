################################################################################
#' Given a dataframe with n timeseries (i.e. n pixels, that is a time series
#' for every pixel), this function finds the zero points of the time series
#' and returns a list.
#' 
#' The number of elements of the list is equal to the number of time series analyzed
#' (i.e. the number of pixels)
#' 
#' Every element of the list is the set of points that identify zero points in 
#' each time series.
#' 
#' @param x a dataframe containing a set of time series identified by a variable id_pixel
#' @param id_date_name name of id_date
#' @return a list as explained above
#' 
find_zero_points <- function(x, id_date_name="id_date", D_mav_name = "D_mav")
{
    # Find all the unique pixel ids
    unique_id_pixel <- x %>%
        select_("id_pixel") %>%
        distinct() %>%
        pull()
    
    # Number of unique pixels that will be examined
    print(paste("Number of pixels to be analysed: ", length(unique_id_pixel)), sep = "")
    
    # List of zero points
    zero_points <- list()
    # This loop finds the zero points
    k <- 1
    for(i in unique_id_pixel)
    {
        #print(paste("Finding zeros of pixel ", i))
        # Get time axis for pixel
        time_axis <- x %>% filter(id_pixel == i) %>% select_(id_date_name) %>% pull()
        # Select derivative smoothed with running average (moving average)
        derivative <- x %>% filter(id_pixel == i) %>% select_(D_mav_name) %>% pull()
        
        # Add a zero to have vector of same length.
        updn <- c(0, diff(sign(derivative)))
        # Find where there are changes. Index of the observation where there is a sign change.
        # if I have a,b,c and a > 0 and b < 0 and c < 0. I get the index of b
        ix <- which(updn != 0)
        # Number of zero points found
        #print(paste("Zero points found:", length(ix)))
        # The bloom is in between these days
        #print(paste(time_axis[ix], time_axis[ix - 1]))
        # Add found points to output list
        zero_points[[k]] <- c(time_axis[ix], time_axis[ix - 1])
        # Increment index
        k <- k + 1
    }
    
    names(zero_points) <- unique_id_pixel
    
    return(zero_points)
}
