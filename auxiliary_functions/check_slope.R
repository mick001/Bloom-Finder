################################################################################
#' This function checks that the slope of the first zero point is positive.
#' 
#' The check is performed on every pixel in climatology
#' 
#' @return vector of pixels that can be further processed
#' 
check_slope <- function()
{
    # Get list of unique pixels
    pixels <- climatology %>%
        select(id_pixel) %>%
        distinct() %>%
        pull()

    # Pixel list to be returned
    out_pixels <- list()
    k <- 1
    # For each pixel, check that the slope in the first zero point is positive,
    # if yes pixel can be further analysed, otherwise no.
    for(i in 1:length(pixels))
    {
        # Pull out the moving average for the current pixel
        df <- climatology %>%
            filter(id_pixel == pixels[i]) %>%
            select(D_mav) %>%
            pull()
        
        # Find zero points
        updn <- c(0, diff(sign(df)))
        # Find where there are changes. Index of the observation where there is a sign change.
        # if I have a,b,c and a > 0 and b < 0 and c < 0. I get the index of b
        ix <- which(updn != 0)
        # The zero is in between these days
        #print(paste(time_axis[ix], time_axis[ix - 1]))
        
        # Find value of moving average before and after found zero point
        deriv_pt1 <- df[ix[1] - 1]
        deriv_pt2 <- df[ix[1]]
        
        # If point 2 is higher than point 1 then slope wrt time is positive!
        if(deriv_pt2 > deriv_pt1)
        {
            out_pixels[[k]] <- pixels[i]
            k <- k + 1
        }
    }

    # Pixels to be further analysed
    return(as.integer(out_pixels))
}
