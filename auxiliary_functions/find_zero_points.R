# Funzione che trova i punti di zero
find_zero_points <- function(x)
{
    # Find all the unique pixel ids
    unique_id_pixel <- x %>%
        select_("id_pixel") %>%
        distinct() %>%
        pull()
    
    # Number of unique pixels that will be examined
    print(paste("Number of pixels to be examined: ", length(unique_id_pixel)), sep = "")
    
    # List of zero points
    zero_points <- list()
    # This loop finds the zero points
    k <- 1
    for(i in unique_id_pixel)
    {
        print(paste("Finding zeros of pixel ", i))
        # Get time axis for pixel
        time_axis <- x %>% filter(id_pixel == i) %>% select(id_date) %>% pull()
        # Select derivative smoothed with running average (moving average)
        derivative <- x %>% filter(id_pixel == i) %>% select(D_mav) %>% pull()
        
        # Add a zero to have vector of same length.
        updn <- c(0, diff(sign(derivative)))
        # Find where there are changes. Index of the observation where there is a sign change.
        # if I have a,b,c and a > 0 and b < 0 and c < 0. I get the index of b
        ix <- which(updn != 0)
        # Number of zero points found
        print(paste("Zero points found:", length(ix)))
        # The bloom is in between these days
        print(paste(time_axis[ix], time_axis[ix - 1]))
        # Add found points to output list
        zero_points[[k]] <- c(time_axis[ix], time_axis[ix - 1])
        # Increment index
        k <- k + 1
    }
    
    names(zero_points) <- unique_id_pixel
    
    return(zero_points)
}
