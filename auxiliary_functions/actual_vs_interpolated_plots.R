################################################################################
#' This function plots the low and high resolution data for the climatology
#' 
#' It plots:
#' - Chl climatology interpolated
#' - Moving average of derivative of cumulative sum of anomalies
#' 
#' @param id_plixel_to_plot vector to be checked for consecutive missing data
#' @param x high resolution climatology
#' @return void
#' @example actual_vs_interpolated_plots(1729)
#' 
actual_vs_interpolated_plots <- function(id_pixel_to_plot, x)
{
   
    
    unique_valid_pixels <- x %>%
        select_("id_pixel") %>%
        distinct() %>%
        pull()
    
    n <- which(unique_valid_pixels == id_pixel_to_plot)
    
    if(length(n) != 0)
    {
        # Set plotting device
        par(mfrow = c(1, 2))
        
        #---------------------------------------------------------------------------
        # Climatology
        
        print(paste("Checking pixel: ", as.character(unique_valid_pixels[n]), sep=""))
        # High resolution avg_chl, points and lines
        plot(0:360, x$avg_chl_interpolated_high_res[((n - 1) * 361 + 1):(n * 361)],
             col = "blue",
             pch = 20,
             ylab = "avg_chl", xlab = "id_date",
             main = paste("Actual vs interpolated avg_chl pixel: ", id_pixel_to_plot, sep = "") )
        lines(0:360, x$avg_chl_interpolated_high_res[((n - 1) * 361 + 1):(n * 361)],
              type = "l",
              col = "blue",
              lwd = 2)
        # Low resolution avg_chl
        points(0:45*8, x$avg_chl_interpolated[((n - 1) * 361+1):(n * 361)][!is.na(x$avg_chl_interpolated[((n-1) * 361+1):(n * 361)])],
               col = "red",
               lty = 10,
               lwd = 2)
        # Legend
        legend("topright", legend = c("Stine intp", "Original data"),
               col=c("blue", "red"), pch = c(20, 20), cex = 0.8)
        
        #---------------------------------------------------------------------------
        # Moving average
        
        plot(0:360, (x$D_mav_high_res_from_stine[((n - 1) * 361 + 1):(n * 361)])[1:361],
             col = "blue",
             pch = 20,
             ylab = "D_mav", xlab = "id_date",
             main = paste("Interpolated D_mav pixel: ", id_pixel_to_plot, sep = "") )
        lines(0:360, (x$D_mav_high_res_from_stine[((n - 1) * 361 + 1):(n * 361)])[1:361],
              col = "blue",
              lwd = 2)
        
        # NON SI PUò PIù FARE IL PLOT PER VIA DELLO SHIFT...
        # Actual moving average
        # points(c(7, (2:43)*8 - 1), x$D_mav[((n-1) * 361+1):(n*361)][!is.na(x$D_mav[((n-1) * 361+1):(n*361)])],
        #        col = "red",
        #        pch = 20)
        
        # Zero line
        abline(0, 0, lwd=2)
        # Legend
        legend("topright", legend = c("Stine intp"),
               col = c("blue"), pch = c(20), cex = 0.8)
        
        # Restore plotting device
        par(mfrow = c(1, 1))
        
    }else
    {
        print(paste("Pixel does not exist: ", id_pixel_to_plot))
    }
}