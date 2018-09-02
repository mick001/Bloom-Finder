################################################################################
#' This function plots the low and high resolution data for the climatology
#' 
#' It plots:
#' - Chl climatology interpolated
#' - Moving average of derivative of cumulative sum of anomalies
#' 
#' @param id_plixel_to_plot vector to be checked for consecutive missing data
#' @return void
#' @example compare_data_interpolation(1729)
#' 
compare_data_interpolation <- function(id_pixel_to_plot)
{
    # Set plotting device
    par(mfrow = c(1, 2))
    
    n <- which(unique_valid_pixels == id_pixel_to_plot)
    
    if(length(n) != 0)
    {
        #---------------------------------------------------------------------------
        # Climatology
        
        print(paste("Checking pixel: ", as.character(unique_valid_pixels[n]), sep=""))
        # High resolution avg_chl, points and lines
        plot(1:328, climatology_high_res$avg_chl_interpolated_high_res[((n - 1) * 328 + 1):(n * 328)],
             col = "blue",
             pch = 20,
             ylab = "avg_chl", xlab = "id_date",
             main = paste("Actual vs interpolated avg_chl pixel: ", as.character(unique(climatology_high_res$id_pixel[((n - 1) * 328 + 1):(n * 328)])), sep = "") )
        lines(1:328, climatology_high_res$avg_chl_interpolated_high_res[((n - 1) * 328 + 1):(n * 328)],
              type = "l",
              col = "blue",
              lwd = 2)
        # Low resolution avg_chl
        points(1:41*8, climatology_high_res$avg_chl_interpolated[((n-1) * 328+1):(n*328)][!is.na(climatology_high_res$avg_chl_interpolated[((n-1) * 328+1):(n*328)])],
               col = "red",
               lty = 10,
               lwd = 2)
        # Legend
        legend("topright", legend = c("Stine intp", "Original data"),
               col=c("blue", "red"), pch = c(20, 20), cex = 0.8)
        
        #---------------------------------------------------------------------------
        # Moving average
        
        plot(32:328, (climatology_high_res$D_mav_high_res_from_stine[((n - 1) * 328 + 1):(n * 328)])[32:328],
             col = "blue",
             pch = 20,
             ylab = "D_mav", xlab = "id_date",
             main = paste("Actual vs interpolated D_mav pixel: ", as.character(unique(climatology_high_res$id_pixel[((n - 1) * 328 + 1):(n * 328)])), sep = "") )
        lines(32:328, (climatology_high_res$D_mav_high_res_from_stine[((n - 1) * 328 + 1):(n * 328)])[32:328],
              col = "blue",
              lwd = 2)
        # Actual moving average
        points(4:41*8, climatology_high_res$D_mav[((n-1) * 328+1):(n*328)][!is.na(climatology_high_res$D_mav[((n-1) * 328+1):(n*328)])],
               col = "red",
               pch = 20)
        # Zero line
        abline(0, 0, lwd=2)
        # Legend
        legend("topright", legend = c("Stine intp", "Original data"),
               col = c("blue", "red"), pch = c(20, 20), cex = 0.8)
        
        # Restore plotting device
        par(mfrow = c(1, 1))
    }else
    {
        warning(paste("Pixel does not exist: ", id_pixel_to_plot))
    }
}