################################################################################
#' This function plots the indeces calculated for the climatology
#' 
#' It plots:
#' - RAW chl climatology
#' - chl climatology interpolated
#' - Anomalies
#' - Cumulative sum of anomalies
#' - Derivative of cumulative sum of anomalies
#' - Moving average of derivative of cumulative sum of anomalies
#' 
#' @param id_plixel_to_plot id of pixel to be plotted
#' @param x climatology dataframe
#' @return void
#' @example plot_calculated_indeces(1729)
#' 
plot_calculated_indexes <- function(id_pixel_to_plot, x)
{
    # Set plot device, 2 rows and 3 columns
    par(mfrow = c(2, 2))
    # Get pixel data
    pixel <- x %>% filter(id_pixel == id_pixel_to_plot)
    
    # Plot
    if( nrow(pixel) > 0)
    {
        print(paste("Plotting data for id_pixel: ", id_pixel_to_plot, sep = ""))
        
        # Raw chl in climatology
        plot(pixel$id_date, pixel$avg_chl, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Avg_chl", main = "Climatology (RAW data)")
        lines(pixel$id_date, pixel$avg_chl, lwd = 2, col = "blue")
        # Interpolated chl (linear interpolation) in climatology
        plot(pixel$id_date, pixel$avg_chl_interpolated, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Avg_chl_interpolated", main = "Climatology with no NAs")
        lines(pixel$id_date, pixel$avg_chl_interpolated, lwd = 2, col = "blue")
        # Anomalies
        plot(pixel$id_date, pixel$A, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Anomalies", main = "Anomalies")
        lines(pixel$id_date, pixel$A, lwd = 2, col = "blue")
        # Filtered (smoothed) anomalies
        plot(pixel$id_date, pixel$D_mav, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Smoothed anomalies", main = "Smoothed A")
        lines(pixel$id_date, pixel$D_mav, lwd = 2, col = "blue")
        abline(0, 0, lwd=2)
        
        # Restore plot device
        par(mfrow = c(1, 1))
    }else
    {
        print(paste("Pixel ", id_pixel_to_plot, " is not in climatology...", sep = ""))
    }
}
