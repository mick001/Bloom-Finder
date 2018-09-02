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
#' @return void
#' @example plot_calculated_indeces(1729)
#' 
plot_calculated_indeces <- function(id_pixel_to_plot)
{
    print(paste("Plotting data for id_pixel: ", id_pixel_to_plot, sep = ""))
    
    # Set plot device, 2 rows and 3 columns
    par(mfrow = c(2,3))
    # Get pixel data
    pixel <- climatology %>% filter(id_pixel == id_pixel_to_plot)
    
    # Plot
    
    plot(pixel$id_date, pixel$avg_chl, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Average chl", main = "Climatology on RAW data")
    lines(pixel$id_date, pixel$avg_chl, lwd = 2, col = "blue")
    
    plot(pixel$id_date, pixel$avg_chl_interpolated, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Average chl intp", main = "Climatology interpolated")
    lines(pixel$id_date, pixel$avg_chl_interpolated, lwd = 2, col = "blue")
    
    plot(pixel$id_date, pixel$A, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Anomalies", main = "Anomalies")
    lines(pixel$id_date, pixel$A, lwd = 2, col = "blue")
    
    plot(pixel$id_date, pixel$C, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "C. sum", main = "Cumulative sum of anomalies")
    lines(pixel$id_date, pixel$C, lwd = 2, col = "blue")
    
    plot(pixel$id_date, pixel$D, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Derivative", main = "Derivative of cumulative sum")
    lines(pixel$id_date, pixel$D, lwd = 2, col = "blue")
    
    plot(pixel$id_date, pixel$D_mav, type="p", pch = 18, col = "red", xlab = "id_date", ylab = "Moving average", main = "Derivative with moving average")
    lines(pixel$id_date, pixel$D_mav, lwd = 2, col = "blue")
    
    # Restore plot device
    par(mfrow = c(1, 1))
}
