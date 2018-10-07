################################################################################
#' This function arranges the information on zero points in a dataframe
#' 
#' @param zero_pts List with zero points found
#' @param n_blooms List with number of blooms found
#' @return A dataframe with collected data
#' 
#' 
build_table_zero_points <- function(zero_pts, n_blooms)
{
    # Nomi dei pixel
    list_names <- names(zero_pts)
    
    # Punti di zero del primo pixel
    data_points <- zero_pts[[1]]
    
    # Dataframe con output
    df_out <- data.frame(id_pixel = rep(list_names[1], length(data_points)),
                         id_date_zero_crossing = data_points,
                         n_blooms = rep(n_blooms[1], length(data_points)),
                         n_zero_points = rep(n_blooms[1], length(data_points))*2,
                         stringsAsFactors = F)
    # Popola dataframe
    for(i in 2:length(zero_pts))
    {
        data_points <- zero_pts[[i]]
        df <- data.frame(id_pixel = rep(list_names[i], length(data_points)),
                         id_date_zero_crossing = data_points,
                         n_blooms = rep(n_blooms[i], length(data_points)),
                         n_zero_points = rep(n_blooms[i], length(data_points))*2,
                         stringsAsFactors = F)
        df_out <- rbind(df_out, df)
    }
    
    # Convert to tibble
    df_out <- df_out %>%
        mutate(id_pixel = as.integer(id_pixel)) %>%
        as_tibble()
    
    return(df_out)
}
