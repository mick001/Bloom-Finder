build_table <- function(zero_pts, n_blooms)
{
    list_names <- names(zero_pts)
    data_points <- zero_pts[[1]]
    df_out <- data.frame(id_pixel = rep(list_names[1], length(data_points)),
                         id_date_zero_crossing = data_points,
                         n_blooms = rep(n_blooms[1], length(data_points)),
                         n_zero_points = rep(n_blooms[1], length(data_points))*2,
                         stringsAsFactors = F)
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
    
    return(df_out)
}