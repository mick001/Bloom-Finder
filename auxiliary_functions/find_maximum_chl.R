################################################################################
#' This function finds maximum of chl and corresponding date for TABELLA_DUE during bloom
#' 
#' @return A dataframe with maximum chl and corresponding date
#' 
find_maximum_chl <- function()
{
    max_chl <- list()
    date_max_chl <- list()
    
    for(i in 1:nrow(TABELLA_DUE))
    {
        current_id_pixel <- TABELLA_DUE[i, c("id_pixel")] %>% pull()
        bloom_start <- TABELLA_DUE[i, c("bloom_start_date")] %>% pull()
        bloom_end <- TABELLA_DUE[i, c("bloom_end_date")] %>% pull()
        
        current_df <- climatology_high_res %>%
            filter(id_pixel == current_id_pixel) %>%
            filter( (id_date_extended >= bloom_start) & (id_date_extended <= bloom_end) ) %>%
            select(avg_chl_interpolated_high_res, id_date_extended) %>%
            filter(avg_chl_interpolated_high_res == max(avg_chl_interpolated_high_res))
        
        max_chl[[i]] <- current_df[1, ] %>%
            pull(avg_chl_interpolated_high_res)
        
        date_max_chl[[i]] <- current_df[1, ] %>%
            pull(id_date_extended)
        
    }
    
    out_df <- data.frame(max_chl = as.numeric(max_chl),
                         id_date_max_chl = as.numeric(date_max_chl)) %>%
        as_tibble()
    
    return(out_df)
}
