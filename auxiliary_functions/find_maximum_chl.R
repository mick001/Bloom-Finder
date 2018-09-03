################################################################################
#' This function finds maximum of chl for TABELLA_DUE during bloom
#' 
#' @return A list of max chl for every pixel and every bloom
#' 
find_max_chl <- function()
{
    max_chl <- list()
    for(i in 1:nrow(TABELLA_DUE))
    {
        current_id_pixel <- TABELLA_DUE[i, c("id_pixel")] %>% pull()
        bloom_start <- TABELLA_DUE[i, c("bloom_start_date")] %>% pull()
        bloom_end <- TABELLA_DUE[i, c("bloom_end_date")] %>% pull()
        
        max_chl[[i]] <- climatology_high_res %>%
            filter(id_pixel == current_id_pixel) %>%
            filter( (id_date_extended >= bloom_start) & (id_date_extended <= bloom_end) ) %>%
            select(avg_chl_interpolated_high_res) %>%
            pull() %>%
            max()
    }
    return(max_chl)
}
