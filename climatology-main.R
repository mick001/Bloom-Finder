#

require(qchlorophyll)
require(dplyr)

#-------------------------------------------------------------------------------
# Set paths

#nc_files_path <- "/mnt/hgfs/SHARE_VM/christian_paper/DATA_TEST"
nc_files_path <- "C:\\users\\michy\\desktop\\christian_paper\\DATA_TEST"

#-------------------------------------------------------------------------------
# Load data

# Carico file .nc ed estraggo CHL1_mean
nc_files_list <- load_all_as_list(path = nc_files_path, variables = c("CHL1_mean"))
# Unisco il tutto in un unico dataframe.
nc_dataframe <- assign_id_and_melt(nc_files_list)
# I dati caricati pronti per le analisi
#View(nc_dataframe)

#-------------------------------------------------------------------------------
# Calculate: 1. Climatology, 2. Required indeces.

# Calcola la media per pixel per data
climatology <- nc_dataframe %>%
    group_by(id_pixel, id_date) %>%
    summarise(avg_chl = mean(CHL1_mean, na.rm=T))

# Set NaN to zero for the purpose of this script
climatology <- climatology %>%
    mutate(avg_chl = case_when(
        is.na(avg_chl) ~ 0,
        TRUE ~ identity(avg_chl)
    ))

# Climatologia
#View(climatology)

# Calcola s, A, D, D_mav (moving average)
climatology <- climatology %>%
    ungroup() %>%
    # Already grouped by id_pixel
    group_by(id_pixel) %>%
    mutate(s = median(avg_chl, na.rm=T)*(1+0.05), # S: threshold median + 5%
           A = avg_chl - s, # anomalie
           C = cumsum(A), # somma cumulata
           D = (C - dplyr::lag(C))/8, # derivata temporale
           D_mav = TTR::runMean(D, 3)) %>% # Applico moving average alla derivata
    ungroup()

# Climatologia + indici
#View(climatology)
#climatology

# Check: should yield TRUE
#identical(diff(climatology$C), (climatology$C - dplyr::lag(climatology$C))[2:length(climatology$C)])

# Curiosità guardiamo le serie storiche del pixel 1.
# pixel1 <- climatology %>% ungroup() %>% filter(id_pixel == 1) %>% select(-id_pixel, -s)
# plot(pixel1$id_date, pixel1$avg_chl, type="l")
# plot(pixel1$id_date, pixel1$A, type="l")
# plot(pixel1$id_date, pixel1$C, type="l")
# plot(pixel1$id_date, pixel1$D, type="l")
# plot(pixel1$id_date, pixel1$D_mav, type="l")


# Funzione che trova i punti di zero
find_zero_points <- function(x)
{
    # Find unique pixels and get id
    unique_id_pixel <- x %>% select(id_pixel) %>% distinct() %>% pull()
    # Number of unique pixels that will be examined
    print(length(unique_id_pixel))
    
    # List of zero points
    zero_points <- list()
    
    # This loop finds the zero points
    for(i in unique_id_pixel)
    {
        print(paste("Finding zeros of ", unique_id_pixel[i]))
        # Time axis
        time_axis <- x %>% filter(id_pixel == i) %>% select(id_date) %>% pull()
        # Derivative smoothed with running average
        derivative <- x %>% filter(id_pixel == i) %>% select(D_mav) %>% pull()
        
        # Add a zero to have vector of same length.
        updn <- c(0, diff(sign(derivative)))
        # Find where there are changes. Index of the observation where there is a sign change.
        # if I have a,b,c and a > 0 and b < 0 and c < 0. I get the index of b
        ix <- which(updn != 0)
        # Number of zero points found
        print(paste("Zero points found:", length(ix)))
        # The bloom is in between these days
        print(paste(time_axis[ix], time_axis[ix-1]))
        #stop(msg = "Stopped by user")
        zero_points[[i]] <- c(time_axis[ix], time_axis[ix-1])
    }
    return(zero_points)
}

zero_pts <- find_zero_points(climatology)

find_number_of_blooms <- function(x)
{
    return( sapply(x, function(x1){ length(x1)/4 }) )
}

n_blooms <- find_number_of_blooms(zero_pts)
n_blooms
table(n_blooms)

# L'intersezione con lo zero si trova indicando due punti quello prima e dopo l'intersezione.
# quindi: punti trovati = n_bloom * 4

# N punti di zero = a n_bloom * 4/2

# n_bloom    1      1.5    2     2.5    3     3.5    4    4.5    5 
# frequency  1507   27    1167   13    403    5      68    1     9 










# This function generates a random cluster
generate_cluster <- function(size=4007)
{
    allowed_lengths <- 75:125
    cluster <- rep(1, sample(allowed_lengths, 1))
    i <- 2
    
    while(length(cluster) < size)
    {
        cluster <- c(cluster, rep(i, sample(allowed_lengths, 1)))
        i <- i + 1
    }
    
    return(cluster[1:size])
}
