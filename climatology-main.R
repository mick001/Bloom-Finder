#-------------------------------------------------------------------------------
# Clean workspace
rm(list=ls())

#-------------------------------------------------------------------------------
# Load required libraries

# Libraries needed
# TTR 0.23-3
require(mice)         # 2.46.0
require(qchlorophyll) # 2.1
require(dplyr)        # 0.7.4

#-------------------------------------------------------------------------------
# PARAMETERS

# Maximum number of allowed consecutive NAs. Mantiene solo i pixel con al più n NA consecutivi
MAX_CONSECUTIVE_NA <- 2
# Threshold percentage needed to calculate s index: median(chl) * (1 + THRESHOLD_PERCENTAGE). Defaults at 5%
THRESHOLD_PERCENTAGE <- 0.05
# Running average (or moving average) parameter
RUNNING_AVERAGE_WINDOW <- 3

#-------------------------------------------------------------------------------
# Set paths

# Path of .nc files
#nc_files_path <- "/mnt/hgfs/SHARE_VM/christian_paper/DATA_TEST"
nc_files_path <- "C:\\users\\michy\\desktop\\christian_paper\\DATA_TEST"
# Auxiliary functions path
aux_functions_path <- "C:\\users\\michy\\desktop\\christian_paper\\SCRIPT\\auxiliary_functions"

#-------------------------------------------------------------------------------
# Load data

# Carico file .nc ed estraggo CHL1_mean
nc_files_list <- load_all_as_list(path = nc_files_path, variables = c("CHL1_mean"))
# Unisco il tutto in un unico dataframe.
nc_dataframe <- assign_id_and_melt(nc_files_list)
# I dati caricati sono pronti per le analisi
#View(nc_dataframe)

rm(nc_files_list, nc_files_path)
#-------------------------------------------------------------------------------
# Calculate:
#           1. Climatology.
#           2. Required indeces.


# Load function to calculate consecutive NAs in climatology
source(file.path(aux_functions_path, "consecutive_na_count.R"))

# Calcola la media per pixel per data (climatologia)
climatology <- nc_dataframe %>%
    # For each pixel, then for each date
    group_by(id_pixel, id_date) %>%
    # Calculate: climatology, i.e. average value for date for pixel (avg_chl)
    #           number of observations used in each date (n_observations_used_per_date)
    summarise(avg_chl = mean(CHL1_mean, na.rm=T),
              n_observations_used_per_date = sum(!is.na(CHL1_mean))) %>%
    # Calculate, for each pixel: how many missing data in the climatology (NA_in_climatology_per_pixel) 
    #                           should the pixel be kept? (keep_pixel_NA_consecutive: TRUE keep, FALSE drop)
    mutate(NA_in_climatology_per_pixel = sum(is.na(avg_chl)),
           # If the number of consecutive NAs is > n then pixel should be removed from analysis
           keep_pixel_NA_consecutive = pixel_consecutive_NA(avg_chl, n = MAX_CONSECUTIVE_NA)) %>%
    # Remove grouping by pixel
    ungroup() %>%
    # Calculate initial number of pixels input in the analysis (pixel_in_analysis)
    mutate(pixels_in_analysis = n_distinct(id_pixel)) %>%
    # Keep only pixels with less than n consecutive NAs
    filter(keep_pixel_NA_consecutive == TRUE) %>%
    # Calculate: final number of pixels to be used for the analysis (pixels_out_analysis)
    #           percentage of pixels kept (pixels_kept_percentage)
    mutate(pixels_out_analysis = n_distinct(id_pixel),
           total_pixels_kept_percentage = pixels_out_analysis / pixels_in_analysis * 100) %>%
    # Remove unused columns
    select(-keep_pixel_NA_consecutive,
           -pixels_in_analysis,
           -pixels_out_analysis)

rm(pixel_consecutive_NA, MAX_CONSECUTIVE_NA)
#-------------------------------------------------------------------------------
# Add info about lon and lat

climatology <- nc_dataframe %>%
    select(lon, lat, id_pixel) %>%
    distinct() %>%
    left_join(climatology, ., by="id_pixel")

#-------------------------------------------------------------------------------
# Interpolate missing data with mice

# Obtain the mice object by interpolating with mice
imputed_df <- climatology %>%
    select(avg_chl, lon, lat) %>%
    mice(method = "pmm", m = 1) # m = 1 for fast testing: default is 5

# Show density plot of imputed vs actual data
densityplot(imputed_df, main = "Density of imputed data vs actual data")

# Complete imputed df # 1
imputed_df <- complete(imputed_df, 1)

# Add interpolated chl to climatology
climatology$avg_chl_interpolated <- imputed_df$avg_chl

rm(imputed_df)
#-------------------------------------------------------------------------------
# Ora in dataframe climatology:

# - id_pixel
# - id_date
# - avg_chl
# - lon
# - lat
# - n_observations_used_per_date: numero di osservazioni utilizzate per calcolare la
#                               climatologia in quella data per quel pixel
# - NA_in_climatology_per_pixel: numero di NA presenti in climatologia per quel pixel
# - total_pixels_kept_percentage: percentuale di pixel mantenuti sul totale di pixel analizzati
# - avg_chl_interpolated: climatologia interpolata con mice.

#-------------------------------------------------------------------------------
# Calculate useful indeces over the interpolated climatology

# Calcola s, A, D, D_mav (moving average)
climatology <- climatology %>%
    # Already ungrouped.
    group_by(id_pixel) %>%
    # For each pixel, calculate: s = median + a % of median
    #                           Anomalies = climatology - s
    #                           C = cumulative sum of anomalies
    #                           D = time derivative of cumulative sum of anomalies
    #                           D_mav = moving average (or running average of derivative)
    mutate(s = median(avg_chl_interpolated)*(1 + THRESHOLD_PERCENTAGE), # S: threshold median + 5%
           A = avg_chl_interpolated - s, # anomalie
           C = cumsum(A), # somma cumulata
           D = (C - dplyr::lag(C)) / 8, # derivata temporale
           D_mav = TTR::runMean(D, RUNNING_AVERAGE_WINDOW)) %>% # Applico moving average alla derivata
    # Remove grouping by pixel
    ungroup()

#Curiosità guardiamo le serie storiche del pixel 1.
pixel1 <- climatology %>% filter(id_pixel == 1730)
plot(pixel1$id_date, pixel1$avg_chl, type="l")
plot(pixel1$id_date, pixel1$avg_chl_interpolated, type="l")
plot(pixel1$id_date, pixel1$A, type="l")
plot(pixel1$id_date, pixel1$C, type="l")
plot(pixel1$id_date, pixel1$D, type="l")
plot(pixel1$id_date, pixel1$D_mav, type="l")

rm(pixel1, RUNNING_AVERAGE_WINDOW, THRESHOLD_PERCENTAGE)
#-------------------------------------------------------------------------------
# Find zero points

source(file.path(aux_functions_path, "find_zero_points.R"))
source(file.path(aux_functions_path, "find_number_of_blooms.R"))
source(file.path(aux_functions_path, "build_table_zero_points_1.R"))

# Find zero points of each climatology
zero_pts <- find_zero_points(climatology)

# Find number of blooms
n_blooms <- find_number_of_blooms(zero_pts)

# L'intersezione con lo zero si trova indicando due punti quello prima e dopo l'intersezione.
# -> N punti di zero = n_bloom * 4/2
# -> Punti trovati = n_bloom * 4
# -> -> Numero di bloom trovati = Punti trovati / 4

# Barplot of frequency of number of blooms found
table(n_blooms)
barplot(table(n_blooms), main = "Frequency of number of blooms found")

# The data is arranged in a dataframe
zero_points_df <- build_table(zero_pts, n_blooms)

rm(n_blooms, zero_pts)
#-------------------------------------------------------------------------------
# Flag pixels with more than 3 blooms and calculate bloom length

zero_points_df <- zero_points_df %>%
    mutate(flagged = n_blooms >= 3 ) %>%
    arrange(id_pixel, id_date_zero_crossing)

# Extract only valid pixels (less than 3 blooms)
valid_pixels_df <- zero_points_df %>%
    filter(flagged == TRUE)

# Interpolate new pixels and then find zero points again

unique_valid_pixels <- valid_pixels_df %>% select(id_pixel) %>% distinct() %>% pull() 

d <- climatology %>% filter(id_pixel == unique_valid_pixels[1]) %>% select(id_date, avg_chl_interpolated)
time_series <- data.frame(id_date_extended = 1:328, id_pixel = rep(unique_valid_pixels[1], 328)) %>% tbl_df()
time_series <- left_join(time_series, d, by = c("id_date_extended"="id_date"))
time_series$avg_chl_interpolated2 <- imputeTS::na.interpolation(time_series$avg_chl_interpolated, option="spline")
for(i in 2:length(unique_valid_pixels))
{
    d <- climatology %>% filter(id_pixel == unique_valid_pixels[i]) %>% select(id_date, avg_chl_interpolated)
    time_series_1 <- data.frame(id_date_extended = 1:328, id_pixel = rep(unique_valid_pixels[i], 328)) %>% tbl_df()
    time_series_1 <- left_join(time_series_1, d, by = c("id_date_extended"="id_date"))
    time_series_1$avg_chl_interpolated2 <- imputeTS::na.interpolation(time_series_1$avg_chl_interpolated, option="spline")
    
    time_series <- rbind(time_series, time_series_1)
}

View(time_series)

n <- 5
plot(1:328, time_series$avg_chl_interpolated2[((n-1) * 328+1):(n*328)], type="l")
lines(1:41*8, time_series$avg_chl_interpolated[((n-1) * 328+1):(n*328)][!is.na(time_series$avg_chl_interpolated[((n-1) * 328+1):(n*328)])], col="red", lty=10)
