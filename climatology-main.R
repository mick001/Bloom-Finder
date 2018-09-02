#-------------------------------------------------------------------------------
# Clean workspace
rm(list = ls())

# Set seed for reproducibility purposes
set.seed(332423)

#-------------------------------------------------------------------------------
# Load required libraries

# TTR                 # 0.23-3
# lazyeval            # 0.2.1
# log4r
require(mice)         # 2.46.0
require(qchlorophyll) # 2.1
require(dplyr)        # 0.7.4
require(scales)       # 0.5.0

#-------------------------------------------------------------------------------
# PARAMETERS

# Maximum number of allowed consecutive NAs in the climatology. Mantiene solo i pixel con al più n NA consecutivi
MAX_CONSECUTIVE_NA <- 2
# Threshold percentage needed to calculate s index = median(chl) * (1 + THRESHOLD_PERCENTAGE). Defaults at 5%
THRESHOLD_PERCENTAGE <- 0.05
# Running average (or moving average) parameter. Number of observations to be used in the moving average.
RUNNING_AVERAGE_WINDOW <- 3
# Percentile interval into which climatology data needs to be squished.
PERCENTILE_SQUISHING_INTERVAL <- c(0.05, 0.95)
# Type of average operation to be used. Set to "geom" for geometric mean.
MEAN_FUNCTION <- "mean"

#-------------------------------------------------------------------------------
# Set paths

# Path of .nc files
#NC_FILES_PATH <- "/mnt/hgfs/SHARE_VM/christian_paper/DATA_TEST"
NC_FILES_PATH <- "C:\\users\\michy\\desktop\\christian_paper\\DATA_TEST"
# Auxiliary functions path
AUX_FUNCTIONS_PATH <- "C:\\users\\michy\\desktop\\christian_paper\\SCRIPT\\auxiliary_functions"
# Output directory
OUTPUT_PATH <- "C:\\users\\michy\\desktop"

#-------------------------------------------------------------------------------
# Setup of logger

logger <- log4r::create.logger(logfile = file.path(OUTPUT_PATH, "log_climatology_main.log"),
                               level = "INFO")

#-------------------------------------------------------------------------------
# Load data

print("Starting to load data...")

# Carico file .nc ed estraggo CHL1_mean
nc_files_list <- load_all_as_list(path = NC_FILES_PATH, variables = c("CHL1_mean"))
# Unisco il tutto in un unico dataframe.
nc_dataframe <- assign_id_and_melt(nc_files_list)

rm(nc_files_list, NC_FILES_PATH)
#-------------------------------------------------------------------------------
# Calculate on RAW chl data the climatology

print("Calculating climatology...")

# Pixels used as input in the analysis
pixels_IN <- length(unique(nc_dataframe$id_pixel))

# Load function to calculate consecutive NAs in climatology
source(file.path(AUX_FUNCTIONS_PATH, "consecutive_na_count.R"))

# Set mean function to be used for calculating climatology
if(MEAN_FUNCTION == "mean")
{
    MEAN_FUNCTION <- function(x){ mean(x, na.rm=T) }
    print("Using arithmetic mean...")
    log4r::info(logger, "Using arithmetic mean for climatology calculation.")
}else if(MEAN_FUNCTION == "geom")
{
    MEAN_FUNCTION <- function(x){ exp(mean(log(x), na.rm = T)) }
    print("Using geometric mean...")
    log4r::info(logger, "Using geometric mean for climatology calculation.")
}else
{
    MEAN_FUNCTION <- function(x){ mean(x, na.rm=T) }
    warning("Mean function specified is not correct... using arithmetic mean")
    log4r::info(logger, "Using arithmetic mean for climatology calculation.")
}

# Calcola la media specificata per pixel per data (climatologia)
climatology <- nc_dataframe %>%
    # For each pixel, then for each date
    group_by(id_pixel, id_date) %>%
    # Calculate: climatology, i.e. average value for date for pixel (avg_chl)
    #           number of observations used in each date (n_observations_used_per_date)
    
    #summarise(avg_chl = mean(CHL1_mean, na.rm=T),
     #         n_observations_used_per_date = sum(!is.na(CHL1_mean))) %>%
    summarise_(.dots = setNames(list(lazyeval::interp( ~ MEAN_FUNCTION(CHL1_mean)),
                                     lazyeval::interp( ~ sum(!is.na(CHL1_mean))) ),
                                c("avg_chl",
                                  "n_observations_used_per_date"))) %>%
    
    # Calculate, for each pixel: how many missing data in the climatology (NA_in_climatology_per_pixel) 
    #                           should the pixel be kept? (keep_pixel_NA_consecutive: TRUE keep, FALSE drop)
    mutate(NA_in_climatology_per_pixel = sum(is.na(avg_chl)),
           # If the number of consecutive NAs is > n then pixel should be removed from analysis
           keep_pixel_NA_consecutive = pixel_consecutive_NA(avg_chl, n = MAX_CONSECUTIVE_NA)) %>%
    # Remove grouping by pixel
    ungroup() %>%
    # Keep only pixels with less than n consecutive NAs
    filter(keep_pixel_NA_consecutive == TRUE) %>%
    # Remove unused columns
    select(-keep_pixel_NA_consecutive)

# Pixels kept for the analysis
pixels_kept <- length(unique(climatology$id_pixel))

# Percentage of pixels kept for the analysis on the total of pixels used as input
print(paste("Percentage of pixels kept: ", pixels_kept/pixels_IN * 100, sep = ""))

# Log the information
log4r::info(logger, paste("Pixels used as input: ", pixels_IN, sep = ""))
log4r::info(logger, paste("Pixels used for analysis: ", pixels_kept, sep = ""))
log4r::info(logger, paste("Percentage of pixels used for analysis over input: ", pixels_kept/pixels_IN * 100, sep = ""))

rm(pixel_consecutive_NA, MAX_CONSECUTIVE_NA, pixels_IN, pixels_kept)
#-------------------------------------------------------------------------------
# Add info on longitude and latitude (lon and lat)

print("Adding info on longitude and latitude...")

climatology <- nc_dataframe %>%
    select(lon, lat, id_pixel) %>%
    distinct() %>%
    left_join(climatology, ., by = "id_pixel")

#-------------------------------------------------------------------------------
# Interpolate missing data with MICE

print("Interpolating climatology with MICE...")

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
# Squish of climatology data in the selected percentile interval

print("Squishing climatology in the selected percentile interval...")

climatology <- climatology %>%
    mutate(avg_chl_interpolated = squish(avg_chl_interpolated,
                                         quantile(avg_chl_interpolated, PERCENTILE_SQUISHING_INTERVAL, na.rm = T)))

rm(PERCENTILE_SQUISHING_INTERVAL)
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
# - avg_chl_interpolated: climatologia interpolata con mice e "squishata" nell'intervallo indicato.

#-------------------------------------------------------------------------------
# Calculate useful indeces over the interpolated climatology (avg_chl_interpolated)

print("Calculating useful indeces...")

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

source(file.path(AUX_FUNCTIONS_PATH, "plot_calculated_indeces.R"))

# Plot pixel 1730 indeces
plot_calculated_indeces(1730)

#-------------------------------------------------------------------------------
# Find zero points

print("Finding zero points and blooms...")

# Load useful functions
source(file.path(AUX_FUNCTIONS_PATH, "find_zero_points.R"))
source(file.path(AUX_FUNCTIONS_PATH, "find_number_of_blooms.R"))
source(file.path(AUX_FUNCTIONS_PATH, "build_table_zero_points_1.R"))

# Find zero points of each climatology
zero_pts_raw <- find_zero_points(climatology)

# Find number of blooms
n_blooms_raw <- find_number_of_blooms(zero_pts_raw)

# L'intersezione con lo zero si trova indicando due punti quello prima e dopo l'intersezione.
# -> N punti di zero = n_bloom * 4/2
# -> Punti trovati = n_bloom * 4
# -> -> Numero di bloom trovati = Punti trovati / 4

# Barplot of frequency of number of blooms found
table(n_blooms_raw)
barplot(table(n_blooms_raw),
        main = "Frequency of number of blooms found",
        col = "blue")

# The calculated data is arranged in a dataframe
zero_points_df <- build_table(zero_pts_raw, n_blooms_raw) %>%
    # Flag pixels with 3 or more blooms
    mutate(flagged = n_blooms >= 3,
           id_week_zero_crossing = ceiling(id_date_zero_crossing / 7)) %>%
    arrange(id_pixel, id_date_zero_crossing)

#rm(n_blooms, zero_pts)
#-------------------------------------------------------------------------------
# Content of zero_points_df:

# - id_pixel
# - id_date_zero_crossing: sequential points before and after zero-crossing found.
# - n_blooms: number of blooms for this pixel
# - flagged: TRUE if for this pixel the number of blooms found is >= 3.

#-------------------------------------------------------------------------------
# Increase resolution of chl from one observation every 8 days to 1 observation per day
# The increase in resolution is obtained by interpolating with Stineman algorithm

#############
# Resolution increase is done only for pixels with less than 3 blooms!!!
#############

print("Interpolating data to increase time resolution from 8 days to 1 day...")

# Extract list of valid pixels (i.e. pixels with less than 3 blooms)
unique_valid_pixels <- zero_points_df %>%
    filter(flagged == FALSE) %>%
    select(id_pixel) %>%
    distinct() %>%
    pull()

# New climatology dataframe with new time axis (with increased time resolution)
climatology_high_res <- data.frame(id_pixel = rep(as.numeric(unique_valid_pixels), each = 328),
                       id_date_extended = rep(1:328, length(unique_valid_pixels))) %>%
    as_tibble() %>%
    arrange(id_pixel, id_date_extended)

# Add known chl values to new climatology time axis
climatology_high_res <- climatology %>%
    select(id_pixel,
           id_date,
           avg_chl_interpolated,
           D_mav) %>%
    filter(id_pixel %in% unique_valid_pixels) %>%
    left_join(climatology_high_res, ., by = c("id_pixel", "id_date_extended" = "id_date"))

# Interpolate using splines (by pixel)
climatology_high_res <- climatology_high_res %>%
    group_by(id_pixel) %>%
    mutate(avg_chl_interpolated_high_res = imputeTS::na.interpolation(avg_chl_interpolated, option="stine"),
           D_mav_high_res_from_stine = imputeTS::na.interpolation(D_mav, option="stine")) %>%
    ungroup()

#-------------------------------------------------------------------------------
# Check interpolation quality on a random pixel

# set.seed(10)
# n <- sample(1:length(unique_valid_pixels), size = 1)
# #n <- 480 # Pixel 2624

source(file.path(AUX_FUNCTIONS_PATH, "actual_vs_interpolated_plots.R"))
compare_data_interpolation(2624)

#-------------------------------------------------------------------------------
# Trovare I PUNTI DI ZERO SULLA MOVING AVERAGE INTERPOLATA, OSSIA SU D_mav_high_res_from_stine

# Find zero points of each climatology. SU D_MAV_HIGH_RES???
zero_pts <- find_zero_points(climatology_high_res, "id_date_extended", D_mav_name = "D_mav_high_res_from_stine")

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
zero_points_df_high_res <- build_table(zero_pts, n_blooms) %>%
    mutate(flagged = n_blooms >= 3,
           # Convert id_date to week
           id_week_zero_crossing = ceiling(id_date_zero_crossing / 7)) %>%
    # Flag pixels with 3 or more blooms
    arrange(id_pixel, id_date_zero_crossing)

#rm(n_blooms, zero_pts)
#-------------------------------------------------------------------------------

