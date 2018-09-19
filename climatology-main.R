#-------------------------------------------------------------------------------
# Bloom finder script V. 1.0

#-------------------------------------------------------------------------------
# Clean workspace
rm(list = ls())

# Set seed for reproducibility purposes
set.seed(332423)

#-------------------------------------------------------------------------------
# Load required libraries

require(qchlorophyll) # 2.1
require(dplyr)        # 0.7.4
require(scales)       # 0.5.0
# TTR                 # 0.23-3
# lazyeval            # 0.2.1
# log4r               # 0.2
# readr               # 1.1.1

#-------------------------------------------------------------------------------
# Set global parameters

# Maximum number of allowed consecutive NAs in the climatology for each pixel
MAX_CONSECUTIVE_NA <- 2
# Threshold percentage to calculate s = median(chl) * (1 + THRESHOLD_PERCENTAGE)
THRESHOLD_PERCENTAGE <- 0.05
# Number of observations to be used in the moving average of C (cumulative sum of anomalies)
RUNNING_AVERAGE_WINDOW <- 3
# Percentile interval into which raw chl data needs to be squished into
PERCENTILE_SQUISHING_INTERVAL <- c(0.05, 0.95)
# Average operation to be used for climatology calculation. Set to "geom" for geometric mean
MEAN_FUNCTION <- "mean"
# Minimum duration in days for a suspected bloom to be considered an actual bloom
MINIMUM_BLOOM_DURATION_DAYS <- 16
# Pixels with number of blooms >= N_BLOOM_MAX are flagged in TABELLA_TRE
N_BLOOM_MAX <- 3

#-------------------------------------------------------------------------------
# Set paths

# Path of .nc files. This path must point to the folder containing the .nc files to be used in the analysis
NC_FILES_PATH <- "C:\\users\\michy\\desktop\\christian_paper\\DATA_TEST"
# Auxiliary functions path. This path must point to the folder "auxiliary_functions"
AUX_FUNCTIONS_PATH <- "C:\\users\\michy\\desktop\\christian_paper\\SCRIPT\\auxiliary_functions"
# Output directory. This path can point to whatever folder you wish
OUTPUT_PATH <- "C:\\users\\michy\\desktop"

#-------------------------------------------------------------------------------
# Setup of logger

logger <- log4r::create.logger(logfile = file.path(OUTPUT_PATH, "log_climatology_main.log"),
                               level = "INFO")

#-------------------------------------------------------------------------------
# Load data

print("Starting to load data...")

# Load .nc files and extract CHL1_mean
nc_files_list <- load_all_as_list(path = NC_FILES_PATH, variables = c("CHL1_mean"))
# Bind all observations in a single dataframe
nc_dataframe <- assign_id_and_melt(nc_files_list)

rm(nc_files_list, NC_FILES_PATH)
#-------------------------------------------------------------------------------
# Squish of raw chl data in the selected percentile interval

print("Squishing climatology in the selected percentile interval...")

nc_dataframe <- nc_dataframe %>%
    group_by(id_pixel, id_date) %>%
    mutate(CHL1_mean = squish(CHL1_mean,
                              quantile(CHL1_mean,
                                       PERCENTILE_SQUISHING_INTERVAL,
                                       na.rm = T))) %>%
    ungroup()

rm(PERCENTILE_SQUISHING_INTERVAL)
#-------------------------------------------------------------------------------
# Calculate climatology on RAW chl data

print("Calculating climatology...")

# Load function to calculate consecutive NAs in climatology
source(file.path(AUX_FUNCTIONS_PATH, "consecutive_na_count.R"))

# Number of pixels used as input
pixels_IN <- length(unique(nc_dataframe$id_pixel))

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

# Calculate climatology
climatology <- nc_dataframe %>%
    # For each pixel, then for each date
    group_by(id_pixel, id_date) %>%
    # Calculate: climatology, i.e. average value for date, for pixel (avg_chl)
    #           number of observations used in each date to calculate avg_chl (n_observations_used_per_date)
    #           number of dates loaded in R (images per pixel loaded)
    summarise_(.dots = setNames(list(lazyeval::interp( ~ MEAN_FUNCTION(CHL1_mean)),
                                     lazyeval::interp( ~ sum(!is.na(CHL1_mean))),
                                     lazyeval::interp( ~ n()) ),
                                c("avg_chl",
                                  "n_observations_used_per_date",
                                  "n_observations_total"))) %>%
    # Calculate for each pixel: how many missing data are in the climatology (NA_in_climatology_per_pixel) 
    #                           should the pixel be kept? (keep_pixel_NA_consecutive: TRUE -> keep, FALSE -> drop)
    mutate(NA_in_climatology_per_pixel = sum(is.na(avg_chl)),
           keep_pixel_NA_consecutive = pixel_consecutive_NA(avg_chl, n = MAX_CONSECUTIVE_NA)) %>%
    # Remove grouping by pixel
    ungroup() %>%
    # If the number of consecutive NAs is > n then that pixel must be removed from analysis
    # Keep only pixels with less than n consecutive NAs
    filter(keep_pixel_NA_consecutive == TRUE) %>%
    # Remove unused columns
    select(-keep_pixel_NA_consecutive)

# Number of pixels kept for the analysis
pixels_kept <- length(unique(climatology$id_pixel))

# Percentage of pixels kept for the analysis on the total of pixels used as input
print(paste("Percentage of pixels kept: ", signif(pixels_kept/pixels_IN * 100, 4), " %", sep = ""))

# Log the information
log4r::info(logger, paste("Pixels used as input: ", pixels_IN, sep = ""))
log4r::info(logger, paste("Pixels used for analysis: ", pixels_kept, sep = ""))
log4r::info(logger, paste("Percentage of pixels used for analysis over input: ", signif(pixels_kept/pixels_IN * 100, 4), " %", sep = ""))

rm(pixel_consecutive_NA, MAX_CONSECUTIVE_NA, pixels_IN, pixels_kept, MEAN_FUNCTION)
#-------------------------------------------------------------------------------
# Add info on longitude and latitude (lon and lat)

print("Adding info on longitude and latitude...")

climatology <- nc_dataframe %>%
    select(lon, lat, id_pixel) %>%
    distinct() %>%
    left_join(climatology, ., by = "id_pixel")

# Keep only a dataframe of unique pixels with their coordinates
nc_dataframe <- nc_dataframe %>%
    select(id_pixel, lon, lat) %>%
    distinct()

#-------------------------------------------------------------------------------
# Imputation of missing data in climatology using linear interpolation

print("Imputing missing data...")

climatology <- climatology %>%
    group_by(id_pixel) %>%
    mutate(avg_chl_interpolated = imputeTS::na.interpolation(avg_chl, option="linear")) %>%
    ungroup()

#-------------------------------------------------------------------------------
# Calculate useful indeces over the climatology without missing values (avg_chl_interpolated)

print("Calculating useful indeces...")

# Calcola s, A, D, D_mav (moving average)
climatology <- climatology %>%
    # Already ungrouped.
    group_by(id_pixel) %>%
    # For each pixel, calculate: s = median + a% of median of avg_chl_interpolated
    #                           Anomalies = avg_chl_interpolated - s
    #                           C = cumulative sum of A
    #                           D = time derivative of C
    #                           D_mav = moving average (or running average) of D
    mutate(s = median(avg_chl_interpolated)*(1 + THRESHOLD_PERCENTAGE),
           A = avg_chl_interpolated - s,
           C = cumsum(A),
           D = (C - dplyr::lag(C)) / 8, # 8 days observations
           D_mav = TTR::runMean(D, RUNNING_AVERAGE_WINDOW)) %>%
    # Remove grouping by pixel
    ungroup()

source(file.path(AUX_FUNCTIONS_PATH, "plot_calculated_indeces.R"))

# Plot pixel 1729 indeces
# plot_calculated_indeces(1729)

rm(RUNNING_AVERAGE_WINDOW, THRESHOLD_PERCENTAGE)
#-------------------------------------------------------------------------------
# Content of dataframe climatology:

# - id_pixel: unique pixel id
# - id_date: day of the year of data with frequency of 1 observation every 8 days
# - avg_chl: climatology on raw data, that is, chl value averaged over the available years by pixel and by date on raw data
# - n_observations_used_per_date: number of observations used for calculating climatology for each date of a given pixel
# - n_observations_total: number of images per pixel loaded in R
# - NA_in_climatology_per_pixel: number of NAs in climatology (avg_chl) for a given pixel
# - lon: longitude
# - lat: latitude
# - avg_chl_interpolated: climatology avg_chl linearly interpolated (to impute missing values)
# - s: median + a % of median of avg_chl_interpolated
# - A: avg_chl_interpolated - s
# - C: cumulative sum of A
# - D: time derivative of C
# - D_mav: Moving average of D

#-------------------------------------------------------------------------------
# Positive slope check

# Since later in the script the following assumption is made:
## Strong hypothesis: D_mav has positive slope in the first zero point
# a check is performed to analyse only those pixels for which the hypothesis hold

# Load function to do the check
source(file.path(AUX_FUNCTIONS_PATH, "check_slope.R"))
# Pixel checked and that will be further processed
pixel_checked <- check_slope()
# Pixel discarded since they do not satisfy hypothesis
pixel_discarded <- unique(climatology$id_pixel)[!unique(climatology$id_pixel) %in% pixel_checked]
# Log discarded pixels
log4r::warn(logger, paste("Pixels discarded due to not satisfying hypothesis: ", pixel_discarded, sep = ""))

# Filter climatology according to checked pixels
climatology <- climatology %>%
    filter(id_pixel %in% pixel_checked)

rm(check_slope, pixel_checked, pixel_discarded)
#-------------------------------------------------------------------------------
# Find zero points on climatology with frequency of 1 observation every 8 days

print("Finding zero points and blooms...")

# Load useful functions
source(file.path(AUX_FUNCTIONS_PATH, "find_zero_points.R"))
source(file.path(AUX_FUNCTIONS_PATH, "find_number_of_blooms.R"))
source(file.path(AUX_FUNCTIONS_PATH, "build_table_zero_points_1.R"))

# Find zero points of each climatology for each pixel
zero_pts_raw <- find_zero_points(climatology)

# Find number of blooms
n_blooms_raw <- find_number_of_blooms(zero_pts_raw)

# NOTE:
# 1. L'intersezione con lo zero si trova indicando due punti: quello prima e quello dopo il cambio di segno del segnale.
# 2. Ogni bloom è caratterizzato da 2 punti di zero.
# 3. Il numero di punti di zero trovati è quindi: n_bloom * 2
# 4. Il numero di punti trovati è quindi: n_bloom * 4
# 5. Il numero di bloom trovati è quindi: n_bloom = Punti trovati / 4

# Barplot of count of number of blooms found
print("Count of number of blooms found")
table(n_blooms_raw)
barplot(table(n_blooms_raw),
        main = "Count of number of blooms found",
        col = "blue")

# The calculated data is arranged in a dataframe
zero_points_df <- build_table(zero_pts_raw, n_blooms_raw) %>%
    # Flag pixels with 3 or more blooms
    mutate(flagged = n_blooms >= N_BLOOM_MAX,
           id_week_zero_crossing = ceiling(id_date_zero_crossing / 7)) %>%
    arrange(id_pixel, id_date_zero_crossing)

# Pixels with a number of blooms equal or higher than 3. Flagged for later insertion in TABELLA_TRE
flagged_pixels <- zero_points_df %>%
    filter(flagged == TRUE) %>%
    select(id_pixel) %>%
    distinct() %>%
    pull()

rm(n_blooms_raw, zero_pts_raw)
#-------------------------------------------------------------------------------
# Content of zero_points_df:

# - id_pixel
# - id_date_zero_crossing: sequential points before and after zero-crossing found
# - n_blooms: number of blooms found for this pixel
# - flagged: TRUE if for this pixel the number of blooms found is >= N_BLOOM_MAX

#-------------------------------------------------------------------------------
# Increase resolution of chl from one observation every 8 days to 1 observation per day

# 1. The increase in resolution is obtained by interpolating with Stineman algorithm
# 2. Interpolation is done on variables D_mav and avg_chl_interpolated

print("Interpolating data to increase time resolution from 8 days to 1 day...")

# Extract list of pixels to interpolate
unique_valid_pixels <- zero_points_df %>%
#    filter(flagged == FALSE) %>%
    select(id_pixel) %>%
    distinct() %>%
    pull()

print(paste("Interpolating ", length(unique_valid_pixels), " pixels...", sep = ""))

# New climatology dataframe with new time axis with finer time resolution
climatology_high_res <- data.frame(id_pixel = rep(unique_valid_pixels, each = 328),
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

# Interpolate (by pixel) using Stineman algorithm
climatology_high_res <- climatology_high_res %>%
    group_by(id_pixel) %>%
    mutate(avg_chl_interpolated_high_res = imputeTS::na.interpolation(avg_chl_interpolated, option="stine"),
           D_mav_high_res_from_stine = imputeTS::na.interpolation(D_mav, option="stine")) %>%
    ungroup()

rm(zero_points_df)
#-------------------------------------------------------------------------------
# Check interpolation quality on a random pixel

source(file.path(AUX_FUNCTIONS_PATH, "actual_vs_interpolated_plots.R"))
# compare_data_interpolation(2624)

#-------------------------------------------------------------------------------
# Find zero points and blooms on the high resolution moving average (i.e. on D_mav_high_res_from_stine)

print("Finding zero points and blooms on finer data...")

# Find zero points of each climatology
zero_pts <- find_zero_points(climatology_high_res,
                             id_date_name = "id_date_extended",
                             D_mav_name = "D_mav_high_res_from_stine")

# Find number of blooms
n_blooms <- find_number_of_blooms(zero_pts)

# Barplot of frequency of number of blooms found
table(n_blooms)
barplot(table(n_blooms),
        main = "Count of number of blooms found",
        col = "blue")

# The data is arranged in a dataframe
zero_points_df_high_res <- build_table(zero_pts, n_blooms) %>%
    # Convert id_date to week
    mutate(id_week_zero_crossing = ceiling(id_date_zero_crossing / 7)) %>%
    arrange(id_pixel, id_date_zero_crossing)

rm(n_blooms, zero_pts)
#-------------------------------------------------------------------------------
# Some pixels may have more blooms now due to interpolation process

# compare_data_interpolation(2807)
# plot_calculated_indeces(2807)

#-------------------------------------------------------------------------------
# Generate TABELLA_DUE

################################################################################
## Strong hypothesis: D_mav has positive slope in the first zero point
################################################################################

# Content of TABELLA_DUE
#
# - id_pixel: unique pixel id
# - bloom_duration_days
# - bloom_duration_weeks
# - bloom_start_date
# - bloom_start_week
# - bloom_end_date
# - bloom_end_week
# - n_blooms: number of blooms found for this pixel
# - flagged: TRUE if for this pixel the number of blooms found is >= 3.
# - lon: longitude
# - lat: latitude
# - max_chl: maximum value of chl during bloom
# - id_date_max_chl: corresponding id_date of max_chl

# Load function to find maximum
source(file.path(AUX_FUNCTIONS_PATH, "find_maximum_chl.R"))

TABELLA_DUE <- zero_points_df_high_res %>%
    
    # For each pixel
    group_by(id_pixel) %>%
    # Assign unique id to zero points
    mutate(id_zero = rep(1:n()/2, each=2, length.out = n())) %>%
    
    # For each pixel, for each zero point
    group_by(id_pixel, id_zero) %>%
    # Take the average of the 2 points which identify the zero point
    # This number (date/week) will IDENTIFY the zero point
    summarise(id_date_zero = mean(id_date_zero_crossing),
              id_week_zero = mean(id_week_zero_crossing)) %>%
    
    # For each pixel, assign a unique id to the blooms. The hypothesis is used here.
    mutate(id_bloom = rep(1:n(), each=2, length.out = n())) %>%
    
    # For each pixel, for each bloom
    group_by(id_pixel, id_bloom) %>%
    # Caclulate duration, start and end
    mutate(bloom_duration_days = id_date_zero - lag(id_date_zero),
           bloom_duration_weeks = id_week_zero - lag(id_week_zero),
           bloom_start_date = lag(id_date_zero),
           bloom_start_week = lag(id_week_zero),
           bloom_end_date = bloom_start_date + bloom_duration_days,
           bloom_end_week = bloom_start_week + bloom_duration_weeks) %>%
    
    # Keep only pixels with duration longer than MINIMUM_BLOOM_DURATION_DAYS
    filter(bloom_duration_days >= MINIMUM_BLOOM_DURATION_DAYS) %>%
    filter(!is.na(bloom_duration_days)) %>%
    
    # For each pixel add the number of blooms found
    group_by(id_pixel) %>%
    mutate(n_blooms = n())

# Remove unused columns
TABELLA_DUE <- TABELLA_DUE %>%
    select(-id_zero, -id_bloom, -id_date_zero, -id_week_zero) %>% 
    # Remove grouping by pixel
    ungroup()

# Add lon and lat to final df
TABELLA_DUE <- TABELLA_DUE %>%
    left_join(nc_dataframe, by = c("id_pixel"))

# Add max chl and date of max chl
TABELLA_DUE <- TABELLA_DUE %>%
    bind_cols(find_max_chl())

rm(zero_points_df_high_res, MINIMUM_BLOOM_DURATION_DAYS, find_max_chl)
#-------------------------------------------------------------------------------
# Generate TABELLA_TRE

# Content of TABELLA_TRE

# - id_pixel: unique pixel id
# - lon: longitude
# - lat: latitude
# - n_blooms: number of blooms found for this pixel
# - too_many_blooms: TRUE if for this pixel the number of blooms found is >= N_BLOOM_MAX.

TABELLA_TRE <- TABELLA_DUE %>%
    select(id_pixel, n_blooms) %>%
    distinct() %>%
    left_join(nc_dataframe, ., by = c("id_pixel")) %>%
    mutate(too_many_blooms = id_pixel %in% flagged_pixels)

rm(nc_dataframe)
#-------------------------------------------------------------------------------
# Save results

readr::write_csv(climatology, file.path(OUTPUT_PATH, "climatology.csv"))
readr::write_csv(climatology_high_res, file.path(OUTPUT_PATH, "climatology_high_res.csv"))
readr::write_csv(TABELLA_DUE, file.path(OUTPUT_PATH, "TABELLA_DUE.csv"))
readr::write_csv(TABELLA_TRE, file.path(OUTPUT_PATH, "TABELLA_TRE.csv"))
