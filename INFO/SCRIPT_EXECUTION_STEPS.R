#-------------------------------------------------------------------------------
#          Execution steps performed by Bloom finder script V. 2.3             #
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# The script executes the following steps
#
# 0. A log file called "log_climatology_main.log" is saved in the specified OUTPUT_PATH.
#
# 1. Data is loaded from the specified data directory NC_FILES_PATH.
#
# 2. After grouping by pixel and by date, data is squished in the required percentile range.
#    Note that this operation is done by pixel and by date, that is, percentile range is calculated
#    on each set of observation for a particular date of a particular pixel.
#
# 3. Climatology is calculated using either arithmetic or geometric mean as specified in MEAN_FUNCTION.
#
# 4. Pixels with more than MAX_CONSECUTIVE_NA consecutive missing data in the climatology are removed from the analysis.
#
# 5. Missing data in climatology is obtained interpolating linearly by pixel.
#
# 6. For each pixel, the following time series are calculated on climatology: s, A, C, D, D_mav
#   
    #   s = median + a % of median
    #   Anomalies = climatology - s
    #   C = cumulative sum of anomalies
    #   D = time derivative of cumulative sum of anomalies
    #   D_mav = moving average (or running average of derivative)
#
# 7. Just before increasing the resolution of the data, the time axis is shifted so that
#    the time series can begin from a point defined by the user.
#    Then, the resolution of both D_mav and chl is increased from one observation every 8 days to 1 observation per day.
#    The increase in resolution is obtained by interpolation using Stineman's algorithm.
#
# 8. Zero points of D_mav and number of blooms of each pixel are searched on the higher resolution data.
#
# 9. Quantities of interest are calculated and stored in TABELLA_DUE
#
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
        # - lon: longitude
        # - lat: latitude
        # - max_chl: maximum value of chl during bloom
        # - id_date_max_chl: corresponding id_date of max_chl
        # - week_max_chl: corresponding week of the year of max_chl
#
# 10. Blooms that last less than the specified number of days in MINIMUM_BLOOM_DURATION_DAYS are removed from TABELLA_DUE.
#
# 11. Quantities of interest are calculated and stored in TABELLA_TRE
#
    # Content of TABELLA_TRE
    
        # - id_pixel: unique pixel id
        # - lon: longitude
        # - lat: latitude
        # - n_blooms: number of blooms found for this pixel
        # - too_many_blooms: TRUE if for this pixel the number of blooms found is >= N_BLOOM_MAX.
#
# 12. The following tables are saved in .csv format in OUTPUT_PATH:
#
        # TABELLA_DUE: results as specified in 11.
        # TABELLA_TRE: results as specified in 13.
        # climatology: low resolution climatology data and indeces (1 observation each 8 days of the year)
        # climatology_high_res: high resolution climatology data and indeces (1 observation per day of the year)
#        
    # Content of climatology_high_res
        
        # - id_pixel: unique pixel id
        # - id_date_extended: day of the year of data with frequency 1 observation every day
        # - avg_chl_interpolated: climatology avg_chl linearly interpolated (to impute missing values)
        # - D_mav: moving average of derivative of cumulative sum of anomalies (or running average of derivative)
        # - avg_chl_interpolated_high_res: avg_chl with higher resolution (1 observation every day)
        # - D_mav_high_res_from_stine: D_mav with higher resolution (1 observation every day)
#
    # Content of climatology
        
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
#        
