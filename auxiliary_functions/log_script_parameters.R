################################################################################
#' This function logs the main parameters of the script
#' 
#' @param logger a logger object
#' @return void
#' 
log_script_parameters <- function(logger)
{
    log4r::info(logger, "**** STARTING SCRIPT EXECUTION ****")
    
    log4r::info(logger, paste("MAX_CONSECUTIVE_NA: ",
                              MAX_CONSECUTIVE_NA,
                              sep = ""))
    
    log4r::info(logger, paste("THRESHOLD_PERCENTAGE: ",
                              THRESHOLD_PERCENTAGE,
                              sep = ""))
    
    log4r::info(logger, paste("Filter function used: ",
                              FILTER_NAME,
                              sep = ""))
    
    log4r::info(logger, paste("RUNNING_AVERAGE_WINDOW: ",
                              RUNNING_AVERAGE_WINDOW,
                              sep = ""))
    
    log4r::info(logger, paste("KERNEL: ",
                              KERNEL,
                              sep = ""))
    
    log4r::info(logger, paste("BANDWIDTH: ",
                              BANDWIDTH,
                              sep = ""))
    
    log4r::info(logger, paste("DF: ",
                              DF,
                              sep = ""))
    
    log4r::info(logger, paste("SPAN: ",
                              SPAN,
                              sep = ""))
    
    log4r::info(logger, paste("PERCENTILE_SQUISHING_INTERVAL: ",
                              paste(PERCENTILE_SQUISHING_INTERVAL, collapse = " "),
                              sep = ""))
    
    log4r::info(logger, paste("MEAN_FUNCTION: ",
                              MEAN_FUNCTION,
                              sep = ""))
    
    log4r::info(logger, paste("MINIMUM_BLOOM_DURATION_DAYS: ",
                              MINIMUM_BLOOM_DURATION_DAYS,
                              sep = ""))
    
    log4r::info(logger, paste("N_BLOOM_MAX : ",
                              N_BLOOM_MAX,
                              sep = ""))
    
    log4r::info(logger, paste("NEW_STARTING_POINT: ",
                              NEW_STARTING_POINT,
                              sep = ""))
    
    log4r::info(logger, paste("SPRING_FALL_ID_DATE_SEPARATOR: ",
                              SPRING_FALL_ID_DATE_SEPARATOR,
                              sep = ""))
    
    log4r::info(logger, paste("STARTING_YEAR: ",
                              STARTING_YEAR,
                              sep = ""))
    
    log4r::info(logger, paste("ENDING_YEAR: ",
                              ENDING_YEAR,
                              sep = ""))
    
    log4r::info(logger, paste("NC_FILES_PATH: ",
                              NC_FILES_PATH,
                              sep = ""))
    
    log4r::info(logger, paste("OUTPUT_PATH: ",
                              OUTPUT_PATH,
                              sep = ""))

}
