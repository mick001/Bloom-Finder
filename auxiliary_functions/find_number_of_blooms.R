################################################################################
#' This function finds the number of blooms and returns a vector with the number
#' of blooms for each set of zero points in the list x
#'
#' @param x a list of points of zero crossing
#' @return a vector with the number of blooms for each set of zero crossing points.
#' 
find_number_of_blooms <- function(x)
{
    # For each bloom you have 2 zero crossings ->
    # For each zero crossing you have 2 points (before zero crossing and after crossing)
    # Therefore: number of blooms = number of points / 4
    
    # L'intersezione con lo zero si trova indicando due punti quello prima e dopo l'intersezione.
    # -> N punti di zero = n_bloom * 4/2
    # -> Punti trovati = n_bloom * 4
    # -> -> Numero di bloom trovati = Punti trovati / 4
    
    return( sapply(x, function(x1){ length(x1) / 4 }) )
}
