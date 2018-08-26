# This function finds the number of blooms
find_number_of_blooms <- function(x)
{
    # For each bloom you have 2 zero crossings ->
    # for each zero crossing you have 2 points (before zero crossing and after crossing).
    return( sapply(x, function(x1){ length(x1)/4 }) )
}
