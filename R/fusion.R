#' @title P-value fusion
#'
#' @description Combine multiple p-values to a single p-value
#'
#' @param x Vector of p-values to combine
#'
#' @param alpha Alpha method to use for the Binomial method
#'
#' @param type Method to use for merging, fisher, binomial, bonferroni, harmonic
#'
#' @return Named list of combined p-values per method
#'
#' @keywords RBPInper, P-value combine
#'
#' @examples See the manuscript.
#'
#' @export
#'
fusion <- function(x = c(0.02,NA, NA), alpha = 0.05,
                   type = c("Fish", "Bino", "Bonf", "Harm")){

 # Replace NA's with 1, that no evidence
  x[is.na(x)] <- 1

  # Check that the p values are within 0 and 1
  if(sum(x < 0 | x > 1) >= 1){
    stop("Error: All p-values must be between 0 and 1") }

 # Check that the alpha level is between 0 and 1
 if(alpha < 0 | alpha > 1){
   stop("Error: alpha must be between 0 and 1") }

  # Check that only one type is requested
  if(length(type) != 1){
    stop("Error: type should only have one input") }

# Get length and check if it is greater than 1
k <- length(x)
if(k <= 1){
  output <- x
  return(output)
}


if(tolower(substr(type, 0,4)) == "fish"){
  # REF: DOI: 10.1007/978-1-4612-4380-9_6
  #****TAG: Independence is assumed****#
  fisher_stat <- -2 * sum(log(x))
  output <- pchisq(fisher_stat, df = 2 * k,
                     lower.tail = FALSE)

} else if(tolower(substr(type, 0,4)) == "bino"){
  # REF: DOI: 10.1037/h0059111
  bino_stat <- sum(x <= alpha)
  output <- sum(dbinom(bino_stat:k, k, alpha))

} else if(tolower(substr(type, 0,4)) == "bonf"){
  # REF: DOI: 10.1002/sim.6082
  #****TAG: No assumption of independence****#
  bon_stat <- min(x)
  output <- min(1, bon_stat * k)


} else if(tolower(substr(type, 0,4)) == "harm"){
  # REF: DOI: 10.1073/pnas.1814092116
  hmp_stat = 1/mean(1/x)
  output <- FMStable::pEstable(1/hmp_stat,
                               FMStable::setParam(alpha = 1,
                                location = (log(length(x)) + 0.874),
                                logscale = log(pi/2), pm = 0),
                               lower.tail = FALSE)

} else {
  stop("Error: only four methods are implemented/n
        including fisher, binomial, bonferroni, harmonic")
}

return(output)
}



