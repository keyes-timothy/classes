# Function to calculate, adjust, and invert the covariance matrix for the reference populations
get_cov <- function(data) {
  Sx <- cov(data)
  Sx <- solve(Sx)
  return(Sx)
}
