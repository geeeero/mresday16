# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# S4 implementation of generalized iLUCK models
# -- Classes for inference on data from a Binomial distribution 
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Class for data: BinomialData
# ---------------------------------------------------------------------------- #

# data class definition
setClass ("BinomialData",
          contains = "LuckModelData")

# constructor function handling as input
# the number of successes s and the sample size n
BinomialData <- function(s, n){
  if(any(floor(c(s,n)) != c(s,n)) || s < 0 || n < 0)
    stop("s and n must be non-negative integers!")
  if(s > n)
    stop("s must be smaller or equal to n!")
  arg1 <- LuckModelData(tau=s, n=n)
  return (new("BinomialData", tauN = tauN(arg1), rawData = rawData(arg1)))
}

# show method for class BinomialData (no plot method)
setMethod("show", "BinomialData", function(object){
  if (is.null(tauN(object))) {
    cat("Default data object with no data specified.\n")
  } else {
    .rawData <- rawData(object)
    .tau <- tau(object)
    .n <- n(object)
    if (is.null(.rawData)) {
      cat("BinomialData object containing number of successes s = ", .tau, " and sample size ", .n, ".\n", sep="")
    } else { # object has rawData
      cat("BinomialData object containing ", .n, " observations: \n",
          as.vector(.rawData), "\n", sep="")
    }
  }
})

# Tests:
BinomialData(s=5,n=10)
#BinomialData(s=5.5,n=-10) # error
#BinomialData(s=5,n=-10)   # error
#