# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# S4 implementation of generalized iLUCK models
# -- Classes for inference on data from a Binomial distribution 
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Class for model: BinomialLuckModel
# ---------------------------------------------------------------------------- #

# class definition: extends class LuckModel
setClass ("BinomialLuckModel",
          contains = "LuckModel")

# ---------------------------------------------------------------------------- #
# Constructor function for BinomialLuckModel
# ---------------------------------------------------------------------------- #
# arguments may either be the arguments for LuckModel()
# or an object of class LuckModel()
BinomialLuckModel <- function (arg1 = NULL, n0 = NULL, y0 = NULL, data = new("BinomialData")){
  if (all(is.null(c(arg1, n0, y0)))) {
    stop("No arguments given for BinomialLuckModel()!")
  } else {
    # .data will stay a default object if no data argument was given.
    #.data <- BinomialData(data)
    .data <- data
    # if arg1 is not a LuckModel, see if n0 and y0 are given.
    # If so, read out n0 and y0.
    if (!is(arg1, "LuckModel")){
      if (any(is.null(c(n0, y0))))
        stop("To define a BinomialLuckModel both arguments n0 and y0 must be given!")
      .n0 <- n0
      .y0 <- y0
    } else { # if arg1 is a LuckModel, read out its slots
      .n0 <- n0(arg1)
      .y0 <- y0(arg1)
      if (!is.null(tauN(data(arg1)))) { # LuckModel contains data    
        if (!is.null(tauN(.data))) # data was given also in the constructor call
          stop("Object to be converted contains already some data.\n Replace this data later if it should be changed.")
        # overwrite default data object with data slot from arg1
        #.data <- BinomialData(data(arg1))
        .data <- data
      }
    }
    # build a LuckModel object from the slots for checking inputs
    .object <- LuckModel(n0 = .n0, y0 = .y0, data = .data)
    # check if y0 is one-dimensional, n0 > 1 and y0 > 0
    if (dim(y0(.object))[1] != 1)
      stop("For inference on Binomial data y0 must be one-dimensional!")
    if (y0(.object)[1] < 0 || y0(.object)[2] > 1)
      stop("For inference on Binomial data, y0 must be in [0, 1]!")
    # now create the WeibullLuckModel object from the LuckModel object
    new("BinomialLuckModel", n0 = n0(.object), y0 = y0(.object), data = .data)
  }
}

# accessor and replacement methods are inherited from LuckModel

# ---------------------------------------------------------------------------- #
# show (= print in S3) method for BinomialLuckModel objects                                                             
# ---------------------------------------------------------------------------- #

# helper function to produce Binomial specific text 
.pItextBin <- function(object) {
  .varenumerator <- function(.y) .y*(1-.y)
  .mirroredy <- function(.ylu) {
    if(.ylu[2] <= 0.5) {
      return(.ylu)
    } else {
      return(c(min(.ylu[1], 1 - .ylu[2]), max(.ylu[1], 1 - .ylu[2])))
    }
  } 
  .minvarenumerator <- function(.ylu) .varenumerator(.mirroredy(.ylu)[1])
  .maxvarenumerator <- function(.ylu) .varenumerator(.mirroredy(.ylu)[2])
  .oneNoneY = paste("corresponding to a Beta prior\n with mean", #
                    y0(object)[1], #                            
                    "and variance", #                             
                    (y0(object)[1]*(1 - y0(object)[1]))/(n0(object)[1] + 1))
  .oneNtwoY = paste("corresponding to a set of Beta priors\n with means in [", #
                    y0(object)[1], ";", y0(object)[2], #                    
                    "] and variances in [", #                                         
                    .minvarenumerator(y0(object))/(n0(object)[1] + 1), ";", #
                    .maxvarenumerator(y0(object))/(n0(object)[1] + 1), "]")    
  .twoNoneY = paste("corresponding to a set of Beta priors\n with mean", #
                    y0(object)[1], #                                    
                    "and variances in [", #                               
                    .varenumerator(y0(object)[1])/(n0(object)[2] + 1), ";", #
                    .varenumerator(y0(object)[1])/(n0(object)[1] + 1), "]")    
  .twoNtwoY = paste("corresponding to a set of Beta priors\n with means in [", #                                   
                    y0(object)[1], ";", y0(object)[2], #                    
                    "] and variances in [", #                                   
                    .minvarenumerator(y0(object))/(n0(object)[2] + 1), ";", #
                    .maxvarenumerator(y0(object))/(n0(object)[1] + 1), "]")    
  .genilucktree(object = object, oneNoneY = .oneNoneY, #
                oneNtwoY = .oneNtwoY, #
                twoNoneY = .twoNoneY, #
                twoNtwoY = .twoNtwoY)
}

.showLuckModels <- function (object, forInference = "", parameterInterpretation = NULL) {
  # define elementary text modules
  .geniluck  <- function (forInference)
    paste("generalized iLUCK model ", forInference, "with prior parameter set:", sep = "")
  .iluck     <- function (forInference)
    paste("iLUCK model ", forInference, "with prior parameter set:", sep = "")
  .luck      <- function (forInference)
    paste("LUCK model ", forInference, "with prior parameters:", sep = "")
  .twoN      <- paste ("\n  lower n0 =", n0(object)[1],"  upper n0 =", n0(object)[2])
  .oneN      <- paste ("\n  n0 =", n0(object)[1])
  .ydim      <- dim(y0(object))[1]
  .twoY      <- function (.dim) {
    if (.ydim == 1) paste ("\n  lower y0 =", y0(object)[.dim,1],"  upper y0 =", y0(object)[.dim,2])
    else paste ("\n  lower y0 =", y0(object)[.dim,1],"  upper y0 =", y0(object)[.dim,2], " for dimension ", .dim)
  }
  .oneY      <- function (.dim) {
    if (.ydim == 1) paste ("\n  y0 =", y0(object)[.dim,1])
    else paste ("\n  y0 =", y0(object)[.dim,1], " for dimension ", .dim)
  }
  .imprec    <- function (.dim) {
    if (.ydim == 1) paste ("\n giving a main parameter prior imprecision of", y0(object)[.dim,2]-y0(object)[.dim,1])
    else paste ("\n giving a main parameter prior imprecision of", y0(object)[.dim,2]-y0(object)[.dim,1], " for dimension", .dim)
  }
  # data switch
  .hasdata <- !is.null(tauN(data(object)))
  # parameterInterpretation switch
  .pI <- !is.null(parameterInterpretation)
  # arguments for the case differentiation function
  .oneNoneY <- function (...) { 
    if (!.hasdata) {                   # without data
      cat(.luck(forInference), .oneN, .oneY(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
    } else {                           # with data
      cat(.luck(forInference), .oneN, .oneY(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
      cat("and ")
      show(data(object))
    } 
  }
  .oneNtwoY <- function (...) {
    if (!.hasdata) {                   # without data
      cat(.iluck(forInference), .oneN, .twoY(1:.ydim), .imprec(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
    } else {                           # with data
      cat(.iluck(forInference), .oneN, .twoY(1:.ydim), .imprec(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
      cat("and ")
      show(data(object))
    }
  }
  .twoNoneY <- function (...) {
    if (!.hasdata) {                   # without data
      cat(.geniluck(forInference), .twoN, .oneY(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
    } else {                           # with data
      cat(.geniluck(forInference), .twoN, .oneY(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
      cat("and ")
      show(data(object))
    }
  }
  .twoNtwoY <- function (...) {
    if (!.hasdata) {                   # without data
      cat(.geniluck(forInference), .twoN, .twoY(1:.ydim), .imprec(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
    } else {                           # with data
      cat(.geniluck(forInference), .twoN, .twoY(1:.ydim), .imprec(1:.ydim), "\n")
      if (.pI) cat(parameterInterpretation, "\n")
      cat("and ")
      show(data(object))
    }
  }
  # use the case differentiation function
  # to produce the output
  .genilucktree(object = object, oneNoneY = .oneNoneY(), oneNtwoY = .oneNtwoY(), #
                twoNoneY = .twoNoneY(), twoNtwoY = .twoNtwoY())
}

.genilucktree <- function (object, oneNoneY, oneNtwoY, twoNoneY, twoNtwoY) {
  if (n0(object)[1] == n0(object)[2]) {         # single n
    if (all(y0(object)[,1] == y0(object)[,2])) { # single n, single y
      oneNoneY
    } else {                                     # single n, interval-valued y
      oneNtwoY
    } 
  } else {                                      # interval-valued n
    if (all(y0(object)[,1] == y0(object)[,2])) { # interval-valued n, single y
      twoNoneY
    } else {                                     # interval-valued n, interval-valued y
      twoNtwoY
    }
  }
}

# show method uses helper function .showLuckModels (see show method for LuckModel)
setMethod("show", "BinomialLuckModel", function(object){
  .showLuckModels(object = object, #
                  forInference = "for inference from Binomial data\n", #
                  parameterInterpretation = .pItextBin(object))
})

# Tests:
#BinomialLuckModel(n0=1, y0=0.1, data=BinomialData(n=10,s=3))

# ---------------------------------------------------------------------------- #
# singleHdi: function returning symmetric credibility intervals for BinomialLuckModel
# ---------------------------------------------------------------------------- #
# shape1 = n0*y0 and 
# shape2 = n0*(1-y0)

if(!isGeneric("singleHdi"))
  setGeneric("singleHdi", function(object, ...) standardGeneric("singleHdi"))

# argument object is used only for method dispatching
setMethod("singleHdi", "BinomialLuckModel", function(object, n, y, gamma) {
  # TODO: proper highest density interval instead of symmetric interval
  # using the quantiles for symmetric credibility interval
  .lqua <- (1-gamma)/2    # lower quantile
  .uqua <- gamma + .lqua  # upper quantile  
  .lhd <- qbeta (.lqua, shape1 = n*y, shape2 = n*(1-y)) # lower border
  .uhd <- qbeta (.uqua, shape1 = n*y, shape2 = n*(1-y)) # upper border
  c(.lhd,.uhd) # return the credibility interval
})


# ---------------------------------------------------------------------------- #
# singleCdf: function returning cdf values for BinomialLuckModel
# ---------------------------------------------------------------------------- #

if(!isGeneric("singleCdf"))
  setGeneric("singleCdf", function(object, ...) standardGeneric("singleCdf"))

# argument object is used only for method dispatching
setMethod("singleCdf", "BinomialLuckModel", function(object, n, y, x) {
  pbeta (x, shape1 = n*y, shape2 = n*(1-y))
})

#Tests:

asdf <- BinomialLuckModel(n0=c(1,2), y0=c(0.4,0.6), data=BinomialData(n=2, s=1))
cdfplot(asdf)
unionHdi(asdf)$borders
cdfplot(asdf, control=controlList(posterior=TRUE))
unionHdi(asdf, posterior=TRUE)$borders


#