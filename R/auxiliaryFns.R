dmvt <- function(x, mu, S, p, df, lg) {
  .Call('iLaplace_dmvt', PACKAGE = 'iLaplaceExtra', x, mu, S, p, df, lg)
}

dhalfCauchy <- function(x, scale, lg = FALSE) {
  .Call('iLaplace_dhalfCauchy', PACKAGE = 'iLaplaceExtra', x, scale, lg)
}

rmvnorm <- function(mu, S, p) {
  .Call('iLaplace_rmvnorm', PACKAGE = 'iLaplaceExtra', mu, S, p)
}

rmvt <- function(mu, S, p, df) {
  .Call('iLaplace_rmvt', PACKAGE = 'iLaplaceExtra', mu, S, p, df)
}

prJeff <- function(nu, lg = FALSE) {
  .Call('iLaplace_prJeff', PACKAGE = 'iLaplaceExtra', nu, lg)
}

grldmvt <- function(x, mu, S, p, df) {
  .Call('iLaplace_grldmvt', PACKAGE = 'iLaplaceExtra', x, mu, S, p, df)
}

hessldmvt <- function(x, mu, S, p, df) {
  .Call('iLaplace_hessldmvt', PACKAGE = 'iLaplaceExtra', x, mu, S, p, df)
}

## setup diagonal tuning matrix for vector parameters

# return name of the calling function
"calling.function" <-
  function(parentheses=TRUE) {
    calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
    if (parentheses){
      calling.function <- paste(calling.function, "()", sep="")
    }
    return(calling.function)
  }

"vector.tune" <- function(mcmc.tune, K){
  if (max(is.na(mcmc.tune))){
    cat("Error: Vector tuning parameter cannot contain NAs.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(mcmc.tune) == 1){
    mcmc.tune <- rep(mcmc.tune, K)
  }
  if (length(mcmc.tune) != K){
    cat("Error: length(vector tuning parameter) != length(theta) or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(mcmc.tune <= 0) != 0) {
    cat("Error: Vector tuning parameter cannot contain negative values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(mcmc.tune)==1){
    return(matrix(mcmc.tune, 1, 1))
  }
  else{
    return(diag(as.double(mcmc.tune)))
  }
}

# Given "q"-dimensional vector "par", "q"-vector "se", and scalars "lengthOut",
# "q" and "delta", this function creates for each of the "q" elements of
# "par" a grid of "lengthOut" values from "par-delta*se", to "par+delta*se". It
# returns a matrix of "lengthOut" times "q" values.
# "seqMat" <- function(par, se, lengthOut, q, delta) {
#   .Call('iLaplace_seqMat', PACKAGE = 'iLaplaceExtra', par, se, lengthOut, q, delta)
# }

# Given the "dimMat" times "dimMat" matrix "hessMat" (being positve definite),
# this function computes the log-determinant of all the diagonal blocks,
# starting from the whole matrix. The function returns a "dimMat"-vector, where
# the first element is the log-determinant of hessMat, the second element is
# the log-determinant of hessMat[-1,-1] and so on, except the last element which
# is simply log(hessMat[dimMat,dimMat])
# "ldetHessBlocks" <- function(hessMat, dimMat) {
#   .Call('iLaplace_ldetHessBlocks', PACKAGE = 'iLaplaceExtra', hessMat, dimMat)
# }

# Given the "dimMat" times "dimMat" matrix "hessMat" (being positve definite),
# this function computes the square root of the diagonal of inverse of blocks
# of "hessMat". The function returns a "dimMat" vector with the first element
# being sqrt(diag(solve(hessMat)))[1], the second element is sqrt(diag(solve(hessMat[-1,-1])))[1],
# and so on, except fo the last element which is sqrt(1/hessMat[dimMat,dimMat]).
# "SEv" <- function(hessMat, dimMat) {
#   .Call('iLaplace_SEv', PACKAGE = 'iLaplace', hessMat, dimMat)
# }
