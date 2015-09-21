##' @name ChibML
##' @title Marginal Likelihood by Chib & Jeliazikov's Method for User-Written Functions
##' @description The function computes the marginal likelihood, i.e. the posterior normalising constant, with the method of Chib & Jeliazikov (2001) for user-written functions, from which an MCMC posterior sample is available.
##' @usage ChibML(logfun, theta.star, tune, V, mcmcsamp, df, verbose)
##' @param logfun The logarithm of the objective function
##' @param theta.star The starting value of the inner MCMC sampling and the value required by the Chib & Jeliazikov's method. This must be a high denstiy point, such as the posterior mean, median or mode.
##' @param tune The tunning value to be used to achieve the desired efficiency
##' @param V The proposal scale matrix
##' @param mcmcsamp The MCMC sample from the joint posterior
##' @param df The degrees of freedom of the proposal
##' @param verbose A switch which determines whether or not the progress of the sampler is printed to the screen. If verbose is greater than 0 the iteration number, and the Metropolis acceptance rate are sent to the screen every \code{verbose}th iteration
##' @details The function produce an approximation of the posterior normalizing constant via the Chib & Jeliazikov method in a single block sampling. The proposal distribution for the block is a Student's \eqn{t}-density with \code{df} degrees of freedom. The proposal is centered at the current value of \eqn{\theta}{theta} and has scale matrix \eqn{H}. \eqn{H} is calculated as: \eqn{H = TVT}{H = T*V*T}, where \eqn{T}{T} is a the diagonal positive definite matrix formed from the \code{tune}.
##'
##' @return double, the logarithm of the posterior normalising constant
##'
##' @references
##' Chib S. & Jeliazikov I. (2001).
##' Marginal likelihood from the Metropolis-Hastings output.
##' \emph{Journal of the American Statistical Association}, \bold{46},
##'  270--281.
##'
##' Robert C. P. & Casella G. (2004).
##'  \emph{Monte Carlo Statistical Methods}. 2nd Edition. New York: Springer.
##'
##' @seealso \code{\link[iLaplaceExamples]{nlpost_gomp}} and \code{\link[iLaplaceExamples]{nlpost_bod2}} for examples; \code{\link[iLaplaceExamples]{MHmcmc}}, \code{\link[iLaplaceExamples]{isML}}
##'
##' @rdname ChibML
##' @export
ChibML <- function(logfun, theta.star, tune, V, mcmcsamp, df, verbose){

  my.env <- environment(fun = logfun)

  tune <- vector.tune(tune, length(theta.star))

  if (nrow(V) != ncol(V) || nrow(V) != length(theta.star)) {
    cat("V not of appropriate dimension.\n")
    stop("Check V and theta.star and call ChibML() again. \n",
         call. = FALSE)
  }
  CC <- NULL
  try(CC <- chol(V), silent = TRUE)
  if (is.null(CC)) {
    cat("V not positive definite.\n")
    stop("Check V and call ChibML() again. \n",
         call. = FALSE)
  }
  V <- tune %*% V %*% tune

  ans = .Call('iLaplaceExamples_mlChib_cpp',
              PACKAGE = 'iLaplaceExamples',
              logfun,
              theta.star,
              V,
              mcmcsamp,
              df,
              verbose,
              my.env)
  return(logfun(theta.star) - ans)
}

##' @name isML
##' @title Marginal Likelihood by Importance Sampling for User-Written Functions
##' @description The function computes the marginal likelihood by importance sampling and from user-written functions.
##' @usage isML(logfun, nsim, theta.hat, tune, V, df, verbose)
##'
##' @param logfun The logarithm of the objective function
##' @param nsim The number of draws form the importance density
##' @param theta.hat The center of the proposal
##' @param tune A tunning value to be used to achieve the desired efficiency
##' @param V The scale matrix of the importance density
##' @param df The degrees of freedom of the importance density
##' @param verbose A switch which determines whether or not the progress of the sampler is printed to the screen. If verbose is greater than 0 the iteration number, and importance sampling approximation are sent to the screen every \code{verbose}th iteration
##'
##' @return double, the logarithm of the marginal likelihood
##'
##' @examples
##' 
##'\dontrun{
##'}
##'
##'
##' @references
##' Chib S. & Jeliazikov I. (2001).
##' Marginal likelihood from the Metropolis-Hastings output.
##' \emph{Journal of the American Statistical Association}, \bold{46},
##'  270--281.
##'
##' Robert C. P. & Casella G. (2004).
##'  \emph{Monte Carlo Statistical Methods}. 2nd Edition. New York: Springer.
##'
##' @seealso \code{\link[iLaplaceExamples]{nlpost_gomp}}, \code{\link[iLaplaceExamples]{nlpost_bod2}} for examples; \code{\link[iLaplaceExamples]{MHmcmc}}, \code{\link[iLaplaceExamples]{ChibML}}
##'
##' @rdname isML
##' @export
isML <- function(logfun, nsim = 10000, theta.hat, tune, V, df = 5, verbose = 1000){
  my.env <- environment(fun = logfun)

  tune <- vector.tune(tune, length(theta.hat))
  if (nrow(V) != ncol(V) || nrow(V) != length(theta.hat)) {
    cat("V not of appropriate dimension.\n")
    stop("Check V and theta.star and call ChibML() again. \n",
         call. = FALSE)
  }
  CC <- NULL
  try(CC <- chol(V), silent = TRUE)
  if (is.null(CC)) { 
    cat("V not positive definite.\n")
    stop("Check V and call ChibML() again. \n",
         call. = FALSE)
  }
  V <- tune %*% V %*% tune
  p <- ncol(V)
  ans = .Call('iLaplaceExamples_mlIS_cpp',
              PACKAGE = 'iLaplaceExamples',
              logfun,
              theta.hat,
              V,
              nsim,
              p,
              df,
              verbose,
              my.env)
  return(ans)
}
