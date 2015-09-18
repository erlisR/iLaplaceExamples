##' @name logitRegCRE
##' @title Logistic regression with crossed random effects with the Salamander data
##' @description  Log-likelihood functions and related quantities for a logistic regression model with crossed random effects based on the \code{\link[iLaplaceExtra]{Salamander}} data. In particular, the function \code{nlogH_salam} gives the joint likelihood for the fixed and random parameters, given the data, \code{grad_salam} gives the gradient of \code{nlogH_salam} with respect to the random effects and \code{hess_salam} gives its Hessian matrix. Finally \code{iLap.nlik_salam} computes the approximate marginal likelihood for the fixed parameters given the data.
##'
##' @aliases nlogH_salam(randPar, fixPar, data)
##' @aliases grad_salam(randPar, fixPar, data)
##' @aliases hess_salam(randPar, fixPar, data)
##' @aliases iLap.nlik_salam(fixPar, data, control)
##'
##' @param randPar The vector of random effects. The first 20 elements are the random effects for femals and the second block of 20 elements are the random effects for males
##' @param fixPar is the vector of fixed parameters made by the 4 regression coefficients and the 2 variance components
##' @param data The data. It must be a named list with appropriate components. See Examples for more details.
##' @param control A named list of parameters. See \code{\link[iLaplace]{iLap_an}} for more details.
##' @return double
##' @rdname logitRegCRE
##' @export
nlogH_salam <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplaceExtra_nlogH', PACKAGE = 'iLaplaceExtra', randPar, fixPar, data)
}
##' @return a 40-dimensional vector
##' @rdname logitRegCRE
##' @export
grad_salam <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplaceExtra_grad_nlogH', PACKAGE = 'iLaplaceExtra', randPar, fixPar, data)
}

##' @details For each experiment, the model can be written as
##' \deqn{Y_{ij} \sim Bernoulli(p_{ij})}{Y_ij \sim Bernoulli(p_ij),} \deqn{logit(p_{ij}) = x_{ij}^T\beta + u_i^f + u_j^m}{logit(p_i) = x_ij^T\beta + u_i^f + u_j^m} \deqn{u_i^f \sim N(0,\sigma_f^2)}{u_i \sim N(0,\sigma_m^2),} \eqn{i,j=1,\ldots,20}{i,j=1,...,20}. There are 20 random effects for female salamnders and 20 random effects for males. Both males' and females' random effects are independent Gaussian random variables but with different variances. There are 4 regression parameters \eqn{\beta}{\beta}, all corresponding to a type of cross. Hence, overall there are 6 fixed parameters of which 4 are regression parameter and 2 are variance parameters, reparametrised in logarithmic scale. If the experiments are analysed separately the marginal likelihood for the fixed parameters entails the integration over a 40-dimensional vector of random effects. But if the experiments are analysed jointly, there are \eqn{3\times 40=120}{3*40=120} random effects to be integrated out. The function \code{iLap.nlik_salam} integrates out the 40 random effects by the improved Laplace method of Ruli et al. (2015) and gives the negative of the logarithm of the integral.
##'
##' @references
##' Ruli E., Sartori N. and Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##'
##' @examples
##' \dontrun{
##' data(Salamander)
##' library(iLaplace)
##' # data for experiment 1 (salam1), 2 (salam2) and 3 (salam3)
##' salam1 = list()
##' salam1$y = salam$y[,1]
##' salam1$x = matrix(0, 120, 4)
##' salam1$x[,1] = rep(1, 120)
##' salam1$x[,2] = c(salam$x[,3]+salam$x[,4])
##' salam1$x[,3] = c(salam$x[,2]+salam$x[,4])
##' salam1$x[,4] = salam1$x[,2]*salam1$x[,3]
##' salam1$z = salam$z
##' salam2 = salam1
##' salam3 = salam1
##' salam2$y = salam$y[,2]
##' salam3$y = salam$y[,3]
##'
##' # maximise the marginal likelihood approximated by the improved
##' # Laplace approximation
##' init <- rep(0, 6)
##' opt1 <- nlminb(init,
##'                obj=function(x) -iLap.nlik_salam(x, data=salam1),
##'                control = list(trace=1)
##'                )
##' opt1
##' }
##' @return a \eqn{40\times 40}{40 by 40} positive definite matrix
##' @rdname logitRegCRE
##' @export
hess_salam <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplaceExtra_hess_nlogH', PACKAGE = 'iLaplaceExtra', randPar, fixPar, data)
}

# ldetHess_nlogH <- function(randPar, fixPar, data) {
#   if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
#   if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
#   if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
#     .Call('iLaplace_ldetHess_nlogH', PACKAGE = 'iLaplace', randPar, fixPar, data)
# }
##' @return a double, the minus logarithm of the integral
##' @rdname logitRegCRE
##' @export
iLap.nlik_salam <- function(fixPar, data,
                            control = list(sp.points = 100,
                                           delta = 20,
                                           nc.cores = detectCores()-1
                                           )
                            )
  {
  fixPar <-  c(fixPar[1:4], exp(fixPar[5:6]))
  ff <- function(x) nlogH_salam(randPar = x, fixPar = fixPar, data = data)
  ff.gr <- function(x) grad_salam(randPar = x, fixPar = fixPar, data = data)
  ff.hess <- function(x) hess_salam(randPar = x, fixPar = fixPar, data = data)

  fullOpt <- nlminb(start = rep(0, 40),
                    objective = ff,
                    gradient = ff.gr,
                    hessian = ff.hess
                    )

  fullOpt$hessian <- ff.hess(fullOpt$par)

  ans <- iLap_an2(fullOpt = fullOpt,
                 ff = ff,
                 ff.gr = ff.gr,
                 ff.hess = ff.hess,
                 control = control
                 )

  return(-ans)
}
