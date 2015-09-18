##' @name gompertz
##' @import Rcpp
##' @useDynLib iLaplaceExamples
##' @title Posterior distribution (and related quantities) for the Gompertz model
##'
##' @description Function \code{nlpost_gomp} gives the non-normalised negative log-posterior distribution for the Gompertz distribution, \code{grad_gomp} gives its gradient and \code{hess_gomp} gives the Hessian matrix. The parameter of interest is (exp(\eqn{\alpha}), exp(\eqn{\beta})), with (\eqn{\alpha,\beta}) being reals. See Details for more information.
##'
##' @aliases grad_gomp
##' @aliases hess_gomp
##'
##' @param param The bivariate vector of parameters.
##' @param data The vector of data.
##'
##' @return double
##'
##' @export
nlpost_gomp <- function(param, data) {
  .Call('iLaplaceExamples_nlpost_gomp', PACKAGE = 'iLaplaceExamples', param, data)
}
##' @return bivariate vector
##' @rdname gompertz
##' @export
grad_gomp <- function(param, data) {
  .Call('iLaplaceExamples_grad_gomp', PACKAGE = 'iLaplaceExamples', param, data)
}

##' @return \eqn{2 \times 2}{2 by 2} matrix
##' @examples
##'\dontrun{
##'
##' library(VGAM)
##' library(iLaplace)
##' library(iLaplaceExamples)
##' set.seed(12)
##' 
##' # generate the data
##' data <- rgompertz(n = 50, shape = exp(2), scale = exp(3))
##' 
##' ff <- function(x, ...) nlpost_gomp(param = x, ...)
##' ff.gr <- function(x, ...) grad_gomp(param = x, ...)
##' ff.hess <- function(x, ...) hess_gomp(param = x, ...)
##' 
##' # find the posterior mode and Hessian at the mode
##' opt.post = nlminb(c(2,3), obj = ff, gradient = ff.gr,
##'                   hessian = ff.hess, data = data)
##' opt.post$hessian = ff.hess(opt.post$par, data = data)
##' 
##' # draw a posterior sample
##' mcmc.gomp = MHmcmc(logfun = function(x) -ff(x, data),
##'                    burnin = 5000, mcmc = 1e+6, thin = 2, tune = 1.1,
##'                    V = solve(opt.post$hessian), df = 5,
##'                    theta.init = opt.post$par, verbose = 10000)
##' # look at the plots
##' plot(mcmc.gomp)
##' 
##' # marginal likelihood by the improved Laplace of Ruli et al. (2015)
##' iLap.logM <- iLaplace(fullOpt = opt.post, ff = ff, ff.gr = ff.gr, ff.hess = ff.hess,
##'                       control = list(sp.points = 100, delta = 13, n.cores = 2),
##'                       clEvalQ = c("iLaplaceExamples"), data = data)
##' 
##' # marginal likelihood by the standard Laplace
##' Lap.logM <- log(2*pi) - opt.post$objective - 0.5*determinant(opt.post$hessian)$mod
##' 
##' # marginal likelihood by improtance sampling
##' is.logM <- isML(logfun = function(x) -ff(x, data), nsim = 1e+6,
##'                 theta.hat = opt.post$par, tune = 1.2, V = solve(opt.post$hessian),
##'                df = 5, verbose = 100000)
##'}
##' @details The prior distribution is bivariate normal with independent components, with mean (0,0) and variance (100, 100).
##'
##' @seealso \code{\link[iLaplaceExamples]{MHmcmc}}
##' @rdname gompertz
##' @export
hess_gomp <- function(param, data) {
  .Call('iLaplaceExamples_hess_gomp', PACKAGE = 'iLaplaceExamples', param, data)
}
