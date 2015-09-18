##' @name gompertz
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
##' set.seed(12)
##'
##' # generate the data
##' data <- rgompertz(n = 50, shape = exp(2), scale = exp(3))
##'
##' #find the posterior mode and Hessian at the mode
##' opt.post = nlminb(c(2,3), obj = nlpost_gomp,
##'                   gradient=grad_gomp,data = data)
##' opt.post$hessian = hessian(nlpost_gomp,
##'                             opt.post$par, data = data)
##'
##'  # draw a posterior sample
##' mcmc.gomp = MHmcmc(logfun = function(x) -nlpost_gomp(x, data=data),
##'                     burnin=5000, mcmc=1e+5, thin=2, tune=1.1,
##'                     V=solve(opt.post$hessian), df=5,
##'                     theta.init=opt.post$par, verbose=10000)
##'
##' plot(mcmc.gomp)
##'}
##' @details The prior distribution is bivariate normal with independent components, with mean (0,0) and variance (100, 100).
##'
##' @seealso \code{\link[iLaplaceExtra]{MHmcmc}}
##' @rdname gompertz
##' @export
hess_gomp <- function(param, data) {
    .Call('iLaplaceExtra_hess_gomp', PACKAGE = 'iLaplaceExtra', param, data)
}
