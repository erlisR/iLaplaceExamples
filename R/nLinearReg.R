##' @name nLinearReg
##'
##' @title Posterior distribution (and related quantities) for two nonlinear regression models
##' @description Posterior distribution, and related quantities, for two non linear regression models each based on a different dataset. For each dataset, two error distributions are considered: normal and Student's \eqn{t}. The functions ending in "lub" refer to the \code{\link[iLaplaceExamples]{Lubricant}} dataset, those with ending "bod2" refer to the \code{\link[iLaplaceExamples]{BOD2}} dataset. The functions with an inner \code{T} refer to the model based on Student's \eqn{t}-error distribution. The functions starting with \code{nlpost} give the negative log-posterior density, those starting with \code{grad} give the gradient of \code{nlpost} and the functions starting with \code{hess} provide the Hessian matrix.
##'
##' @aliases nlpostT_lub, grad_lub, gradT_lub, hess_lub, hessT_lub, nlpostT_bod2, grad_bod2, gradT_bod2, hess_bod2, hessT_bod2
##'
##' @param beta The vector of regression parameters. For the Lubricant dataset beta is a 9-variate vector, for the BOD2 data beta is a bivariate vector
##' @param lsig The logarithm of the scale parameter
##' @param lnu For the models with Student's \eqn{t}-errors only: the logarithm of the degrees of freedom
##' @param y The response vector
##' @param x1 For Lubricant data only: the temperature of the lubricant
##' @param x2 For Lubricant data only: the pressure of the lubricant
##' @param Time For BOD2 data only: the time
##' @param n The sample size
##' @param muBeta The center of the prior for beta. See Details
##' @param SigBeta The scale matrix of the prior for beta
##' @param sigScale For the models with Student's \eqn{t}-errors only: the value of the scale of the prior for the scale parameter.
##'
##'
##' @return double
##' @export
nlpost_lub <- function(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_nlpost_lub', PACKAGE = 'iLaplaceExamples', beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)
}
##' @return double
##' @rdname nLinearReg
##' @export
nlpostT_lub <- function(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_nlpostT_lub', PACKAGE = 'iLaplaceExamples', beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

##' @return 10-variate vector
##' @rdname nLinearReg
##' @export
grad_lub <- function(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_grad_lub', PACKAGE = 'iLaplaceExamples', beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

##' @return 11-variate vector
##' @rdname nLinearReg
##' @export
gradT_lub <- function(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_gradT_lub', PACKAGE = 'iLaplaceExamples', beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

##' @return \eqn{10\times 10}{10 by 10} matrix
##' @rdname nLinearReg
##' @export
hess_lub <- function(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_hess_lub', PACKAGE = 'iLaplaceExamples', beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)
}
##' @return \eqn{11\times 11}{11 by 11} matrix
##' @rdname nLinearReg
##' @export
hessT_lub <- function(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_hessT_lub', PACKAGE = 'iLaplaceExamples', beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

##' @return double
##' @rdname nLinearReg
##' @export
nlpost_bod2 <- function(beta, lsig, y, Time, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_nlpost_bod2', PACKAGE = 'iLaplaceExamples', beta, lsig, y, Time, n, muBeta, SigBeta, sigScale)
}

##' @return double
##' @rdname nLinearReg
##' @export
nlpostT_bod2 <- function(beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_nlpostT_bod2', PACKAGE = 'iLaplaceExamples', beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale)
}

##' @return 3-variate vector
##' @rdname nLinearReg
##' @export
grad_bod2 <- function(beta, lsig, y, Time, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_grad_bod2', PACKAGE = 'iLaplaceExamples', beta, lsig, y, Time, n, muBeta, SigBeta, sigScale)
}

##' @return 4-variate vector
##' @rdname nLinearReg
##' @export
gradT_bod2 <- function(beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale) {
  .Call('iLaplaceExamples_gradT_bod2', PACKAGE = 'iLaplaceExamples', beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale)
}

##' @return \eqn{3\times 3}{3 by 3} matrix
##' @rdname nLinearReg
##' @export
hess_bod2 <- function(beta, lsig, y, Time, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_hess_bod2', PACKAGE = 'iLaplaceExamples', beta, lsig, y, Time, n, muBeta, SigBeta, sigScale)
}
##' @examples
##' \dontrun{
##' 
##' library(iLaplace)
##' 
##' ## regression models with BOD2 data
##' # set the data
##' data(BOD2)
##' y = BOD2$demand
##' Time = BOD2$Time
##' n = length(y)
##' muBeta = rep(0, 2)
##' SigBeta = matrix(0, 2,2)
##' diag(SigBeta) = 10
##' sigScale = 10
##' dfprop = 3
##' 
##' # the objective function and derivatives
##' ff = function(pp, ...) nlpost_bod2(beta = pp[1:2], lsig = pp[3], ...)
##' ff.gr = function(pp, ...) grad_bod2(beta = pp[1:2], lsig = pp[3], ...)
##' ff.hess = function(pp, ...) hess_bod2(beta = pp[1:2], lsig = pp[3], ...)
##' 
##' # find the mode
##' init <- c(5, 2, 0.5)
##' opt.post = nlminb(init, obj = ff, gradient = ff.gr, hessian = ff.hess,
##'                   y = y, Time = Time, n =  n, muBeta = muBeta,
##'                   SigBeta = SigBeta, sigScale = sigScale, control = list(trace = 1))
##' opt.post$hessian = ff.hess(opt.post$par, y = y, Time = Time, n =  n, muBeta = muBeta,
##'                            SigBeta = SigBeta, sigScale = sigScale)
##' 
##' # simulate from the posterior
##' mcmc.post = MHmcmc(logfun = function(x) -ff(x, y = y, Time = Time, n =  n, muBeta = muBeta,
##'                                             SigBeta = SigBeta, sigScale = sigScale),
##'                    burnin = 20000, mcmc = 1e+6,
##'                    thin = 5, tune = 0.98, V = solve(opt.post$hessian),
##'                    df = dfprop, theta.init = init, verbose = 50000)
##' 
##' 
##' # marginal likelihood by Chib's method
##' logM.Chib <- ChibML(logfun = function(x) -ff(x, y = y, Time = Time, n =  n,
##'                                             muBeta = muBeta, SigBeta = SigBeta,
##'                                             sigScale = sigScale),
##'                    theta.star = opt.post$par,
##'                    tune = 1.5, V = solve(opt.post$hessian), t(mcmc.post),
##'                    df = dfprop, verbose = 100000)
##' 
##' # marginal likelihood by the improved Laplace
##' logM.iLap <- iLaplace(fullOpt = opt.post, ff = ff, ff.gr = ff.gr, ff.hess = ff.hess,
##'                       control = list(sp.points = 100, delta = 13, n.cores = 3),
##'                       clEvalQ = c("iLaplaceExamples"), y = y, Time = Time, n =  n,
##'                       muBeta = muBeta, SigBeta = SigBeta,
##'                       sigScale = sigScale)
##' 
##' logM.Lap <- 3*log(2*pi)/2 - opt.post$objective - 0.5*determinant(opt.post$par)$mod
##' # marginal likelihood by importance samling
##' logM.IS <- isML(logfun = function(x) -ff(x, y = y, Time = Time, n =  n,
##'                                          muBeta = muBeta, SigBeta = SigBeta,
##'                                          sigScale = sigScale), nsim = 1e+6,
##'                 theta.hat = opt.post$par, tune = 1.5,
##'                 V = solve(opt.post$hessian),df = dfprop, verbose = 1000)
##' logM.iLap
##' logM.Lap
##' logM.IS
##' logM.Chib
##' 
##' }
##' 
##' 
##' @seealso \code{\link[iLaplaceExamples]{MHmcmc}}
##'
##'  @details The prior for \code{beta} is the multivariate Studnet's \eqn{t}{t} with the appropriate dimension, with location \code{muBeta}, scale matrix \code{SigBeta} and degrees of freedom equal to 5. The prior for the \code{scale} is the HalfCauchy distribution with scale \code{sigScale} and the prior for the degress of freedom is the Jeffreys' rule prior of Fonseca et al. (2008).

##' @references
##' Fonseca T. C., Ferreira M. A. R. & Migon H. S. (2008)
##' Objective Bayesian analysis for the Student-\eqn{t}{t} regression model.
##' \emph{Biometrika} \bold{95}, 325--333.
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}

##' @return \eqn{4\times 4}{4 by 4} matrix
##' @rdname nLinearReg
##' @export
hessT_bod2 <- function(beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplaceExamples_hessT_bod2', PACKAGE = 'iLaplaceExamples', beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale)
}

