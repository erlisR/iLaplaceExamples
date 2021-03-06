% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/nLinearReg.R
\name{nLinearReg}
\alias{gradT_bod2}
\alias{gradT_bod2,}
\alias{gradT_lub}
\alias{gradT_lub,}
\alias{grad_bod2}
\alias{grad_bod2,}
\alias{grad_lub}
\alias{grad_lub,}
\alias{hessT_bod2}
\alias{hessT_lub}
\alias{hessT_lub,}
\alias{hess_bod2}
\alias{hess_bod2,}
\alias{hess_lub}
\alias{hess_lub,}
\alias{nLinearReg}
\alias{nlpostT_bod2}
\alias{nlpostT_bod2,}
\alias{nlpostT_lub}
\alias{nlpostT_lub,}
\alias{nlpost_bod2}
\alias{nlpost_lub}
\title{Posterior distribution (and related quantities) for two nonlinear regression models}
\usage{
nlpost_lub(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)

nlpostT_lub(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)

grad_lub(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)

gradT_lub(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)

hess_lub(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)

hessT_lub(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)

nlpost_bod2(beta, lsig, y, Time, n, muBeta, SigBeta, sigScale)

nlpostT_bod2(beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale)

grad_bod2(beta, lsig, y, Time, n, muBeta, SigBeta, sigScale)

gradT_bod2(beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale)

hess_bod2(beta, lsig, y, Time, n, muBeta, SigBeta, sigScale)

hessT_bod2(beta, lsig, lnu, y, Time, n, muBeta, SigBeta, sigScale)
}
\arguments{
\item{beta}{The vector of regression parameters. For the Lubricant dataset beta is a 9-variate vector, for the BOD2 data beta is a bivariate vector}

\item{lsig}{The logarithm of the scale parameter}

\item{y}{The response vector}

\item{x1}{For Lubricant data only: the temperature of the lubricant}

\item{x2}{For Lubricant data only: the pressure of the lubricant}

\item{n}{The sample size}

\item{muBeta}{The center of the prior for beta. See Details}

\item{SigBeta}{The scale matrix of the prior for beta}

\item{sigScale}{For the models with Student's \eqn{t}-errors only: the value of the scale of the prior for the scale parameter.}

\item{lnu}{For the models with Student's \eqn{t}-errors only: the logarithm of the degrees of freedom}

\item{Time}{For BOD2 data only: the time}
}
\value{
double

double

10-variate vector

11-variate vector

\eqn{10\times 10}{10 by 10} matrix

\eqn{11\times 11}{11 by 11} matrix

double

double

3-variate vector

4-variate vector

\eqn{3\times 3}{3 by 3} matrix

\eqn{4\times 4}{4 by 4} matrix
}
\description{
Posterior distribution, and related quantities, for two non linear regression models each based on a different dataset. For each dataset, two error distributions are considered: normal and Student's \eqn{t}. The functions ending in "lub" refer to the \code{\link[iLaplaceExamples]{Lubricant}} dataset, those with ending "bod2" refer to the \code{\link[iLaplaceExamples]{BOD2}} dataset. The functions with an inner \code{T} refer to the model based on Student's \eqn{t}-error distribution. The functions starting with \code{nlpost} give the negative log-posterior density, those starting with \code{grad} give the gradient of \code{nlpost} and the functions starting with \code{hess} provide the Hessian matrix.
}
\details{
The prior for \code{beta} is the multivariate Studnet's \eqn{t}{t} with the appropriate dimension, with location \code{muBeta}, scale matrix \code{SigBeta} and degrees of freedom equal to 5. The prior for the \code{scale} is the HalfCauchy distribution with scale \code{sigScale} and the prior for the degress of freedom is the Jeffreys' rule prior of Fonseca et al. (2008).
}
\examples{
\dontrun{

library(iLaplace)

## regression models with BOD2 data
# set the data
data(BOD2)
y = BOD2$demand
Time = BOD2$Time
n = length(y)
muBeta = rep(0, 2)
SigBeta = matrix(0, 2,2)
diag(SigBeta) = 10
sigScale = 10
dfprop = 3

## Model with normal errors

# the objective function and derivatives
ff = function(pp, ...) nlpost_bod2(beta = pp[1:2], lsig = pp[3], ...)
ff.gr = function(pp, ...) grad_bod2(beta = pp[1:2], lsig = pp[3], ...)
ff.hess = function(pp, ...) hess_bod2(beta = pp[1:2], lsig = pp[3], ...)

# find the mode
init <- c(5, 2, 0.5)
opt.post = nlminb(init, obj = ff, gradient = ff.gr, hessian = ff.hess,
                  y = y, Time = Time, n =  n, muBeta = muBeta,
                  SigBeta = SigBeta, sigScale = sigScale, control = list(trace = 1))
opt.post$hessian = ff.hess(opt.post$par, y = y, Time = Time, n =  n, muBeta = muBeta,
                           SigBeta = SigBeta, sigScale = sigScale)

# simulate from the posterior
mcmc.post = MHmcmc(logfun = function(x) -ff(x, y = y, Time = Time, n =  n, muBeta = muBeta,
                                            SigBeta = SigBeta, sigScale = sigScale),
                   burnin = 20000, mcmc = 1e+6,
                   thin = 5, tune = 0.98, V = solve(opt.post$hessian),
                   df = dfprop, theta.init = init, verbose = 50000)


# marginal likelihood by Chib's method
logM.Chib <- ChibML(logfun = function(x) -ff(x, y = y, Time = Time, n =  n,
                                            muBeta = muBeta, SigBeta = SigBeta,
                                            sigScale = sigScale),
                   theta.star = opt.post$par,
                   tune = 1.5, V = solve(opt.post$hessian), t(mcmc.post),
                   df = dfprop, verbose = 100000)

# marginal likelihood by the improved Laplace
logM.iLap <- iLaplace(fullOpt = opt.post, ff = ff, ff.gr = ff.gr, ff.hess = ff.hess,
                      control = list(sp.points = 100, delta = 13, n.cores = 3),
                      clEvalQ = c("iLaplaceExamples"), y = y, Time = Time, n =  n,
                      muBeta = muBeta, SigBeta = SigBeta,
                      sigScale = sigScale)

logM.Lap <- 3*log(2*pi)/2 - opt.post$objective - 0.5*determinant(opt.post$hessian)$mod
# marginal likelihood by importance samling
logM.IS <- isML(logfun = function(x) -ff(x, y = y, Time = Time, n =  n,
                                         muBeta = muBeta, SigBeta = SigBeta,
                                         sigScale = sigScale), nsim = 1e+6,
                theta.hat = opt.post$par, tune = 1.5,
                V = solve(opt.post$hessian),df = dfprop, verbose = 1000)
logM.iLap
logM.Lap
logM.IS
logM.Chib

## Model with Student's t errors

# the objective function and derivatives
tff = function(pp, ...) nlpostT_bod2(beta = pp[1:2], lsig = pp[3], lnu = pp[4], ...)
tff.gr = function(pp, ...) gradT_bod2(beta = pp[1:2], lsig = pp[3], lnu = pp[4], ...)
tff.hess = function(pp, ...) hessT_bod2(beta = pp[1:2], lsig = pp[3], lnu = pp[4], ...)

# find the mode
init <- c(2, 5, 0.5, 0)
opt.postt = nlminb(init, obj = tff, gradient = tff.gr, hessian = tff.hess,
                  y = y, Time = Time, n =  n, muBeta = muBeta,
                  SigBeta = SigBeta, sigScale = sigScale, control = list(trace = 1))
opt.postt$hessian = tff.hess(opt.postt$par, y = y, Time = Time, n =  n, muBeta = muBeta,
                           SigBeta = SigBeta, sigScale = sigScale)

# simulate from the posterior
mcmc.postt = MHmcmc(logfun = function(x) -tff(x, y = y, Time = Time, n =  n, muBeta = muBeta,
                                            SigBeta = SigBeta, sigScale = sigScale),
                   burnin = 20000, mcmc = 1e+6,
                   thin = 5, tune = 1.2, V = solve(opt.postt$hessian),
                   df = dfprop, theta.init = opt.postt$par, verbose = 50000)

logM.Lapt <- 4*log(2*pi)/2 - opt.postt$objective - 0.5*determinant(opt.postt$hessian)$mod

# marginal likelihood by Chib's method
logM.Chibt <- ChibML(logfun = function(x) -tff(x, y = y, Time = Time, n =  n,
                                            muBeta = muBeta, SigBeta = SigBeta,
                                            sigScale = sigScale),
                   theta.star = opt.postt$par,
                   tune = 1.5, V = solve(opt.postt$hessian), t(mcmc.postt),
                   df = dfprop, verbose = 100000)

logM.iLapt <- iLaplace(fullOpt = opt.postt, ff = tff, ff.gr = tff.gr, ff.hess = tff.hess,
                      control = list(sp.points = 100, delta = 9, n.cores = 4),
                      clEvalQ = c("iLaplaceExamples"), y = y, Time = Time, n =  n,
                      muBeta = muBeta, SigBeta = SigBeta,
                      sigScale = sigScale)

logM.ISt <- isML(logfun = function(x) -tff(x, y = y, Time = Time, n =  n,
                                         muBeta = muBeta, SigBeta = SigBeta,
                                         sigScale = sigScale), nsim = 1e+6,
                theta.hat = opt.postt$par, tune = 1.5,
                V = solve(opt.postt$hessian),df = dfprop, verbose = 10000)

logM.iLapt
logM.Lapt
logM.ISt
logM.Chibt
}
}
\references{
Fonseca T. C., Ferreira M. A. R. & Migon H. S. (2008)
Objective Bayesian analysis for the Student-\eqn{t}{t} regression model.
\emph{Biometrika} \bold{95}, 325--333.

Ruli E., Sartori N. & Ventura L. (2015)
Improved Laplace approximation for marignal likelihoods.
\url{http://arxiv.org/abs/1502.06440}
}
\seealso{
\code{\link[iLaplaceExamples]{MHmcmc}}
}

