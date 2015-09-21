##' @name MHmcmc
##' @importFrom coda mcmc
##' @title Metropolis-Hastings sampling from User-Written R function
##' @description This function allows user to contruct a sample from a user-defined continous distribution using a random walk Metropolis algorithm.
##'
##' @usage MHmcmc(logfun, burnin, mcmc, thin, tune, V, df, theta.init, verbose)
##'
##' @param logfun The unnormalised log-density of the distribution from which to take a sample. This must be a function defined in R whose first argument is a continous (possibly vector) variable. This first argument is the point in the state space at which the (log)density is to be evaluated. See the Details section and the examples below for more information.
##' @param burnin The number of burn-in iterations for the sampler
##' @param mcmc The number of MCMC iterations after burnin
##' @param thin The thinning interval used in the simulation. The number of MCMC iterations must be divisible by this value
##' @param tune The tuning paramete.r for the Metropolis sampling. Can be either a positive scalar or a k-vector, where k is the length of theta
##' @param V The scale matrix for the Student's \eqn{t} proposal distribution. Must be a square matrix with dimension equal to the length of \code{theta.init}
##' @param df The degrees of freedom of the Student's \eqn{t} proposal distribution
##' @param theta.init Starting values for the sampling. Must be of the appropriate dimension. It must also be the case that \code{logfun(theta.init)} gives finite values, e.g. \code{-Inf} are not allowed.
##' @param verbose A switch which determines whether or not the progress of the sampler is printed to the screen. If verbose is greater than 0 the iteration number, and the Metropolis acceptance rate are sent to the screen every \code{verbose}th iteration
##'
##' @details \code{MHmcmc} produces a sample from a user-defined distribution using a random walk Metropolis algorithm with multivariate multivariate Student's \eqn{t}{t} proposal distribution. See Gelman et al. (2003) and Robert & Casella (2004) for details of the random walk Metropolis algorithm.
##'
##' The proposal distribution is centered at the current value of \eqn{\theta}{theta} and has scale matrix \eqn{H}. \eqn{H} is calculated as: \eqn{H = TVT}{H = T*V*T}, where \eqn{T}{T} is a the diagonal positive definite matrix formed from the \code{tune}.
##'
##' @return An mcmc object that contains the posterior sample. This object can be summarized by functions provided by the \code{coda} package.
##'
##' @examples
##' \dontrun{
##'
##' data(Lubricant)
##' y <- Lubricant$viscos
##' x1 <- Lubricant$tempC
##' x2 <- Lubricant$pressure/1000
##' n <- length(y)
##' muBeta <-  rep(0, 9)
##' SigBeta <-  matrix(0, 9, 9)
##' diag(SigBeta) <-  50
##' sigScale <-  10
##'
##' # re-define nlpost_lub, its gradient and Hessian matrix
##' ff <- function(x) nlpost_lub(beta = x[1:9], lsig = x[10],
##'                              y = y, x1 = x1, x2 = x2, n = n,
##'                              muBeta = muBeta, SigBeta = SigBeta,
##'                              sigScale = sigScale)
##' ff.gr <- function(x) grad_lub(beta = x[1:9], lsig = x[10], y = y,
##'                               x1 = x1, x2 = x2, n = n, muBeta = muBeta,
##'                               SigBeta = SigBeta, sigScale = sigScale)
##' ff.h <- function(x) hess_lub(beta = x[1:9], lsig = x[10], y = y, x1 = x1,
##'                              x2 = x2, n = n, muBeta = muBeta,
##'                              SigBeta = SigBeta, sigScale = sigScale)

##' # find minimum and the Hessian at the min value
##' init <-  c(10.522, 20.6, 1.464, -0.258, 0.022, 0.394, 0.036,
##'            46.566, -0.455, -3.195)

##' opt.post <-  nlminb(init, obj=ff, gradient = ff.gr,
##'                     hessian = ff.h, control = list(trace = 1))
##' opt.post$hessian = ff.h(opt.post$par)

##' # the mcmc sample
##' mcmc.post = MHmcmc(logfun = function(x) -ff(x),
##'                burnin = 40000, mcmc = 1e+5, thin = 1, tune = 0.5,
##'                V = solve(opt.post$hessian), df = dfprop,
##'                theta.init = init, verbose = 10000)

##' # the marginal posteriors
##' labl <- c("beta1", "beta2", "beta3", "beta4", "beta5",
##'           "beta6", "beta7", "beta8", "beta9", "log sgima")

##' panel.hist <- function(x, ...)
##' {
##' usr <- par("usr"); on.exit(par(usr))
##' par(usr = c(usr[1:2], 0, 1.5) )
##' h <- hist(x, plot = FALSE)
##' breaks <- h$breaks; nB <- length(breaks)
##' y <- h$counts; y <- y/max(y)
##' rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
##' }
##'
##' Lab.palette1 <- colorRampPalette(c("lightblue", "yellow",
##'                                     "orange", "red"),
##'                                      space = "Lab")
##'
##' pdf("pairPost.pdf", width=10, height=10)
##' par(mai=c(0.08,0.08,0.08,0.08))
##' pairs(as.matrix(mcmc.post[,1:10]),
##'       panel = function(...)
##'                 smoothScatter(...,
##'                               nrpoints = 50,
##'                               colramp = Lab.palette1,
##'                               add = TRUE
##'                               ),
##'       diag.panel = panel.hist,
##'       labels = labl)
##' dev.off()
##' }
##'
##' @references
##' Chib S. & Greenberg E. (1995).
##' Understanding the Metropolis-Hastings Algorithm. \emph{The American Statistician}, \bold{49}, 327-335.
##'
##' Gelman A., Carlin J. B., Stern H. S. & Rubin, D. B. (2003).
##' \emph{Bayesian Data Analysis}. 2nd Edition. Boca Raton: Chapman & Hall/CRC.
##'
##' Plummer M., Best, N., Cowles K. & Vines K. (2006).
##' CODA: Convergence Diagnosis and Output Analysis for MCMC. \emph{R News} \bold{6}, 7--1.
##'
##' Robert C. P. & Casella G. (2004).
##'  \emph{Monte Carlo Statistical Methods}. 2nd Edition. New York: Springer.
##'
##' @seealso \code{\link[iLaplaceExamples]{Lubricant}}, \code{\link[iLaplaceExamples]{nlpost_lub}},  \code{\link[iLaplaceExamples]{nlpost_gomp}}
##'
##' @export
MHmcmc <- function(logfun, burnin, mcmc, thin, tune, V, df, theta.init, verbose){
  my.env <- environment(fun = logfun)
  tune <- vector.tune(tune, length(theta.init))

  if (nrow(V) != ncol(V) || nrow(V) != length(theta.init)) {
    cat("V not of appropriate dimension.\n")
    stop("Check V and theta.init and call MHmcmc() again. \n",
         call. = FALSE)
  }
  CC <- NULL
  try(CC <- chol(V), silent = TRUE)
  if (is.null(CC)) {
    cat("V not positive definite.\n")
    stop("Check V and call MHmcmc() again. \n",
         call. = FALSE)
  }
  V <- tune %*% V %*% tune

  sample <- .Call('iLaplaceExamples_MCMCmetrop_cpp',
                  PACKAGE = 'iLaplaceExamples',
                  logfun,
                  theta.init,
                  V,
                  mcmc+burnin,
                  burnin,
                  df,
                  verbose,
                  my.env)
  sample <- mcmc(data = t(sample), start = burnin + 1, end = mcmc, thin = thin)
  return(sample)
}
