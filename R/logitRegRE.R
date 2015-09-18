##' @name logitRegRE
##' @title Logistic regression with random effects
##' @description The following functions give the posterior distribution and related quantities for a logistic regression model with random effects. Only the case with one observation for each group is considered, e.g. there is no replication within groups and the fixed effect is an intercept. In particular, \code{nlpost_rebin} gives the negative log-posterior for the fixed parameter (\code{theta}) and random effects (\code{u}), \code{grU_rebin} give the gradient of \code{nlpost_rebin} with respect to \code{u} and \code{hessU_rebin} gives the Hessian matrix.
##'
##' @aliases grU_rebin
##' @aliases hessU_rebin
##'
##' @param u The vector of random effects. There is one random effect for each observation.
##' @param theta The vector of fixed parameters made by a fixed effect parameter and the logartihm of the precision.
##' @param y The response.
##' @param n The sample size
##'
##' @return double
##' @export
nlpost_rebin <- function(u, theta, y, n) {
  .Call('iLaplaceExamples_nlpost_rebin', PACKAGE = 'iLaplaceExamples', u, theta, y, n)
}
##' @return an n+2 dimensional vector
##' @rdname logitRegRE
##' @export
grU_rebin <- function(u, theta, y, n) {
  .Call('iLaplaceExamples_grU_rebin', PACKAGE = 'iLaplaceExamples', u, theta, y, n)
}

# grAll_rebin <- function(u, theta, y, n) {
#   .Call('iLaplace_grAll_rebin', PACKAGE = 'iLaplace', u, theta, y, n)
# }
##' @examples
##' \dontrun{
##' #generate the  data
##' n <- 100
##' s2 <- 1
##' b <- 2
##' # model: logit(p_i) = b + u_i; y_i ~ dbinom(1, p_i)
##' set.seed(2)
##' u <- rnorm(n, 0, sqrt(s2))
##' p  <- 1/(1+exp(-(b + u)))
##' y <- rbinom(n,1, prob=p)
##'
##' #MCMCM with JAGS (needs JAGS)
##' library(rjags)
##' dataList <- list(y = y, N=n)
##'
##' # Define the model:
##' modelString = "
##' model {
##' for ( i in 1:N ) {
##' u[i] ~ dnorm(0.0, tau)
##' logit(p[i]) <- beta + u[i]
##' y[i] ~ dbern(p[i])
##' }
##' beta ~ dnorm(0.0, 1.0)
##' tau ~ dgamma(1.0,1.0)
##' ltau <- log(tau)
##' }
##' "
##'
##' writeLines(modelString, con="JAGSmodel.txt" )
##' initsList <- list(tau=1, beta=2, u = u)
##'
##' # Run the chains:
##' jagsModel = jags.model(file="JAGSmodel.txt",
##'                        data=dataList, inits=initsList,
##'                        n.chains=1, n.adapt=5000)
##' update(jagsModel, n.iter=1e+5)
##' codaSamples = coda.samples(jagsModel,
##'                            variable.names=c("beta", "ltau"),
##'                            n.iter = 1e+5)
##'
##' plot(coda.samples)
##' }
##'
##' @return an \eqn{(n+2)\times(n+2)}{(n+2)by(n+2)} positive definite matrix
##' @details The model can be written as \deqn{Y_i \sim Bernoulli(p_i)}{Y_i \sim Bernoulli(p_i),} \deqn{logit(p_i) = \beta + u_i}{logit(p_i) = \beta + u_i,} \deqn{u_i \sim N(0,\tau^-1)}{u_i \sim N(0,\tau^-1),} \deqn{\beta \sim N(0,1)}{\beta \sim N(0,1),} \deqn{\tau \sim Gamma(1,1)}{\tau \sim Gamma(1,1).} Note the \eqn{\tau}{\tau} is implemented in log.
##' @rdname logitRegRE
##' @export
hessU_rebin <- function(u, theta, y, n) {
  .Call('iLaplaceExamples_hessU_rebin', PACKAGE = 'iLaplaceExamples', u, theta, y, n)
}

# hessAll_rebin <- function(u, theta, y, n) {
#   .Call('iLaplace_hessAll_rebin', PACKAGE = 'iLaplace', u, theta, y, n)
# }

# grUB_rebin <- function(u, theta, y, n) {
#   .Call('iLaplace_grUB_rebin', PACKAGE = 'iLaplace', u, theta, y, n)
# }

# hessUB_rebin <- function(u, theta, y, n) {
#   .Call('iLaplace_hessUB_rebin', PACKAGE = 'iLaplace', u, theta, y, n)
# }
