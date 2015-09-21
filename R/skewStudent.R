nlDen_mvt <- function(x, m, df) {
    .Call('iLaplaceExamples_nlDen_mvt', PACKAGE = 'iLaplaceExamples', x, m, df)
}

nlDen_mvskt <- function(x, a, c, m, df) {
    .Call('iLaplaceExamples_nlDen_mvskt', PACKAGE = 'iLaplaceExamples', x, a, c, m, df)
}

grad_mvt <- function(x, m, df) {
    .Call('iLaplaceExamples_grad_mvt', PACKAGE = 'iLaplaceExamples', x, m, df)
}

grad_mvskt <- function(x, a, c, m, df) {
    .Call('iLaplaceExamples_grad_mvskt', PACKAGE = 'iLaplaceExamples', x, a, c, m, df)
}

hess_mvt <- function(x, m, df) {
    .Call('iLaplaceExamples_hess_mvt', PACKAGE = 'iLaplaceExamples', x, m, df)
}

hess_mvskt <- function(x, a, c, m, df) {
    .Call('iLaplaceExamples_hess_mvskt', PACKAGE = 'iLaplaceExamples', x, a, c, m, df)
}

##' @name skewStudent
##' @title Improved and standard Laplace approximations for multivariate skewed Student's \eqn{t}{t} densities
##' @description These functions give the normalising constant of (skewed) multivariate Student's t densities approximated by the standard and the improved Laplace of Ruli et al. (2015). In particular, \code{iLap_mvt} provides the normalising constant of the \code{m}-variate Student's \eqn{t}{t} density with \code{df} degrees of freedom, centered at zero and indentiy scale matrix. \code{iLap_mvskt} provides the normalising constanto of the \code{m}-variate t/skew-t density of Jones (2002).
##'
##' @aliases iLap_mvt(m, df)
##' @aliases iLap_mvskt(m, df)
##'
##' @param m The dimensionality
##' @param df The degrees of freedom
##' @param a The parameter a of the t/skew-t density of Jones(2002)
##' @param c The parameter c of the t/skew-t density of Jones(2002)
##' @return a list with components \code{Lap} and \code{iLap}
##' @details These functions use the improved and the standard Laplace approximation for computing the normalising constant of the multivariate t/skew-t density. The functions are similar to \code{\link[iLaplace]{iLaplace}}.
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##'
##' Jones M. C. (2002).
##' Marginal replacement in multivariate densities, with application to skewing spherically symmetric distributions.
##' \emph{Journal of Multivariate Analysis} \bold{81}, 85--99.
##'
##' @seealso \code{\link[iLaplace]{iLaplace}}
##'
##' @examples 
##' \dontrun{
##' 
##' library(iLaplaceExamples)
##' 
##' # improved (black) and standard(red) Laplace approximations for the multivariate t
##' dims = seq(3,100, by = 5)
##' lap3 = sapply(dims, function(i) iLap_mvt(i, 3))
##' lap5 = sapply(dims, function(i) iLap_mvt(i, 5))
##' lap10 = sapply(dims, function(i) iLap_mvt(i, 10))
##' lap20 = sapply(dims, function(i) iLap_mvt(i, 20))
##' 
##' # improved (black) and standard (red) Laplace approximations for the multivariate t/skewt
##' # a=4, c=1
##' lap3s2 = sapply(dims, function(i) iLap_mvskt(a=4, c=1, i, df=3))
##' lap5s2 = sapply(dims, function(i) iLap_mvskt(a=4, c=1, i, df=5))
##' lap10s2 = sapply(dims, function(i) iLap_mvskt(a=4, c=1, i, df=10))
##' lap20s2 = sapply(dims, function(i) iLap_mvskt(a=4, c=1, i, df=20))
##' 
##' postscript("fig1.eps", width=5, height=5)
##' # red = standard Laplace
##' # black = improved Laplace
##' par(par(mai=c(0.85,0.85,0.2,0.03)))
##' plot(dims, lap3[1,], pch=1, cex=1.2, ylim=c(-55, 0), xlab=NA, main=NA, ylab=NA)
##' points(dims, lap5[1,], pch=2, cex=1.2, col=1)
##' points(dims, lap10[1,], pch=3, cex=1.2, col=1)
##' points(dims, lap20[1,], pch=4, cex=1.2, col=1)
##' points(dims, lap3[2,], pch=1, cex=1.2, col=2, lwd = 2)
##' points(dims, lap5[2,], pch=2, cex=1.2, col=2,lwd = 2)
##' points(dims, lap10[2,], pch=3, cex=1.2, col=2,lwd = 2)
##' points(dims, lap20[2,], pch=4, cex=1.2, col=2,lwd = 2)
##' mtext(text="dimension", side=1,line=2.2, cex = 1.3)
##' mtext(text="log normalising constant", side=2,line=2.5, cex = 1.3)
##' mtext(text="a=c=3/2", side=3,line=0.4, cex = 1.3)
##' # points(dims, lap3[2,], pch="°", cex=1.3)
##' abline(h=0, lwd = 2)
##' dev.off()
##' 
##' postscript("fig2.eps", width=5, height=5)
##' # red = standard Laplace
##' # black = improved Laplace
##' # par(par(mai=c(0.85,0.85,0.2,0.03)))
##' plot(dims, lap3s2[1,], pch=1, cex=1.2, ylim=c(-50, 0), xlab=NA, main=NA, ylab=NA)
##' points(dims, lap5s2[1,], pch=2, cex=1.2, col=1)
##' points(dims, lap10s2[1,], pch=3, cex=1.2, col=1)
##' points(dims, lap20s2[1,], pch=4, cex=1.2, col=1)
##' points(dims, lap3s2[2,], pch=1, cex=1.2, col=2, lwd = 2)
##' points(dims, lap5s2[2,], pch=2, cex=1.2, col=2,lwd = 2)
##' points(dims, lap10s2[2,], pch=3, cex=1.2, col=2,lwd = 2)
##' points(dims, lap20s2[2,], pch=4, cex=1.2, col=2,lwd = 2)
##' mtext(text="dimension", side=1,line=2.2)
##' mtext(text="log normalising contant", side=2.5,line=1.7)
##' mtext(text="a=4, c=1", side=3,line=0.4)
##' # points(dims, lap3[2,], pch="°", cex=1.3)
##' abline(h=0, lwd = 2)
##' dev.off()
##' }
##' @rdname skewStudent
##' @export
iLap_mvt = function(m, df) {
  # full optimization and hessian
  fullOpt = nlminb(rep(0,m), objective = nlDen_mvt, gradient = grad_mvt, m = m, df = df)
  fullOpt$hessian = hess_mvt(fullOpt$par, m, df)
  fullOpt$ldet = as.double(determinant(fullOpt$hessian)$mod)

  lo = fullOpt$par - 100*sqrt(diag(solve(fullOpt$hessian)))
  up = fullOpt$par + 100*sqrt(diag(solve(fullOpt$hessian)))
  lnc = lden = 0.0

  # approximate marginal denstiy for the first component
  marg = function(y) {
    tmp = function(xx) nlDen_mvt(c(y, xx), m=m, df=df)
    gr.tmp = function(xx) grad_mvt(c(y, xx), m=m, df=df)[-1]
    tmpOpt = nlminb(fullOpt$par[2:m], objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvt(c(y, tmpOpt$par), m=m, df=df)[-1,-1]
    ldetc = determinant(tmpOpt$hessian)$mod
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc) - tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }

  # normalizes the marginal density and computes the density at the mode
  nc.marg = integrate(Vectorize(function(x) marg(x), "x"), lower=lo[1], upper=up[1])$value
  lnc = log(nc.marg) + lnc
  lden = log(marg(fullOpt$par[1])) + lden

  # conditionals form the second to the p-1 th component
  middleConds <- function(y, index){
    # index = 2, ..., m -1.
    no = c(1:index)
    tmp = function(xx) nlDen_mvt(c(fullOpt$par[1:(index-1)], y, xx), m=m, df=df)
    gr.tmp = function(xx) grad_mvt(c(fullOpt$par[1:(index-1)], y, xx), m=m, df=df)[-no]
    tmpOpt = nlminb(c(fullOpt$par[-no]), objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvt(c(fullOpt$par[1:(index-1)], y, tmpOpt$par), m=m, df=df)[-no,-no]
    if(index==(m-1)) {
      ldetc = log(tmpOpt$hessian)
    } else {
      ldetc = determinant(tmpOpt$hessian)$mod
    }
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc)- tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }

  # normalizes the conditionals and evaluates their density at the mode
  for(i in 2:(m-1)){
    nc.midcond = integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lo[i], upper=up[i])$value
    lnc = lnc + log(nc.midcond)
    lden = log(middleConds(fullOpt$par[i], index=i)) + lden
  }

  # The conditional of pth given the rest
  lastCond = function(y) {
    tt = -nlDen_mvt(c(fullOpt$par[1:(m-1)], y), m=m, df=df) + fullOpt$objective
    return(exp(tt))
  }

  # Normalized the last conditional and evaluates the density at the mode
  nc.lastcond = integrate(Vectorize(function(x) lastCond(x), "x"), lower=lo[m], upper=up[m])$value

  lnc = log(nc.lastcond) + lnc
  lden = log(lastCond(fullOpt$par[m])) + lden

  ans = c(-fullOpt$obj - lden + lnc, 0.5*m*log(2*pi) - fullOpt$obj - 0.5*fullOpt$ldet)
  return(ans)
}

##' @return a list with components \code{Lap} and \code{iLap}
##' @rdname skewStudent
##' @export
iLap_mvskt = function(a, c, m, df) {
  # full optimization and hessian
  fullOpt = nlminb(rep(0,m), objective = nlDen_mvskt, gradient=grad_mvskt, a=a, c=c, m=m, df=df)
  fullOpt$hessian = hess_mvskt(fullOpt$par, a=a, c=c, m=m, df=df)
  fullOpt$ldet = as.double(determinant(fullOpt$hessian)$mod)

  lo = fullOpt$par - 100*sqrt(diag(solve(fullOpt$hessian)))
  up = fullOpt$par + 100*sqrt(diag(solve(fullOpt$hessian)))
  lnc = lden = 0.0

  # approximate marginal denstiy for the first component
  marg = function(y) {
    tmp = function(xx) nlDen_mvskt(c(y, xx), a=a, c=c, m=m, df=df)
    gr.tmp = function(xx) grad_mvskt(c(y, xx), a=a, c=c, m=m, df=df)[-1]
    tmpOpt = nlminb(fullOpt$par[2:m], objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvskt(c(y, tmpOpt$par), a=a, c=c, m=m, df=df)[-1,-1]
    ldetc = determinant(tmpOpt$hessian)$mod
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc) - tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }

  # normalizes the marginal density and computes the density at the mode
  nc.marg = integrate(Vectorize(function(x) marg(x), "x"), lower=lo[1], upper=up[1])$value
  lnc = log(nc.marg) + lnc
  lden = log(marg(fullOpt$par[1])) + lden

  # conditionals form the second to the p-1 th component
  middleConds <- function(y, index){
    # index = 2, ..., m -1.
    no = c(1:index)
    tmp = function(xx) nlDen_mvskt(c(fullOpt$par[1:(index-1)], y, xx), a=a, c=c, m=m, df=df)
    gr.tmp = function(xx) grad_mvskt(c(fullOpt$par[1:(index-1)], y, xx), a=a, c=c, m=m, df=df)[-no]
    tmpOpt = nlminb(c(fullOpt$par[-no]), objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvskt(c(fullOpt$par[1:(index-1)], y, tmpOpt$par), a=a, c=c, m=m, df=df)[-no,-no]
    if(index==(m-1)) {
      ldetc = log(tmpOpt$hessian)
    } else {
      ldetc = determinant(tmpOpt$hessian)$mod
    }
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc)- tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }

  # normalizes the conditionals and evaluates their density at the mode
  for(i in 2:(m-1)){
    nc.midcond = integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lo[i], upper=up[i])$value
    lnc = lnc + log(nc.midcond)
    lden = log(middleConds(fullOpt$par[i], index=i)) + lden
  }

  # The conditional of pth given the rest
  lastCond = function(y) {
    tt = -nlDen_mvskt(c(fullOpt$par[1:(m-1)], y), a=a, c=c, m=m, df=df) + fullOpt$objective
    return(exp(tt))
  }

  # Normalized the last conditional and evaluates the density at the mode
  nc.lastcond = integrate(Vectorize(function(x) lastCond(x), "x"), lower=lo[m], upper=up[m])$value

  lnc = log(nc.lastcond) + lnc
  lden = log(lastCond(fullOpt$par[m])) + lden

  ans = c(-fullOpt$obj - lden + lnc, 0.5*m*log(2*pi) - fullOpt$obj - 0.5*fullOpt$ldet)
  return(ans)
}

