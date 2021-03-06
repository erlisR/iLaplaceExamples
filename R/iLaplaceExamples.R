#' iLaplaceExamples: A package for examples with iLaplace.
#'
#' The iLaplaceExamples package provides utilites for running examples of approximating marignal likelihoods by the improved Laplace method of Ruli et al. (2015), implemented in the package iLaplace. Other approximation method such as the standard Laplace and some Monte Carlo methods are aslo provided. There are five examples:
#' \enumerate{
#'  \item Gompertz model
#'  \item multivariate t/skew-t distribtuion
#'  \item nonlinear regression models with two the \code{\link[iLaplaceExamples]{BOD2}} and the \code{\link[iLaplaceExamples]{Lubricant}} dataset
#'  \item binary logistic regression with random effects
#'  \item binary logistic regression with crossed random effects with the \code{\link[iLaplaceExamples]{Salamander}} data}
#' which functions are listed below.
#' foo, bar and baz.
#'
#' @section Example 1 functions:
#' \code{\link[iLaplaceExamples]{nlpost_gomp}}, \code{\link[iLaplaceExamples]{grad_gomp}}, \code{\link[iLaplaceExamples]{hess_gomp}}
#'
#' @section Example 2 functions:
#' \code{\link[iLaplaceExamples]{iLap_mvt}}, \code{\link[iLaplaceExamples]{iLap_mvskt}}
#' 
#' @section Example 3 functions:
#' \code{\link[iLaplaceExamples]{nlpost_lub}}, \code{\link[iLaplaceExamples]{nlpostT_lub}}, \code{\link[iLaplaceExamples]{grad_lub}}, \code{\link[iLaplaceExamples]{gradT_lub}}, \code{\link[iLaplaceExamples]{hess_lub}}, \code{\link[iLaplaceExamples]{hessT_lub}}, \code{\link[iLaplaceExamples]{nlpost_bod2}}, \code{\link[iLaplaceExamples]{nlpostT_bod2}}, \code{\link[iLaplaceExamples]{grad_bod2}}, \code{\link[iLaplaceExamples]{gradT_bod2}}, \code{\link[iLaplaceExamples]{hess_bod2}}, \code{\link[iLaplaceExamples]{hessT_bod2}}
#' 
#' @section Example 4 functions:
#' \code{\link[iLaplaceExamples]{nlpost_rebin}}, \code{\link[iLaplaceExamples]{grU_rebin}}, \code{\link[iLaplaceExamples]{hessU_rebin}}
#' 
#' @section Example 5 functions:
#' \code{\link[iLaplaceExamples]{nlogH_salam}}, \code{\link[iLaplaceExamples]{grad_salam}}, \code{\link[iLaplaceExamples]{hess_salam}}
#' 
#' @docType package
#' @name iLaplaceExamples
NULL
