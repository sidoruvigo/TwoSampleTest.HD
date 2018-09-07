#' Package ‘TwoSampleTest.HD’
#'
#' Documentation for package ‘TwoSampleTest.HD’ version 1.0
#' @aliases TwoSampleTest.HDpackage
#' @title A two-sample test for the equality of distributions for high-dimensional data
#'
#' @description This package implements four different tests proposed in Cousido-Rocha et al. (2018). These methods test the (global) null hypotheses of equality of the univariate marginals of the p-variate distributions in the two populations.
#' In other words, the null hypotheses is an intersection of the p null hypotheses corresponding to p different two-sample problems.
#' These methods are particularly well suited to the low sample size, high dimensional setting (n << p).
#' The sample size can be as small as 2. The test accounts for the possibility that the p variables in each data set can be weakly dependent.
#' Three of the methods arise from different approaches to estimate the variance of the same statistic.
#' This statistic averages p individual statistics based on comparing the empirical characteristic functions computed from the two samples.
#' The last method is an alternative global test whose statistic averages the p-values derived from applying permutation tests to the individual statistics mentioned above.
#' When the global null hypotheses is rejected such permutation p-values can also be used to identify which variables contribute to this significance.
#'  The standarized version of each test statistic and its p-value are computed among other things.
#'
#'
#' @details
#' \itemize{
#' \item{Package: ‘TwoSampleTest.HD’}
#' \item{Version: 1.0}
#' \item{Maintainer: Marta Cousido Rocha \email{martacousido@@uvigo.es}}
#' \item{License: GPL-2}
#' }
#'
#' @return
#' \itemize{
#' \item{‘Equaldis.TStest.HD’}
#' }
#'
#' @author
#' \itemize{
#' \item{Cousido Rocha, Marta.}
#' \item{Soage González, José Carlos.}
#' \item{de Uña-Álvarez, Jacobo.}
#' \item{D. Hart, Jeffrey.}
#' }
#'
#' @references Cousido-Rocha, M., de Uña-Álvarez J., and Hart, J. (2018). A two-sample test for the equality of distributions for high-dimensional data. Preprint.
#'
#'
"_PACKAGE"
#> [1] "_PACKAGE"
