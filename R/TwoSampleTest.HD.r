#' @title A two-sample test for the equality of distributions for high-dimensional data
#' @aliases TwoSampleTest.HD
#' @description Performs the four tests of equality of the p marginal distributions
#'  for two groups proposed in Cousido- Rocha et al.(2018). The methods have been
#'  designed for the low sample size and high dimensional setting. Furthermore,
#'  the possibility that the p variables in each data set can be weakly dependent is considered.
#'  The function also reports a set of p permutation p-values, each of which is derived
#'  from testing the equality of distributions in the two groups for each of the variables
#'  separately. These p-values are useful when a proposed test rejects the global null
#'  hypothesis since it makes it possible to identify which variables have contributed to this significance.
#' @details   The function implements the two-sample tests proposed by Cousido-Rocha, et al. (2018).
#'  The methods “spect”, “boot” and “us” are based on a global statistic which is the average of p
#'  individual statistics corresponding to each of the p variables. Each of these individual statistics measures the difference between
#'  the empirical characteristic functions computed from the two samples.
#'  An alternative expression shows that each statistic is essentially the intergrated
#'  squared difference between kernel density estimates. The global statistic (average) is
#'  standardized using one of three different variance estimators. The method “spect”
#'  uses a variance estimator based on spectral analysis, the method “boot” implements
#'  the block bootstrap to estimate the variance and the method “us” employs a variance estimator
#'  derived from U-statistic theory (more details in Cousido-Rocha et al., 2018).
#'  The methods “spect” and “boot” are suitable under some assumptions including that the
#'  sequence of individual statistics that define the global statistic is strictly stationary,
#'  whereas the method “us” avoids this assumption. However the methods “spect” and “boot”
#'  have been checked in simulations and they perform well even when the stationarity assumption
#'  is violated. The methods “spect” and “us” have their corresponding versions for independent
#'  data (“spect ind” and “us ind”), for which the variance estimators are simplified taking into
#'  acount the independence of the variables. The asymptotic normality (when p tends to infinity)
#'  of the standardized version of the statistic is used to compute the corresponding p-value.
#'  On the other hand, Cousido-Rocha et al. (2018) also proposed the method “perm” whose global
#'  statistic is the average of the permutation p-values corresponding to the individual statistics
#'  mentioned above. This method assumes that the sequence of p-values is strictly stationary,
#'  however in simulations it seems that it performs well when this assumption does not hold.
#'  In addition to providing an alternative global test, these p-values can be used when the
#'  global null hypothesis is rejected and one wishes to identify which of the p variables have
#'  contributed to that rejection. The global statistic depends on a parameter which plays
#'  a role similar to that of a smoothing parameter or bandwidth in kernel density estimation.
#'  For the four global tests this parameter is estimated using the information from all the variables or features.
#'  For the individual statistics from which the permutation
#'  p-values are computed, there are two possibilities: (i) use the value employed in the global test
#'  (b_I.permutation.p.values=“global”), (ii) estimate this parameter for each variable
#'  separately using only its sample information (b_I.permutation.p.values=“individual”)).
#' @param X A matrix where each row is one of the p-samples in the first group.
#' @param Y A matrix where each row is one of the p-samples in the second group.
#' @param method the two-sample test. By default the “us” method is computed. See details.
#' @param I.permutation.p.values  Logical. Default is FALSE. A variable indicating whether to compute the permutation p-values or not when the selected method is not “perm”. See details.
#' @param b_I.permutation.p.values The method used to compute the individual statistics on which are based the permutation p-values. Default is “global”. See details.
#'
#' @return A list containing the following components:
#' \item{standarized statistic: }{the value of the standarized statistic.}
#' \item{p.value: }{the p-value for the test.}
#' \item{statistic: }{the value of the statistic.}
#' \item{variance: }{the value of the variance estimator.}
#' \item{p: }{number of samples or populations.}
#' \item{n: }{sample size in the first group.}
#' \item{m: }{sample size in the second group.}
#' \item{method: }{a character string indicating which two sample test is performed.}
#' \item{I.statistics: }{the p individual statistics.}
#' \item{I.permutation.p.values: }{the p individual permutation p-values.}
#' \item{data.name: }{a character string giving the name of the data.}
#'
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{José Carlos Soage González}
#' \item{Jacobo de Uña-Álvarez}
#' \item{Jeffrey D. Hart}
#' }
#' @references Cousido-Rocha, M., de Uña-Álvarez J., and Hart, J. (2018).
#'  A two-sample test for the equality of marginal distributions for high-dimensional data.
#'  Preprint.
#'
#' @examples
#' \dontshow{
#' # Example
#' ### Data set to check the performance of the code
#'
#' p <- 100
#' n <- 5
#' m <- 5
#'
#' X <- matrix(rnorm(p * n), ncol = n)
#' Y <- matrix(rnorm(p * n), ncol = n)
#' system.time(res <- TwoSampleTest.HD(X, Y, method = "perm"))
#' }
#'
#'
#' \donttest{
#' # We consider a simulated data example. We have simulated the following situation.
#' # We have two groups, for example, 7 patients with tumor 1 and 7 patients with tumor 2.
#' # For each patient 1000 variables are measured, for example, gene expression levels.
#' # Besides, the distributions of 100 of the variables are different in the two groups,
#' # and the differences are in terms of location. The variables are independent to
#' # simplify the generation of the data sets.
#' p <- 1000
#' n = m = 7
#' inds <- sample(1:4, p, replace = TRUE)
#' X <- matrix(rep(0, n * p), ncol = n)
#' for (j in 1:p){
#'   if (inds[j] == 1){
#'     X[j, ] <- rnorm(n)
#'     }
#'   if (inds[j] == 2){
#'     X[j, ] <- rnorm(n, sd = 2)
#'   }
#'   if (inds[j] == 3){
#'     X[j, ] <- rnorm(n, mean = 1)
#'   }
#'   if (inds[j] == 4){
#'     X[j, ] <- rnorm(n, mean = 1, sd = 2)
#'   }
#'  }
#' rho <-  0.1
#' ind <- sample(1:p, rho * p)
#' li <- length(ind)
#' indsy <- inds
#' for (l in 1:li){
#'   if (indsy[ind[l]]==1){
#'      indsy[ind[l]]=3
#'      } else {
#'        if (indsy[ind[l]]==2){
#'           indsy[ind[l]]=4
#'            } else {
#'             if (indsy[ind[l]]==3){
#'                indsy[ind[l]]=1
#'                } else {
#'                 indsy[ind[l]] = 2
#'                 }
#'                }
#'              }
#'            }
#' Y <- matrix(rep(0, m * p), ncol = m)
#' for (j in 1:p){
#'   if (indsy[j] == 1){
#'     Y[j,] <- rnorm(m)}
#'   if (indsy[j] == 2){
#'     Y[j, ] <- rnorm(m, sd = 2)
#'   }
#'   if (indsy[j]==3){
#'     Y[j, ] <- rnorm(m, mean = 1)
#'   }
#'   if (indsy[j] == 4){
#'     Y[j,] <- rnorm(m, mean = 1, sd = 2)
#'   }
#'  }
#'
#' # Our interest is to test the null hypothesis that the distribution of each of the 1000 variables
#' # is the same in the two groups.
#'
#' # We use for this purpose the four methods proposed in Cousido-Rocha et al. (2018).
#'
#' res1 <- TwoSampleTest.HD(X, Y, method = "spect")
#' res1
#' res2 <- TwoSampleTest.HD(X, Y, method = "boot")
#' res2
#' res3 <- TwoSampleTest.HD(X, Y, method = "us")
#' res3
#' res4 <- TwoSampleTest.HD(X, Y, method = "perm")
#' res4
#' # The four methods reject the global null hypothesis.
#' # Hence, we use the individual permutation p-values
#' # to identify which variables are not equally distributed in the two groups.
#' pv<-res4$I.permutation.p.values
#'
#' # Applying a multiple testing procedure to these p-values
#' # we can detect the variables with different distributions for the two groups.
#' # The following plot of the individual permutation p-values is also informative.
#' # We remark in red the 100 smallest p-values.
#'
#' pv_sort <- sort(pv)
#' cri <- pv_sort[100]
#' ind <- which(pv <= cri)
#' plot(1:p, pv, main = "Individual permutation p-values",
#'      xlab = "Variables", ylab = "p-values")
#' points(ind, pv[ind], col = "red")
#' }
#' @importFrom utils combn
#' @export


################################################################################
# Two-sample problem
################################################################################

TwoSampleTest.HD <- function(X, Y, method = c("spect", "spect_ind", "boot", "us", "us_ind", "perm"),
                             I.permutation.p.values = FALSE, b_I.permutation.p.values = c("global", "individual")) {

  cat("Call:", "\n")
  print(match.call())

  if (missing(method)) {
    method <- "us"
    cat("'us' method used by default\n")
  }

  if (missing(b_I.permutation.p.values)) {
    b_I.permutation.p.values <- "global"
    cat("'global' bandwith used by default\n")
  }

  method <- match.arg(method)
  DNAME  <- deparse(substitute(c(X, Y)))
  METHOD <- "A two-sample test for the equality of distributions for high-dimensional data"

  match.arg(method)
  match.arg(b_I.permutation.p.values)

  p <- nrow(X)
  n <- ncol(X)
  m <- ncol(Y)

  c1 <- 1 / p
  c2 <- 1.114 * mean(c(n, m)) ^ (-1 / 5)

  sa <- apply(X, 1, var)
  sb <- apply(Y, 1, var)

  si <- ((n - 1) * unlist(sa) + (m - 1) * unlist(sb)) / (n + m - 2)

  spool <- sqrt(c1 * sum(si))
  h <- spool * c2

  y <- c(rep(1, 5), rep(2, 5))


  #=============================================================================
  ## Functions to compute the statistic
  #=============================================================================

  ### Method to choose m
  mval <- function(rho, lagmax, kn, rho.crit) {

    ## Compute the number of insignificant runs following each rho(k),
    ## k = 1, ..., lagmax.
    num.ins <- sapply(1:(lagmax - kn + 1),
                      function(j) sum((abs(rho) < rho.crit)[j:(j + kn - 1)]))

    ## If there are any values of rho(k) for which the kn proceeding
    ## values of rho(k + j), j = 1, ..., kn are all insignificant, take the
    ## smallest rho(k) such that this holds (see footnote c of
    ## Politis and White for further details).
    if(any(num.ins == kn)){
      return(which(num.ins == kn)[1])
    } else {
      ## If no runs of length kn are insignificant, take the smallest
      ## value of rho(k) that is significant.
      if(any(abs(rho) > rho.crit)) {
        lag.sig <- which(abs(rho) > rho.crit)
        k.sig <- length(lag.sig)

        if(k.sig == 1) {
          ## When only one lag is significant, mhat is the sole
          ## significant rho(k).
          return(lag.sig)
        } else {
          ## If there are more than one significant lags but no runs
          ## of length kn, take the largest value of rho(k) that is
          ## significant.
          return(max(lag.sig))
        }
      } else {

        ## When there are no significant lags, mhat must be the
        ## smallest positive integer (footnote c), hence mhat is set
        ## to one.
        return(1)
      }
    }
  }

  Lval <- function(x, method = mean) {

    x <- matrix(x)
    n <- nrow(x)
    d <- ncol(x)

    ## Parameters for adapting the approach of Politis and White (2004)
    kn <- max(5, ceiling(log10(n)))
    lagmax <- ceiling(sqrt(n)) + kn
    rho.crit <- 1.96 * sqrt(log10(n) / n)

    m <- numeric(d)
    for (i in 1:d) {
      rho <- stats::acf(x[,i], lag.max = lagmax, type = "correlation", plot = FALSE)$acf[-1]
      m[i] <- mval(rho, lagmax, kn, rho.crit)
    }
    return(2 * method(m))
  }


  ##########################################
  ### Variance estimators (stationary) #####
  ##########################################

  ### Spectral variance estimator
  variance_spectral <- function(J) {

    part2 <- 0
    k <- Lval(matrix(J), method = min)
    c <- stats::acf(J, type = "covariance", lag.max = k, plot = FALSE)$acf
    c0 <- c[1]
    c  <- c[-1]
    for (i in 1:k) {
      part2 <- part2 + (1 - (i / (k + 1))) * c[i]
    }

    statistic <- c0 + 2 * part2
    return(statistic)
  }

  # Particular case of independence
  variance_spectral_ind <- function(J) {
    k <- 1
    c <- stats::acf(J, type = "covariance", lag.max = k, plot = FALSE)$acf
    statistic <- c[1]
    return(statistic)
  }

  ### Variance block bootstrap
  variance <- function(pv) {

    p <- length(pv)
    # print(h)
    # h = 0.8
    k <- Lval((pv), method = min)
    # bootstats = 1:(p - k + 1)
    # bootstats[1] = sum(J[1:k])
    # for(j in 1:(p - k)) {
    # bootstats[j + 1] = bootstats[j] - J[j] + J[j + k]
    # }
    bootstats <- 1:(p / k)
    for (j in 1:(p / k)) {
      bootstats[j] <- sum(pv[(k * (j - 1) + 1):(j * k)])
    }
    varest <- var(bootstats / sqrt(k))
    return(varest)
  }


  #############################################
  ### Variance estimator (Non-stationary)  ####
  #############################################

  ### Variance estimator based on serfling and time series
  VarJhat <- function(x, y, h) {

    E1 <- E1hat(c(x, y), h)
    E2 <- E2hat(c(x, y), h)
    E3 <- E3hat(c(x, y), h)
    n  <- length(x)
    m  <- length(y)
    c1 <- (n - 2) * (n - 3) / (n * (n - 1))
    c1 <- c1 + (m - 2) * (m - 3) / (m * (m - 1))
    c1 <- c1 + 4 * (n - 1) * (m - 1) / (n * m)
    c1 <- c1 + 2
    c1 <- c1 - 4 * (n - 2) / n
    c1 <- c1 - 4 * (m - 2) / m
    c2 <- -4 / (n * (n - 1)) - 4 / (m * (m - 1)) - 8 / (m * n)
    c3 <- 2 / (n * (n - 1)) + 2 / (m * (m - 1)) + 4 / (m * n)
    vhat <- c1 * E1 + c2 * E2 + c3 * E3
    return(vhat * (4 * pi * h ^ 2))
  }

  E1hat <- function(Z, h) {

    N <- length(Z)
    Zmat <- t(matrix(Z, N, N))
    Zmat <- Z - Zmat
    Zmat <- stats::dnorm(Zmat, sd = sqrt(2) * h)
    diag(Zmat) <- 0
    E1 <- 0
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        t1 <- stats::dnorm(Z[i] - Z[j], sd = sqrt(2) * h)
        t1 <- sum(t1 * Zmat) - sum(t1 * Zmat[i,] + t1 * Zmat[j,] + t1 * Zmat[, i] +
                                     t1 * Zmat[, j]) + 2 * t1 * Zmat[i, j]
        E1 <- E1 + t1
      }
    }
    Number <- N * (N - 1) * (N - 2) * (N - 3) / 2
    E1 / Number
  }


  E2hat <- function(Z, h) {
    N  <- length(Z)
    E2 <- 0
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        t1  <- stats::dnorm(Z[i] - Z[j], sd = sqrt(2) * h)
        vec <- setdiff(1:N, c(i, j))
        t2  <- stats::dnorm(Z[i] - Z[vec], sd = sqrt(2) * h)
        t3  <- stats::dnorm(Z[j] - Z[vec], sd = sqrt(2) * h)
        t1  <- sum(t1 * t2) + sum(t1 * t3)
        E2  <- E2 + t1
      }
    }
    Number <- N * (N - 1) * (N - 2)
    E2 / Number
  }


  E3hat <- function(Z, h) {
    N  <- length(Z)
    E3 <- 0
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        t1 <- stats::dnorm(Z[i] - Z[j], sd = sqrt(2) * h) ^ 2
        E3 <- E3 + t1
      }
    }
    Number <- N * (N - 1) / 2
    E3 / Number
  }


  ### Dirichlet kernel
  variance_est_Dirichlet <- function(X, Y, h, J) {
    VJ <- 1:p

    for (i in 1:p) {
      VJ[i] <- VarJhat(X[i,], Y[i,], h)
    }

    varest1 <- sum(VJ)
    part2   <- 0
    k <- Lval(matrix(J), method = min)
    c <- 1:k

    for (i in 1:k) {
      sum <- 0
      for (j in 1:(p - i)) {
        sum <- sum + J[j] * J[j + i]
      }
      c[i] <- sum
    }
    for (i in 1:k) {
      part2 <- part2 + c[i]
    }

    statistic <- (1 / p) * (varest1) + (2 / p) * part2
    return(statistic)
  }


  ### Particular case of independence
  variance_est_ind <- function(X, Y, h) {
    VJ <- 1:p

    for (i in 1:p) {
      VJ[i] <- VarJhat(X[i,], Y[i,], h)
    }

    varest1 <- sum(VJ)
    statistic <- (1 / p) * varest1
    return(statistic)
  }


  ### Function Ji
  Ji <- function(x, y, h) {

    n <- length(x)
    m <- length(y)
    X <- t(matrix(x, n, n))
    X <- x - X
    X <- stats::dnorm(X, sd = sqrt(2) * h)
    diag(X) <- 0
    t1 <- sum(X) / (n * (n - 1))
    Y  <- t(matrix(y, m, m))
    Y  <- y - Y
    Y  <- stats::dnorm(Y, sd = sqrt(2) * h)
    diag(Y) <- 0
    t2 <- sum(Y) / (m * (m - 1))
    XY <- x[1] - y
    for (j in 2:n) {
      XY <- rbind(XY, x[j] - y)
    }
    XY <- stats::dnorm(XY, sd = sqrt(2) * h)
    t3 <- mean(XY)
    (t1 + t2 - 2 * t3) * (sqrt(2 * pi * 2 * h ^ 2))
  }


  A <- function(X, Y, h) {
    p <- nrow(X)
    Jivalue <- 1:p
    for (i in 1:p) {
      Jivalue[i] <- Ji(X[i, ], Y[i, ], h)
    }
    return(Jivalue)
  }


  # Here we compute the statistic Ji, it is necessary to show it in the package
  J <- A(X, Y, h)

  suma <- function(J) {
    p <- length(J)
    sum(J) / sqrt(p)
  }


  # The test statistic before the standardization (also necessary to report in the package)
  e <- sum(J) / sqrt(p)

  if (method == "boot") {
    ### Bootstrap
    ### Variance of the statistic (save in memory and access to it)
    var <- variance(J)

    ### Standarized test statistic using the bootstrap method (show)
    s <- e / sqrt((var))

    ### Corresponding p-value (show)
    pvalor <- 1 - stats::pnorm(s)
    s
  }

  if(method == "spect") {
    ### Spectral

    ### Variance of the statistic (save in memory and access to it)
    var_spectral <- variance_spectral(J)

    ### Standarized test statistic using the spectral method (show)
    s_spectral <- (e) / sqrt((var_spectral))

    ### Corresponding p-value (show)
    pvalor_spectral <- 1 - stats::pnorm(s_spectral)

    s_spectral
  }

  if(method == "us" | method == "us_ind") {
    ### Serfling Dirichlet

    ### Variance of the statistic (save in memory and access to it)
    var_est_Dirichlet <- variance_est_Dirichlet(X, Y, h, J)

    ### Standarized test statistic using the Serfling Dirichlet method (show)
    s_est_Dirichlet <- e / sqrt((var_est_Dirichlet))

    ### Corresponding p-value (show)
    pvalor_Dirichlet <- 1 - stats::pnorm(s_est_Dirichlet)

    s_est_Dirichlet
  }

  #-----------------------------------------------------------------------------
  # Versions for independent data
  #-----------------------------------------------------------------------------

  if (method == "us_ind") {
    ### Serfling for independent data (us_ind)

    var_est_ind <- variance_est_ind(X, Y, h) ### save in memory
    s_est_ind <- e / sqrt(unlist(var_est_ind)) ### show
    pvalor_est_ind <- 1 - stats::pnorm(s_est_Dirichlet) ### show
  }

  if (method == "spect_ind") {
    ### Spectral for independent data  (spect_ind)
    var_spectral_ind <- variance_spectral_ind(J)### save in memory

    s_spectral_ind <- e / sqrt((var_spectral_ind)) ### show

    pvalor_spectral_ind <- 1 - stats::pnorm(s_spectral_ind)  ### show
  }

  ### Furthermore than these three methods we also reported one based on p-values
  ### instead of Ji statistics


  ###############################
  ### Permutation test        ###
  ###############################

  if(method == "perm" | I.permutation.p.values == TRUE) {

    # if(method != "perm" & I.permutation.p.values == TRUE){
    #   print("OK")
    # }
    #
    # if(method != "perm" & I.permutation.p.values == FALSE){
    # } else {
    # if(method == "perm" | I.permutation.p.values == TRUE) {

    kern.permute <- function(x, y, b) {
      m <- length(x)
      n <- length(y)
      S <- 1:(m + n)
      stat <- Ji(x, y, b)
      X <- c(x, y)
      if (m == n) {
        M <- utils::combn(1:(2 * m), m)
        Ml <- ncol(M) / 2
        stats <- 1:Ml
        for (j in 1:Ml) {
          xstar <- X[M[, j]]
          ystar <- X[M[, 2 * Ml + 1 - j]]
          stats[j] <- Ji(xstar, ystar, b)
        }
      }
      if (m != n) {
        M  <- utils::combn(S, m)
        Ml <- ncol(M)
        stats <- 1:Ml
        for (j in 1:Ml) {
          xstar <- X[M[, j]]
          S1 <- setdiff(S, M[, j])
          ystar <- X[S1]
          stats[j] <- Ji(xstar, ystar, b)
        }
      }
      vec <- (1:Ml)[stats >= stat]
      list(stat, stats, length(vec) / Ml)
    }

    kern.permute_1 <- function(x, y, b) {
      kern.permute(x, y, b)[3]
    }

    if(b_I.permutation.p.values == "global"){
      permutation_diff_ind <- function(X, Y, h) {
        p <- nrow(X)
        o <- rep(1, p)
        for (i in 1:p) {
          o[i] <- as.numeric(kern.permute_1(X[i,], Y[i,], h))
        }
        return(o)
      }

      ### P-values corresponding to each null hypothesis (in memory and access to them)

      pv <- permutation_diff_ind(X, Y, h)
    }


    if(b_I.permutation.p.values == "individual") {

      ### The p-values computing with indvidual bandwidth
      sa <- apply(X, 1, var)
      sb <- apply(Y, 1, var)
      si <- ((n - 1) * unlist(sa) + (m - 1) * unlist(sb)) / (n + m - 2)
      h  <- sqrt(si) * c2


      permutation_diff_ind <- function(X, Y, h) {
        p <- nrow(X)
        o <- rep(1, p)
        for (i in 1:p) {
          o[i] <- as.numeric(kern.permute_1(X[i,], Y[i,], h[i]))
        }
        return(o)
      }

      ### P-values corresponding to each null hypothesis (in memory and access to them)
      pv <- permutation_diff_ind(X, Y, h)
    }



    ### Graph in memory and acess if it selected
    # graphics::plot(1:p, pv, main = "Individual p-values", xlab = "Number of null hypothesis",
    #      ylab = "p-values")


    ############################
    ### Variance estimators ####
    ############################

    variance_spectralR <- function(J) {
      part2 <- 0
      part3 <- 0
      k <- Lval(matrix(J), method = min)
      c <- stats::acf(J, type = "covariance", lag.max = k, plot = FALSE)$acf
      c0 <- c[1]
      c <- c[-1]
      for (i in 1:k) {
        part2 <- part2 + cos(pi * i / (k + 1)) * c[i]
        part3 <- part3 + c[i]

      }

      statistic <- c0 + 2 * part2
      stat <- c0 + 2 * part3
      return(max(statistic, stat))
    }


    ### Variance of the statistic
    var_pv_sR <- variance_spectralR(pv)


    ### Non standarized statistic
    pv_ <- mean(pv)
    N <- n + m
    if (m == n) {
      nm <- ((factorial(N)) / (factorial(n) * factorial(N - n))) / 2
    } else {
      nm <- ((factorial(N)) / (factorial(n) * factorial(N - n)))
    }
    ### Number of possible permutations
    nm
    mean_p <- ((nm + 1) / nm) * 0.5

    ### Standarized statistic (show)
    s_sR <- ((pv_ - mean_p) * sqrt(p)) / sqrt((var_pv_sR))

    ### Corresponding p-value
    pvalor_s_sR <- stats::pnorm(s_sR)

  }

  statistic <- switch(method, spect = s_spectral, spect_ind = s_spectral_ind,
                      boot = s, us = s_est_Dirichlet, us_ind = s_est_ind, perm = s_sR)
  names(statistic) <- "standarized statistic"

  statistic2 <- switch(method, spect = e, spect_ind = e,  boot = e, us = e, us_ind = e,
                       perm = (pv_)*sqrt(p))

  p.value <- switch(method, spect = pvalor_spectral, spect_ind = pvalor_spectral_ind,
                    boot = pvalor, us = pvalor_Dirichlet, us_ind = pvalor_est_ind,
                    perm = pvalor_s_sR)

  met <- switch(method, spect = "spect", spect_ind = "spect_ind",  boot = "boot",
                us = "us", us_ind = "us_ind", perm = "perm")

  variance <- switch(method, spect = var_spectral, spect_ind = var_spectral_ind,
                     boot = var, us = var_est_Dirichlet, us_ind = var_est_ind,
                     perm = var_pv_sR)

  RVAL <- list(statistic = statistic, p.value = p.value, method = METHOD,
               data.name = DNAME, sample.size = n, method1 = met)


  if(method == "perm"){
    RVAL2 <- list(standarized.statistic = statistic, p.value = p.value,
                  statistic = (pv_)*sqrt(p), variance = variance, p = p, n = n, m = m,
                  method = met, I.statistics = J, I.permutation.p.values = pv,
                  data.name = DNAME)
  }

  if (method != "perm" & I.permutation.p.values == FALSE) {
    RVAL2 <- list(standarized.statistic = statistic, p.value = p.value,
                  statistic = e, variance = variance, p = p, n = n, m = m,
                  method = met, I.statistics = J, data.name = DNAME)
  }


  if (method != "perm" & I.permutation.p.values == TRUE) {
    RVAL2 <- list(standarized.statistic = statistic, p.value = p.value,
                  statistic = e, variance = variance, p = p, n = n, m = m,
                  method = met, I.statistics = J, I.permutation.p.values = pv,
                  data.name = DNAME)
  }

  class(RVAL) <- "htest"

  print(RVAL)
  return(invisible(RVAL2))
}




# Example
#
# ### Data set to check the performance of the code
#
# p <- 100
# n <- 5
# m <- 5
#
# X <- matrix(rnorm(p * n), ncol = n)
# Y <- matrix(rnorm(p * n), ncol = n)
# system.time(res <- TwoSampleTest.HD(X, Y, method = "perm"))
