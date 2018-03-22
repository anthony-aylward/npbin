#===============================================================================
# npbin.R
#===============================================================================

# functions for non-parametric binomial imbalance testing




# Imports ----------------------------------------------------------------------

#' @import data.table
#' @import parallel
#' @import nloptr
#' @import VGAM




# Function definitions ---------------------------------------------------------

#' @title Beta binomial density.
#'
#' @details
#' \code{dbetabinom.vec} returns a vector giving the beta binomial density.
#'
#' @param x Vector of quantiles
#' @param m Vector giving numbers of trials
#' @param rate1 First shape parameter
#' @param rate2 Second shape parameter
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return vector describing the density
#' @export
#' @seealso \code{\link{estNull2}}, \code{\link{betaTrim_mle}}
dbetabinom.vec <- function(x, m, rate1, rate2, log=TRUE) {
  n <- length(x)
  lvec <- sapply(
    1:n,
    function(ii) dbetabinom.ab(
      x[ii],
      size = m[ii],
      shape1 = rate1,
      shape2 = rate2,
      log = log
    )
  )
  unlist(lvec)
}

#' @title Create B-spline.
#'
#' @details
#' Create a B-spline.
#'
#' @param xrange Support of the spline
#' @param breaks Break points of the spline
#' @param k Spline order
#' @param ncores Number of cores to use
#' @return A list containing breaks and spline coefficients
#' @seealso \code{\link{bsplfun.updt}}, \code{\link{emBinBspl}},
#'   \code{\link{emBspl}}
bsplfun <- function(
  xrange = c(0, 1),
  breaks = seq(xrange[1], xrange[2], length.out = 100),
  k = 4,
  ncores = 1
) {
  nbreaks <- length(breaks)
  nbasis <- nbreaks + k - 2
  step.break <- mean(abs(diff(breaks)))
  breaks.extend <- c(
    breaks[1] - step.break * ((k-1):1),
    breaks,
    breaks[nbreaks] + step.break * (1:(k-1))
  )
  coef0 <- matrix(1,1,1)
  bspl <- lapply(
    1:length(breaks.extend),
    function(ii) list(breaks=breaks.extend[ii:(ii+1)],coef=coef0)
  )
  for(kk in 2:k){
    bspl.kk <- mclapply(
      1:(length(bspl)-1),
      function(ii) bsplfun.updt(ii, bspl),
      mc.cores = ncores,
      mc.preschedule = FALSE
    )
    bspl <- bspl.kk
  }
   out <- bspl[1:nbasis]
   for(ii in 1:(k-1)){
      outii <- out[[ii]]
      breaksii <- outii[["breaks"]][(k+1-ii):(k+1)]
      coefii <- outii[["coef"]][,(k-ii+1):k]
      if(ii==1){coefii <- matrix(coefii,ncol=1)}
      out[[ii]] <- list(breaks=breaksii,coef=coefii)
    }
   for(ii in 1:(k-1)){
      outii <- out[[nbasis+1-ii]]
      breaksii <- outii[["breaks"]][1:(ii+1)]
      coefii <- outii[["coef"]][,1:ii]
      if(ii==1){coefii <- matrix(coefii,ncol=1)}
      out[[nbasis+1-ii]] <- list(breaks=breaksii,coef=coefii)
    }
   for(ii in 1:nbasis){
      outii <- out[[ii]]
      breaksii <- outii[["breaks"]]
      coefii <- outii[["coef"]]
      intii <- sum(
        sapply(
          1:(length(breaksii) - 1),
          function(jj) sum(
            (breaksii[jj + 1]^(1:k) - breaksii[jj]^(1:k)) * coefii[,jj] / (1:k)
          )
        )
      )
      coefii <- coefii/intii
      out[[ii]] <- list(breaks=breaksii,coef=coefii)
    }
   out
}

#' @title Update B-spline.
#'
#' @details
#' Updates the B-spline at a certain index and returns a new spline.
#'
#' @param ii Index to update
#' @param bspl Spline to update
#' @return A list containing breaks and coefficients for a new spline
#' @seealso \code{\link{bsplfun}}
bsplfun.updt <- function(ii, bspl) { 
  bs1 <- bspl[[ii]]
  breaks1 <- bs1[["breaks"]]
  coef1 <- bs1[["coef"]]
  bs2 <- bspl[[ii+1]]
  breaks2 <- bs2[["breaks"]]
  coef2 <- bs2[["coef"]]
  kk <- length(breaks1)       # we assume that length(breaks1) and length(breaks2) are the same
  breaks <- c(breaks1[1], breaks2) # we assume that breaks1[-1] == breaks2[1:(kk-1)]
  coef <- (
    rbind(0, cbind(coef1, 0)) / (breaks1[kk] - breaks1[1])
    - rbind(cbind(coef1, 0), 0) * breaks1[1] / (breaks1[kk] - breaks1[1])
    - rbind(0, cbind(0, coef2)) / (breaks2[kk] - breaks2[1])
    + rbind(cbind(0, coef2), 0) * breaks2[kk] / (breaks2[kk] - breaks2[1])
  )
  list(breaks = breaks, coef = coef)
}

#' @title Spline basis function integration.
#'
#' @details
#' Integration over one spline basis function
#'
#' @param bk Pair of points
#' @param x x
#' @param m m
#' @param cf Coefficient
#' @seealso \code{\link{intBinBspl}}
iBiBsFun <- function(bk, x, m, cf) {
  k <- length(cf)
  sapply(
    1:k,
    function(jj) (
      cf[jj]
      * exp(lgamma(m + 1) - lgamma(x + 1) + lgamma(x + jj) - lgamma(m + jj + 1))
      * (
        pbeta(bk[2], shape1 = jj + x, shape2 = m - x + 1)
        - pbeta(bk[1], shape1 = jj + x, shape2 = m - x + 1)
      )
    )
  )
}

#' @title Binomial and B-spline integration.
#'
#' @details
#' Integration of Binomial + B spline.
#'
#' @param bs A B-spline
#' @param x x
#' @param m m
#' @return The integral value.
#' @seealso \code{\link{iBiBsFun}}, \code{\link{getDesignMtx}}
intBinBspl <- function(bs, x, m) { 
  breaks <- bs[['breaks']]
  cf <- bs[['coef']]
  out <- sapply(
    1:(length(breaks) - 1),
    function(ii) sum(iBiBsFun(breaks[ii:(ii + 1)], x, m, cf[, ii]))
  )
  sum(out)
}

#' @title Design matrix extraction
#'
#' @details
#' get the design matrix c_ij
#'
#' @param bs A B-spline
#' @param x x
#' @param m m
#' @param ncores Number of cores to use
#' @seealso \code{\link{ntBinBspl}}, \code{\link{emBinBspl}},
#'   \code{\link{emBspl}}
getDesignMtx <- function(bs, x, m, ncores = 1) {
  nbasis <- length(bs)
  ndt <- length(x)
  out <- mclapply(
    1:ndt,
    function(ii) sapply(bs, function(jj) intBinBspl(jj, x[ii], m[ii])),
    mc.cores = ncores,
    mc.preschedule=FALSE
  )
  out <- t(as.matrix(as.data.frame(out)))
  dimnames(out) <- NULL
  out
}

#' @title NPBin model.
#'
#' @details
#' We remove the basis that only covers bin at the boundary
#'
#' @param x x
#' @param m m
#' @param breaks break points for the spline
#' @param k Spline order
#' @param ncores Number of cores to use
#' @param err.max Max error for the EM algorithm
#' @param iter.max Max number of iterations
#' @return A list giving model details
#' @export
#' @seealso \code{\link{bsplfun}}, \code{\link{getDesignMtx}},
#'   \code{\link{npbin}}
emBinBspl <- function(
  x,
  m,
  breaks = seq(0, 1, length.out = 101),
  k = 4,
  pi.init = rep(1, length(breaks) + k - 4) / (length(breaks) + k - 4),
  ncores = 1,
  err.max = 0.000001,
  iter.max = 200
) {
  n <- length(x)
  bspl <- bsplfun( ## create the break points and coefficients of the B spline functions
    range(breaks),
    breaks = breaks,
    k = k,
    ncores = ncores
  )
  bspl <- bspl[-c(1, length(bspl))]
  nb <- length(bspl)
  dmtx <- getDesignMtx(bspl, x, m, ncores = ncores) ### get the design matrix
  err <- Inf
  iter <- 0
  ll.init <- -Inf
  ll.all <- ll.init
  err.all <- err
  while (err > err.max & iter < iter.max) {
    dtot <-  dmtx %*% matrix(pi.init, length(pi.init), 1) 
    post <- (
      dmtx
      * (matrix(pi.init, n, nb, byrow = TRUE))
      / matrix(dtot, n, nb, byrow = FALSE)
    )
    pi <- colSums(post)
    pi <- pi / sum(pi)
    ll <- sum(apply(dmtx, 1, function(ii) log(sum(ii * pi))))
    err <- ll - ll.init
    err <- max(err, err / abs(ll))
    ll.all <- c(ll.all, ll)
    err.all <- c(err.all, err)
    pi.init <- pi
    ll.init <- ll
    iter <- iter + 1
  } 
  list(
    pi = as.numeric(pi.init),
    post = post,
    bspl = bspl,
    dmtx = dmtx,
    f = as.numeric(dmtx %*% pi.init),
    ll.all = ll.all,
    err.all = err.all,
    convergence = list(
      err = err,
      err.max = err.max,
      converged = (err < err.max),
      niter = iter,
      ll = ll.init
    ),
    controls = list(k = k, nbasis = nb, breaks = breaks)
  )
}

#' @title B-spline evaluation.
#'
#' @details
#' evaluate B-spline function
#'
#' @param p p
#' @param bspl A B-spline
#' @param ncores Number of cores to use
#' @return A matrix giving the evaluation result
#' @seealso \code{\link{evBsplDrv}}, \code{\link{estNull1}},
#'   \code{\link{emBspl}}
evBspl <- function(p, bspl, ncores = 1) {
  nb <- length(bspl)
  k <- dim(bspl[[1]][["coef"]])[1]
  evBspl.sg <- function(p, bsplsg) {
    coef <- bsplsg[["coef"]]
    bk <- bsplsg[["breaks"]]
    id <- which(p < bk)
    if (length(id) == 0 | id[1] == 1) {
      out <- 0
    } else {
      kk <- id[1] - 1
      out <- sum(sapply(1:k, function(jj) coef[jj, kk] * p^(jj - 1)))
    }
    out
  }
  out <- mclapply(
    p,
    function(pp) sapply(bspl, function(ii) evBspl.sg(pp, ii)),
    mc.cores = ncores,
    mc.preschedule = FALSE
  )
  out <- t(as.matrix(as.data.frame(out)))
  dimnames(out) <- NULL
  out
}

#' @title B-spline derivative evaluation.
#'
#' @details
#' evaluate the derivative of B-spline
#'
#' @param p p
#' @param bspl A B-spline
#' @param ncores Number of cores to use
#' @return A matrix giving the evaluation result
#' @seealso \code{\link{evBspl}}, \code{\link{estNull1}}
evBsplDrv <- function(p, bspl, ncores = 1) {
  nb <- length(bspl)
  k <- dim(bspl[[1]][["coef"]])[1]
  evBsplDrv.sg <- function(p, bsplsg){
    coef <- bsplsg[["coef"]]
    bk <- bsplsg[["breaks"]]
    id <- which(p < bk)
    if (length(id) == 0 | id[1] == 1) {
      out <- 0
    } else {
      kk <- id[1] - 1
      out <- sum(sapply(2:k, function(jj) (jj - 1) * coef[jj, kk] * p^(jj - 2)))
    }
    out
  }
  out <- mclapply(
    p,
    function(pp) sapply(bspl, function(ii) evBsplDrv.sg(pp, ii)),
    mc.cores = ncores,
    mc.preschedule = FALSE
  )
  out <- t(as.matrix(as.data.frame(out)))
  dimnames(out) <- NULL
  out
}

#' @title Null model estimation
#'
#' @details
#' estimate the null model, part 1
#'
#' @param mod mod
#' @param pseq pseq
#' @param ncores Number of cores to use
#' @return List describing the model
#' @export
#' @seealso \code{\link{evBspl}}, \code{\link{evBsplDrv}}
#' @family estNulls
estNull1 <- function(
  mod,
  pseq = (1:9999) / 1e4,
  ncores = 1
) {
  pseq <- (1:9999) / 10000
  bspl <- mod[["bspl"]]
  pi <- mod[["pi"]]
  drvbspl <- evBsplDrv(pseq, bspl, ncores = ncores)
  df <- apply(drvbspl, 1, function(ii) sum(ii * pi))
  evbspl <- evBspl(pseq, bspl, ncores = ncores)
  fseq <- apply(evbspl, 1, function(ii) sum(ii * pi))
  list(fseq = fseq, df = df)
}

#' @title Null model estimation
#'
#' @details
#' estimate the null model, part 2
#'
#' @param x x
#' @param m m
#' @param mod mod
#' @param prep Info prepared by estNull1
#' @param init init
#' @param iter.max Maximum number of iterations for EM algorithm
#' @param err.max err.max
#' @param algorithm Numerical library to use
#' @param pseq pseq
#' @param lb Lower bounds
#' @param ub Upper bounds
#' @return mod
#' @export
#' @seealso \code{\link{dbetabinom.vec}}
#' @family estNulls
estNull2 <- function(
  x,
  m,
  mod,
  prep,
  init = NULL,
  iter.max = 200,
  err.max = 1e-6,
  algorithm = "NLOPT_GN_DIRECT_L",
  pseq = (1:9999)/1e4,
  lb = c(0, 0),
  ub = rep(log(1e4), 2)
) {
  if (is.null(init)) {
    phat <- x / m
    m1p <- mean(phat)
    n <- length(x)
    ss <- n * var(phat)
    mi <- sum(1 / m)
    m2p <- (ss - (m1p - m1p^2) * mi) / (n - mi)
    sc <- m1p^2 / m2p - 1
    init <- c(m1p * sc, (1 - m1p) * sc)
  }
  fseq <- prep[["fseq"]]
  df <- prep[["df"]]
  ell <- function(ipt){
    shape1 <- exp(ipt[1])
    shape2 <- exp(ipt[2])
    df0 <- (
      (shape1 - 1 - (shape1 + shape2 - 2) * pseq)
      * dbeta(pseq, shape1 - 1, shape2 - 1, log = FALSE)
      * exp(lbeta(shape1 - 1, shape2 - 1) - lbeta(shape1, shape2))
    )
    f0seq <- dbeta(pseq, shape1, shape2, log = FALSE)
    ovec <- (fseq * df0 - df * f0seq)^2 / fseq^3 ##  
    ovec[is.na(ovec)] <- 0
    mean(ovec)
  }
  optout <- nloptr(
    log(init),
    ell,
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = algorithm,
      maxeval = iter.max,
      ftol_rel = err.max,
      xtol_rel = sqrt(err.max)
    )
  )
  coef.opt <- exp(optout[["solution"]])
  f <- mod[["f"]]
  f0 <- dbetabinom.vec(x, m, coef.opt[1], coef.opt[2], log = FALSE)
  pi0 <- min(1, 1 / quantile(f0 / f, probs = 0.975))
  maxi <- c(x[which.min(f / f0)], m[which.min(f / f0)])
  coef <- list(shape1 = coef.opt[1], shape2 = coef.opt[2], pi0 = pi0)
  out <- mod
  out[["coef.null"]] <- coef
  out[["pi0"]] <- pi0
  out[["f0"]] <- f0
  out[["locfdr"]] <- pi0*f0/f
  out[["convergence.null"]] <- list(opt.out=optout,p.maxlr=maxi)
  out
}

#' @title Null model estimation
#'
#' @details
#' estimate the null model
#'
#' @param x x
#' @param m m
#' @param mod mod
#' @param prep Info prepared by estNull1
#' @param init init
#' @param iter.max Maximum number of iterations for EM algorithm
#' @param err.max err.max
#' @param algoritm Numerical library to use
#' @param pseq pseq
#' @param lb Lower bounds
#' @param ub Upper bounds
#' @return mod
#' @export
#' @seealso \code{\link{npbin}}
#' @family estNulls
estNull  <- function(
  x,
  m,
  mod,
  init = NULL,
  ncores = 1,
  iter.max = 200,
  err.max = 1e-6,
  algorithm = "NLOPT_GN_DIRECT_L",
  pseq = (1:9999) / 1e4,
  lb = c(0,0),
  ub = rep(log(1e4),2)
) {
  prep <- estNull1(mod, pseq, ncores = ncores)
  estNull2(
    x,
    m,
    mod,
    prep,
    init = init,
    iter.max = iter.max,
    err.max = err.max,
    algorithm = algorithm,
    pseq = pseq,
    lb = lb,
    ub = ub
  )
}

#' @title Convert locfdr to FDR.
#'
#' @details
#' Convert locfdr to FDR.
#'
#' @param locfdr Local FDR value
#' @return FDR value
#' @export
#' @seealso \code{\link{npbin}}
locfdr2FDR <- function(locfdr) {
  n <- length(locfdr)
  sapply(1:n, function(ii) mean(locfdr[locfdr <= locfdr[ii]]))
}

#' @title Convert ranking to the number of discoveries
#'
#' @details
#' Convert ranking to the number of discoveries
#'
#' @param r Vector of ranks
#' @param id id
#' @return vector indicating discoveries
rank2nhit <- function(r, id) {sapply(r, function(y) sum((r <= y) & id))}

#' @title NPBin
#'
#' @details
#' Perform an NPBin analysis
#'
#' @param dt.ct Data table containing allele counts
#' @param n_breaks Number of spline breaks
#' @param spline_order Spline order
#' @param pi_init Initial weights
#' @param n_cores Number of cores to use
#' @return Data table containing model information
#' @export
#' @seealso \code{\link{initialize_weights}}, \code{\link{emBinBspl}},
#'   \code{\link{estNull}}
npbin <- function(dt.ct, n_breaks, spline_order, pi_init, n_cores) {
  n <- nrow(dt.ct)
  breaks <- seq(0, 1, length.out = n_breaks)
  
  # estimate the overall model
  overall_model_estimate <- emBinBspl(
    dt.ct[, xm],
    dt.ct[, m],
    breaks = breaks,
    k = spline_order,
    pi.init = pi_init,
    ncores = n_cores,
    err.max = 1e-3,
    iter.max = 200
  )  

  # estimate the null model
  null_model_estimate <- estNull(
    dt.ct[, xm],
    dt.ct[, m],
    overall_model_estimate,
    init = NULL,
    iter.max = 200,
    ncores = n_cores,
    ub = rep(log(1e4), 2),
    err.max = 1e-4
  )
  
  dt.ct[,
    fnp := null_model_estimate[["f"]]
  ][,
    f0np := null_model_estimate[["f0"]]
  ][,
    locfdrnp := null_model_estimate[["locfdr"]]
  ][,
    fdrnp := locfdr2FDR(locfdrnp)
  ][,
    ranknp := rank(locfdrnp, ties.method = "max")
  ]
  
  dt.ct
}




# functions for methods to be compared with ------------------------------------

#' @title Efron null model
#'
#' @details
#' Estimate the null model following Efron's approach.
#'
#' @param x x
#' @param m m
#' @param p p
#' @param pct0 pct0
#' @param init init
#' @param iter.max Maximum number of iterations
#' @param err.max err.max
#' @param lb.opt Lower bound
#' @param ub.opt Upper bound
#' @return List describing optimized model
#' @export
#' @family comparison functions
betaTrim_mle <- function(
  x,
  m,
  p,
  pct0 = 0.25,
  init=c(1, 1, 0.5),
  iter.max = 200,
  err.max = 1e-6,
  lb.opt = c(0, 0, 0),
  ub.opt = c(log(1e4), log(1e4), 1)
) {
  if(is.null(init)) {
     m1p <- mean(p)
     m2p <- var(p)
     sc <- m1p * (1 - m1p) / m2p - 1
     init <- c(log(m1p * sc), log((1 - m1p) * sc), 0.9)
  }
  beta_bulk <- function(input) {
    rate1 <- exp(input[1])
    rate2 <- exp(input[2])
    pi0 <- input[3]
    if((rate1 > 0) & (rate2 > 0) & (pi0 > 0) & (pi0 < 1)) {
      n <- length(p)
      ub <- quantile(p, probs = 0.5 + pct0)
      lb <- quantile(p, probs = 0.5 - pct0)
      ix <- (p <= ub) & (p >= lb)
      dbt <- dbeta(p, shape1 = rate1, shape2 = rate2, log = TRUE)
      dbm <- log(
        1 - pi0 * (
          pbeta(ub, shape1 = rate1, shape2 = rate2, log = FALSE)
          - pbeta(lb, shape1 = rate1, shape2 = rate2, log = FALSE)
        )
      )
      betavec <- ix * (dbt + log(pi0)) + (1 - ix) * dbm
      betavec[is.na(betavec)] <- 0
      -sum(betavec)
    } else {
      Inf
    }
  }
  optout <- nloptr(
    init,
    beta_bulk,
    lb = lb.opt,
    ub = ub.opt,
    opts = list(
      algorithm = "NLOPT_GN_DIRECT_L",
      maxeval = iter.max,
      ftot_abs = err.max,
      ftol_rel = err.max
    )
  )
  coefout <-  optout[["solution"]]
  rate1 <- exp(coefout[1])
  rate2 <- exp(coefout[2])
  pi0 <- coefout[3]
  f0 <- dbetabinom.vec(x, m, rate1 = rate1, rate2 = rate2, log = FALSE)
  out <- list(coef = c(rate1, rate2, pi0), pi0 = pi0, f0 = f0, optout = optout)
}

#' @title Direct B-spline estimation.
#'
#' @details
#' Estimate the overall model directly via B-spline. Remove the basis that only
#' covers bin at the boundary.
#'
#' @param x x
#' @param m m
#' @param p p
#' @param breaks breaks
#' @param k Spline order
#' @param pi.init Initial weights
#' @param ncores Number of cores to use
#' @param err.max err.max
#' @param iter.max Max iterations
#' @return List describing optimized model
#' @export
#' @family comparison functions
emBspl <- function(
  x,
  m,
  p,
  breaks = seq(0, 1, length.out = 101),
  k = 4,
  pi.init = rep(1, length(breaks) + k - 4) / (length(breaks) + k - 4),
  ncores = 1,
  err.max = 1e-5,
  iter.max = 200
) {
    p[p == 0] <- 1 / max(m)^3
    p[p == 1] <- 1 - 1 / max(m)^3
  n <- length(p)
  bspl <- bsplfun( ## create the break points and coefficients of the B spline functions
    range(breaks),
    breaks = breaks,
    k = k,
    ncores = ncores
  )
  bspl <- bspl[-c(1, length(bspl))]
  nb <- length(bspl)
  dmtx <- evBspl(p, bspl, ncores = ncores) ### get the design matrix
  err <- 1e5 * err.max
  iter <- 0
  ll.init <- -Inf
  ll.all <- ll.init
  err.all <- err
  while ((is.na(err)|err > err.max) & iter < iter.max) {
    dtot <-  dmtx %*% matrix(pi.init, length(pi.init), 1) 
    dtot[dtot == 0] <- min(dtot[dtot > 0])
    post <- (
      dmtx
      *
      (matrix(pi.init, n, nb, byrow = TRUE))
      /
      matrix(dtot, n, nb, byrow = FALSE)
    )
    pi <- colSums(post)
    pi <- pi / sum(pi)
    ell <- apply(dmtx, 1, function(ii) sum(ii * pi))
    ll <- sum(log(ell))
    err <- ll - ll.init
    err <- max(err, err / max(abs(ll.init), abs(ll)))
    ll.all <- c(ll.all, ll)
    err.all <- c(err.all, err)
    pi.init <- pi
    ll.init <- ll
    iter <- iter + 1
  }
  binmtx <- getDesignMtx(bspl, x, m, ncores = ncores) ### get the design matrix
  list(
    pi = as.numeric(pi.init),
    post = post,
    bspl = bspl,
    dmtx = dmtx,
    binmtx = binmtx,
    f = as.numeric(binmtx %*% pi.init),
    ll.all = ll.all,
    err.all = err.all,
    convergence = list(
      err = err,
      err.max = err.max,
      converged = (err < err.max),
      niter = iter,
      ll = ll.init
    ),
    controls = list(k = k, nbasis = nb, breaks = breaks)
  )
}

#' @title Wrap alternative models.
#'
#' @details
#' wrapper of EBO and EBE
#' @param x x
#' @param m m
#' @param p p
#' @param breaks breaks
#' @param pi.init Initial weights
#' @param pct0 pct0
#' @param init init
#' @param iter.max Max iterations
#' @param err.max err.max
#' @param ncores Number of cores to use
#' @param lb.opt Lower bound
#' @param ub.opt Upper bound
#' @return List describing model
#' @export
#' @family comparison functions
ebBeta <- function(
  x,
  m,
  p,
  breaks = seq(0, 1, length.out = 101),
  k = 4,
  pi.init = rep(1, length(breaks) + k - 4) / (length(breaks) + k - 4),
  pct0 = 0.45,
  init = NULL,
  iter.max = 200,
  err.max = 1e-5,
  ncores = 1,
  lb.opt = c(0, 0, 0),
  ub.opt = c(log(1e4), log(1e4), 1)
) {
  mod <- emBspl(
    x,
    m,
    p,
    breaks = breaks,
    k = k,
    pi.init = pi.init,
    ncores = ncores,
    err.max = err.max,
    iter.max = iter.max
  )
  f <- mod[["f"]]  
  null <- betaTrim_mle(
    x,
    m,
    p,
    pct0,
    init = init,
    iter.max = iter.max,
    err.max = err.max,
    lb.opt,
    ub.opt
  )
  f0 <- null[["f0"]]
  pi0 <- null[["pi0"]]
  locfdr <- pi0 * f0 / f
  locfdrnorm <- locfdr / max(locfdr)
  coefnull <- null[["coef"]]
  out <- mod
  out[["coef.null"]] <- list(
    shape1 = coefnull[1],
    shape2 = coefnull[2],
    pi0 = pi0
  )
  out[["f0"]] <- f0
  out[["f"]] <- f
  out[["pi0"]] <- pi0
  out[["locfdr"]] <- locfdr
  out[["locfdrnorm"]] <- locfdrnorm
  out[["null"]] <- null[["optout"]]
  out
}