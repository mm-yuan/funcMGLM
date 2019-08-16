##########################################################################################
summary.mvglmer <- function (object, ...) {
  families <- object$families
  n_outcomes <- length(families)
  components <- object$components
  extract_components <- function (nam) {
    components[grep(nam, names(components), fixed = TRUE)]
  }
  respVars <- unlist(extract_components("respVar"), use.names = FALSE)
  descrpt <- data.frame(" " = unlist(extract_components("N"), use.names = FALSE),
                        row.names = respVars, check.rows = FALSE, check.names = FALSE)
  out <- list(n = components$n1, descrpt = descrpt, D = object$postMeans$D,
              families = families, respVars = respVars,
              control = object$control, mcmc.info = object$mcmc.info,
              DIC = object$DIC, pD = object$pD, call = object$call, engine = object$engine)
  for (i in seq_len(n_outcomes)) {
    out[[paste0("Outcome", i)]] <- data.frame("PostMean" = object$postMeans[[paste0("betas", i)]],
                                              "StDev" = object$StDev[[paste0("betas", i)]],
                                              "StErr"= object$StErr[[paste0("betas", i)]],
                                              "2.5%" = object$CIs[[paste0("betas", i)]][1, ],
                                              "97.5%" = object$CIs[[paste0("betas", i)]][2, ],
                                              "P" = object$Pvalues[[paste0("betas", i)]],
                                              "Rhat" = object$Rhat[[paste0("betas", i)]],
                                              row.names = names(object$postMeans[[paste0("betas", i)]]),
                                              check.names = FALSE)
    if (families[[i]][["family"]] == "gaussian") {
      D <- data.frame("PostMean" = object$postMeans[[paste0("sigma", i)]],
                      "StDev" = object$StDev[[paste0("sigma", i)]],
                      "StErr"= object$StErr[[paste0("sigma", i)]],
                      "2.5%" = object$CIs[[paste0("sigma", i)]][1],
                      "97.5%" = object$CIs[[paste0("sigma", i)]][2],
                      "P" = object$Pvalues[[paste0("sigma", i)]],
                      "Rhat" = object$Rhat[[paste0("sigma", i)]],
                      row.names = "sigma", check.names = FALSE)
      out[[paste0("Outcome", i)]] <- rbind(out[[paste0("Outcome", i)]], D)
    }
  }
  class(out) <- "summary.mvglmer"
  out
}

print.summary.mvglmer <- function (x, digits = max(4, getOption("digits") - 4), ...) {
  cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
  cat("Data Descriptives:")
  cat("\nNumber of Groups:", x$n)
  cat("\nNumber of Observations:\n")
  print(x$descrpt)
  cat("\n")
  if (!is.null(x$DIC)){
    model.sum <- data.frame(DIC = x$DIC, pD = x$pD, row.names = "")
    print(model.sum)
  }
  cat("\nRandom-effects covariance matrix:\n")
  D <- x$D
  ncz <- nrow(D)
  diag.D <- ncz != ncol(D)
  sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
  if (ncz > 1) {
    if (diag.D) {
      dat <- as.data.frame(round(rbind(sds), digits))
      names(dat) <- "StdDev"
    } else {
      corrs <- cov2cor(D)
      corrs[upper.tri(corrs, TRUE)] <- 0
      mat <- round(cbind(sds, corrs[, -ncz]), digits)
      mat <- rbind(mat)
      mat <- apply(mat, 2, sprintf, fmt = "% .4f")
      mat[mat == mat[1, 2]] <- ""
      mat[1, -1] <- abbreviate(colnames(D)[-ncz], 6)
      colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
      dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
      names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
      row.names(dat) <- c(dimnames(D)[[1]])
    }
  } else {
    dat <- data.frame("StdDev" = c(sds, x$sigma),
                      row.names = if (!is.null(x$sigma)) c(rownames(D), "Residual") else rownames(D),
                      check.rows = FALSE, check.names = FALSE)
  }
  print(dat)
  n_outcomes <- length(x$families)
  for (i in seq_len(n_outcomes)) {
    cat("\nOutcome:", x$respVars[i],"\n")
    print(round(x[[paste0("Outcome", i)]], digits))
  }
  cat("\nMCMC summary:\n")
  tt <- x$mcmc.info$elapsed.mins
  cat("engine:", x$engine, 
      "\niterations:", x$control$n.iter, 
      if (x$engine == "JAGS") 
        paste("\nadapt:", x$control$n.adapt,
              "\nburn-in:", x$control$n.burnin)
      else 
        paste("\nwarmup:", x$control$n.warmup), 
      "\nthinning:", x$control$n.thin,
      "\ntime:", if (tt > 60) round(tt/60, 1) else round(tt, 1),
      if (tt > 60) "hours" else "min")
  cat("\n\n")
  invisible(x)
}

plot.mvglmer <- function (x, which = c("trace", "autocorr", "density"),
                          param = c("betas", "sigma", "D"),
                          ask = TRUE, ...) {
  if (!inherits(x, "mvglmer"))
    stop("Use only with 'mvglmer' objects.\n")
  which <- match.arg(which)
  if (which %in% c("trace", "density", "autocorr")) {
    param <- match.arg(param, several.ok = TRUE)
    if (any(param == "D")) {
      keepD <- lower.tri(x$postMeans$D, TRUE)
      x$mcmc$D <- t(apply(x$mcmc$D, 1, c))[, c(keepD)]
      dnams <- which(keepD, arr.ind = TRUE)
      colnames(x$mcmc$D) <- paste0("D[", dnams[, 1], ", ", dnams[, 2], "]")
    }
    if (any(param == "tauBs")) {
      colnames(x$mcmc$tauBs) <- "tauBs"
    }
    which_parms <- unlist(sapply(param,
                                 function (pat) grep(paste0("^", pat), names(x$mcmc))),
                          use.names = FALSE)
    pp <- do.call(cbind, x$mcmc[which_parms])
    nams <- colnames(pp)
    op <- if (ask) par(mfrow = c(2, 2), ask = ask) else par(mfrow = c(4, 2))
    if (which == "trace") {
      for (i in 1:ncol(pp))
        plot(pp[, i], type = "l", xlab = "iterations", ylab = nams[i])
    } else if (which == "density") {
      for (i in 1:ncol(pp)) {
        bw <- bw.SJ(pp[, i]) * 1.5
        plot(density(pp[, i], bw = bw), xlab = nams[i],
             main = paste("Density of", nams[i]))
      }
    } else {
      for (i in 1:ncol(pp))
        acf(pp[, i], ylab = nams[i], main = paste("Series", nams[i]))
    }
    par(op)
  }
  invisible()
}

fixef.mvglmer <- function (object, ...) {
  if (!inherits(object, "mvglmer"))
    stop("Use only with 'mvglmer' objects.\n")
  comps <- object$components
  nams_outcomes <- unlist(comps[grep("respVar", names(comps), fixed = TRUE)], 
                          use.names = FALSE)
  pMeans <- object$postMeans
  betas <- pMeans[grep("betas", names(pMeans), fixed = TRUE)]
  names(betas) <- nams_outcomes
  betas
}

bind_chains <- function (ar) {
  d <- dim(ar)
  e <- seq_len(d[2L]) * d[1L]
  s <- c(1, head(e, -1) + 1)
  ind <- mapply(seq, from = s, to = e, SIMPLIFY = FALSE)
  m <- array(0.0, c(d[1L] * d[2L], d[3L]))
  for (i in seq_len(d[2L])) {
    m[ind[[i]], ] <- ar[, i, ]
  }
  colnames(m) <- dimnames(ar)[[3]]
  m
}

fix_D <- function (D) {
  d <- dim(D)
  k <- round(sqrt(d[2L]))
  m <- array(0.0, c(d[1L], k, k))
  for (i in seq_len(d[1L]))
    m[i, , ] <- matrix(D[i, ], k, k)
  m
}

fix_b <- function (b, n_RE) {
  d <- dim(b)
  n <- round(d[2L] / n_RE)
  m <- array(0.0, c(d[1L], n, n_RE))
  for (i in seq_len(d[1L]))
    m[i, , ] <- matrix(b[i, ], n, n_RE)
  m
}

stdErr <- function (x) {
  x <- as.matrix(x)
  vars <- apply(x, 2L, var)
  ess <- effectiveSize(x)
  sqrt(vars / ess)
}

effectiveSize <- function (x) {
  # copied (and made a bit more efficient) from the coda package
  spectrum0.ar <- function (x) {
    d <- dim(x)
    nrx <- d[1L]
    ncx <- d[2L]
    v0 <- numeric(ncx)
    res <- as.matrix(lm.fit(cbind(1, seq_len(nrx)), cbind(x, x))$residuals)
    for (i in seq_len(ncx)) {
      if (identical(all.equal(sd(res[, i]), 0), TRUE)) {
        v0[i] <- 0
      } else {
        ar.out <- ar(x[, i], aic = TRUE)
        v0[i] <- ar.out$var.pred / (1 - sum(ar.out$ar))^2
      }
    }
    v0
  }
  x <- as.matrix(x)
  spec <- spectrum0.ar(x)
  ifelse(spec == 0, 0, nrow(x) * apply(x, 2L, var) / spec)
}

computeP <- function (x) {
  above <- mean(x >= 0)
  below <- mean(x < 0)
  2 * min(above, below)
}

modes <- function (y) {
  test <- try(d <- density(y, bw = "nrd", adjust = 3, n = 1000), silent = TRUE)
  if (!inherits(test, "try-error")) d$x[which.max(d$y)] else NA
}

printCall <- function (call) {
  d <- deparse(call)
  if (length(d) <= 3) {
    paste(d, sep = "\n", collapse = "\n")
  } else {
    d <- d[1:3]
    d[3] <- paste0(d[3], "...")
    paste(d, sep = "\n", collapse = "\n")
  }
}
