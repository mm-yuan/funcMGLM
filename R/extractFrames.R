
# Adapted extractFrames function for HC(compute inverse of Z matrice)

#####################################################################################
extractFrames <- function (formula, data) {
  Terms <- terms(formula)
  term_labels <- attr(Terms, "term.labels")
  which_RE <- grep("|", term_labels, fixed = TRUE)
  namesVars <- all.vars(formula)
  respVar <- as.character(formula)[2L]
  # Fixed Effects
  formYx <- paste(term_labels[-which_RE], collapse = " + ")
  formYx <- as.formula(paste(respVar, "~", formYx))
  TermsX <- terms(formYx, data = data)
  mfX <- model.frame(TermsX, data)
  TermsX <- terms(mfX)
  X <- model.matrix(TermsX, data)
  # Random Effects
  spl <- unlist(strsplit(term_labels[which_RE], " | ", fixed = TRUE))
  idVar <- spl[2L]
  data <- data[complete.cases(data[namesVars]), ]
  id <- data[[idVar]]
  id <- match(id, unique(id))
  #**************
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  offset_sft = offset - 1
  ZrowsStart <- offset[1:(length(offset)-1)]
  ZrowsEnd = offset_sft[2:length(offset_sft)]
  Zrows = as.matrix(data.frame(ZrowsStart, ZrowsEnd))
  #**************
  formYz <- paste(spl[1], collapse = " + ")
  formYz <- as.formula(paste(respVar, "~", formYz))
  TermsZ <- terms(formYz, data = data)
  mfZ <- model.frame(TermsZ, data = data)
  TermsZ <- terms(mfZ)
  Z <- model.matrix(TermsZ, data)
  #**************
  Zt = t(Z)
  Zinv = ginv(Z)
  Ztinv = ginv(Zt)
  # individual Z inverse
  Zv = matrix(0, nrow(Zinv), ncol(Zinv))
  for (i in 1:length(ZrowsStart)) {
    index1 = ZrowsStart[i]
    index2 = ZrowsEnd[i]
    Z_sub = Z[index1:index2, ]
    Z_sub_inv = ginv(Z_sub)
    Zv[, index1:index2] = Z_sub_inv
  }
  #**************
  # response variable
  y <- model.response(mfX)
  if (is.factor(y))
    y <- as.vector(unclass(y) - 1)
  # hierarchical centering
  find_positions <- function (nams1, nams2) {
    nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
    vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
              glob2rx(paste0("*:", nams1)))
    out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
    out
  }
  check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  has_interceptX <- attr(TermsX, "intercept")
  has_interceptZ <- attr(TermsZ, "intercept")
  performHC <- has_interceptX && (has_interceptX == has_interceptZ)
  if (performHC) {
    terms.labs_X <- attr(TermsX, "term.labels")
    terms.labs_Z <- attr(TermsZ, "term.labels")
    # check for time-varying covariates
    timeTerms <- if (length(terms.labs_Z))
      unlist(lapply(terms.labs_Z, FUN = function(x) grep(x, colnames(X), fixed = TRUE)))
    which_td <- unname(which(apply(X, 2, check_td, id = id)))
    all_TDterms <- unique(c(timeTerms, which_td))
    baseline <- seq_len(ncol(X))[-all_TDterms]
    ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions, 
                                           nams2 = colnames(X)))
    ind_colmns2 <- seq_len(ncol(X))
    ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
    data.id <- data[!duplicated(id), ]
    Xhc <- if (length(terms.labs_Z)) {
      mfHC <- model.frame(TermsX, data = data.id)
      which.timevar <- unique(unlist(lapply(terms.labs_Z, 
                                            FUN = function (x) grep(x, names(mfHC), fixed = TRUE))))
      mfHC[which.timevar] <- lapply(mfHC[which.timevar], 
                                    function (x) { x[] <- 1; x })
      model.matrix(formYx, mfHC)
    } else {
      model.matrix(formYx, model.frame(TermsX, data = data.id))
    }
  }
  environment(TermsX) <- environment(TermsZ) <- NULL
  #***************************
  Xc = scale(X[, -1], center = TRUE, scale = FALSE) # except intercept
  Xc = cbind(X[,1], Xc)
  Xs = scale(X[, -1], center = TRUE, scale = TRUE) # except intercept
  Xs = cbind(X[,1], Xs)
  
  XhcC = scale(Xhc[, unlist(ind_colmns[1])[-1]], center = TRUE, scale = FALSE) # covariates except time 
  XhcC = cbind(Xhc[, -unlist(ind_colmns[1])[-1]], XhcC)
  XhcS = scale(Xhc[, unlist(ind_colmns[1])[-1]], center = TRUE, scale = TRUE) # covariates except time 
  XhcS = cbind(Xhc[, -unlist(ind_colmns[1])[-1]], XhcS)
  
  #if (length(terms.labs_Z) != 0) { 
  Zc = scale(Z[, unlist(ind_colmns[-1])], center = TRUE, scale = FALSE) # Zc/Zs only involves time 
  Zc = cbind(Z[,1], Zc)
  Zc_inv = ginv(Zc)
  
  Zs = scale(Z[, unlist(ind_colmns[-1])], center = TRUE, scale = TRUE) # Zc/Zs only involves time 
  Zs = cbind(Z[,1], Zs)
  Zs_inv = ginv(Zs)
  #} else{ Zc_inv = c(); Zs_inv = c()}
  #***************************
  if (ncol(X) >2) {
    means_X = apply(X[, -1], 2, mean)
    SDs_X = apply(X[, -1], 2, sd)
  } else{
    means_X = mean(X[, -1])
    SDs_X = sd(X[, -1])
  }
  mean_sd_X = means_X/SDs_X
  
  if (ncol(Z) >2) {
    means_Z = apply(Z[, -1], 2, mean)
    SDs_Z = apply(Z[, -1], 2, sd)
  } else{
    means_Z = mean(Z[, -1])
    SDs_Z = sd(Z[, -1])
  }
  mean_sd_Z = means_Z/SDs_Z
  
  #ind_excludeZ = sort(c(ind_colmns2, unlist(ind_colmns[1])[-1]))
  ind_excludeZ = unlist(ind_colmns[1])[-1]
  
  if (length(ind_excludeZ) > 1) {
    means_Xhc = apply(Xhc[, ind_excludeZ], 2, mean)
    SDs_Xhc = apply(Xhc[, ind_excludeZ], 2, sd) 
  } else if (length(ind_excludeZ) == 1) {
    means_Xhc = mean(Xhc[, ind_excludeZ])
    SDs_Xhc = sd(Xhc[, ind_excludeZ]) 
  }
  mean_sd_Xhc = means_Xhc/SDs_Xhc
  
  #***************************
  # extract results
  list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
       id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd, 
       y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,
       Z_ = Z, Zinv = Zinv, Ztinv = Ztinv, Zv = Zv, Zc = Zc,  Zs = Zs,
       means_X = means_X, SDs_X = SDs_X, mean_sd_X = mean_sd_X,
       means_Z = means_Z, SDs_Z = SDs_Z, mean_sd_Z = mean_sd_Z,
       means_Xhc = means_Xhc, SDs_Xhc = SDs_Xhc, mean_sd_Xhc = mean_sd_Xhc,
       TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
       Xhc = Xhc, colmns_HC = ind_colmns, colmns_nHC = ind_colmns2,
       ncx = ncol(X), ncz = ncol(Z))
}


