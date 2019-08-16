
mvglmer <- function (formulas, data, families, engine = c("JAGS", "STAN"), 
                     overdispersion = FALSE, priors = NULL, init = NULL, 
                     control = NULL, optionHC = c("NHC","HC"), scaling = c("Non","center","standardize")) {
  
########
# Data #
########
  
  cl <- match.call()
  #engine <- match.arg(engine)
  if (!is.list(families))
    stop("'families' must be a list of family objects.")
  # depending on the input of the user, set families to the corresponding
  # base R functions
  families[] <- lapply(families, function (family) {
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    family
  })
  # using the formulas and the data extract and constuct the objects required to
  # pass to JAGS
  components <- lapply(unname(formulas), extractFrames, data = data)
  # perform some checks
  if (!all((ns <- sapply(components, `[[`, 'n')) == components[[1L]][['n']])) {
    stop("it seems that the number of subjects differ between the responses. More ",
         "specifically, according to the data the number of subjects per response is ",
         paste(paste(sapply(components, `[[`, 'respVar'), ns, sep = " = "),
               collapse = ", "), ".\n")
  }
  components <- unlist(components, recursive = FALSE)
  n_outcomes <- length(formulas)
  names(components) <- paste0(names(components),
                              rep(seq_len(n_outcomes),
                                  each = length(components) / n_outcomes))
  colmns_HC <- components[grep("colmns_HC", names(components), fixed = TRUE)]
  colmns_nHC <- components[grep("colmns_nHC", names(components), fixed = TRUE)]
  ncol_X <- components[grep("ncx", names(components), fixed = TRUE)] 
  seq_outcomes <- seq_len(n_outcomes)
  nams_vars <- c("N", "id", "Z_", "Zinv", "Zv", "Ztinv", "X", "Xhc", "ncx", "y",
                 "ZrowsStart", "ZrowsEnd","Xc", "Xs", "Zc", "Zs", "XhcC", "XhcS",
                 "means_X", "SDs_X", "mean_sd_X", "means_Z", "SDs_Z", "mean_sd_Z",
                 "means_Xhc", "SDs_Xhc", "mean_sd_Xhc") ## add "X" for nHC ##
  vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
  if (any(ind_td <- sapply(colmns_nHC, length))) {
    vars <- c(vars, paste0("X", which(ind_td > 0)))
  }
  Data <- c(list(n = components$n1), components[vars])
  Data$n_RE <- sum(unlist(components[grep("ncz", names(components), fixed = TRUE)]))
  RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                    sq = components[grep("ncz", names(components), fixed = TRUE)],
                    incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                    SIMPLIFY = FALSE)
  names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
  Data <- c(Data, RE_inds, unlist(colmns_HC, recursive = FALSE), colmns_nHC)
  ###*************
  listN = unlist(components[grep("N", names(components), fixed = TRUE)])
  Zrow_ind =  mapply(function (sq, incr) seq_len(sq) + incr,
                     sq = components[grep("N", names(components), fixed = TRUE)],
                     incr = cumsum(c(0, head(listN, -1))),
                     SIMPLIFY = FALSE)
  totalM = sum(unlist(components[grep("N", names(components), fixed = TRUE)]))
  ZZmatrix = matrix(0, nrow = totalM, ncol = Data$n_RE)
  ZZtinv = matrix(0, nrow = totalM, ncol = Data$n_RE)
  ZZinv = matrix(0, ncol = totalM, nrow = Data$n_RE)
  ZZc = matrix(0, nrow = totalM, ncol = Data$n_RE)  
  ZZs = matrix(0, nrow = totalM, ncol = Data$n_RE)  
  ZZct = matrix(0, ncol = totalM, nrow = Data$n_RE)  
  ZZst = matrix(0, ncol = totalM, nrow = Data$n_RE)  
  for (i in 1:length(listN)){
   ZZtinv[Zrow_ind[[i]], RE_inds[[i]]] = ginv(t(as.matrix(components[grep("Z_", names(components), fixed = TRUE)][[i]])))
   ZZinv[RE_inds[[i]], Zrow_ind[[i]]] = ginv(as.matrix(components[grep("Z_", names(components), fixed = TRUE)][[i]]))
   ZZc[Zrow_ind[[i]], RE_inds[[i]]] = as.matrix(components[grep("Zc", names(components), fixed = TRUE)][[i]])
   ZZs[Zrow_ind[[i]], RE_inds[[i]]] = as.matrix(components[grep("Zs", names(components), fixed = TRUE)][[i]])
   ZZct[RE_inds[[i]], Zrow_ind[[i]]] = t(as.matrix(components[grep("Zc", names(components), fixed = TRUE)][[i]]))
   ZZst[RE_inds[[i]], Zrow_ind[[i]]] = t(as.matrix(components[grep("Zs", names(components), fixed = TRUE)][[i]]))
  }
  Data <- c(Data, list(ZZinv = ZZinv, ZZtinv = ZZtinv, 
                       ZZc = ZZc, ZZs = ZZs, ZZct = ZZct, ZZst = ZZst))
  ###**************
  # control
  con <- list(n.processors = parallel::detectCores() - 1, n.chains = 2,
              working.directory = getwd(), clear.model = TRUE,
              seed = 1L, optimize_only = FALSE, verbose = FALSE)
  if (engine == "JAGS") {
    con$n.iter <- 28000L
    con$n.burnin <- 3000L
    con$n.thin <- 50L
    con$n.adapt <- 3000L 
  } else {
    con$n.iter <- 1000
    con$n.warmup <- floor(con$n.iter / 2)
    con$n.thin <- 1
    con$adapt_delta <- 0.8
  }
  #control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (!any(namc == "n.thin")) {
    con$n.thin <- if (engine == "JAGS") {
      max(1, floor((con$n.iter - con$n.burnin) * con$n.chains / 1000))
    } else {
      max(1, floor((con$n.iter - con$n.warmup) * con$n.chains / 1000))
    }
  }
  if (length(noNms <- namc[!namc %in% namC]) > 0)
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ######################################################################################
  # Priors
  if (engine == "JAGS") {
    prs <- list(priorR_D = diag(rep(as.numeric(NA), Data$n_RE), Data$n_RE),
                priorK_D = Data$n_RE + 1, A_RD = 0.5, B_RD = 0.01,
                tau_half_cauchy = 0.1)
    pr_taus_betas <- rep(list(0.01), n_outcomes)
    names(pr_taus_betas) <- paste0("tau_betas", seq_len(n_outcomes))
    prs <- c(prs, pr_taus_betas)
    if (any(sapply(families, function (x) x$family == "gaussian")) || overdispersion) {
      prs$A_tau <- 0.01
      prs$B_tau <- 0.01
    }
  } else {
    prs <- list(scale_sigmas = 5, scale_diag_D = 3, lkj_shape = 2,
                priorK_D = Data$n_RE + 1)
    pr_scale_betas <- rep(list(10), n_outcomes)
    names(pr_scale_betas) <- paste0("scale_betas", seq_len(n_outcomes))
    prs <- c(prs, pr_scale_betas)
  }
  if (!is.null(priors)) {
    lngths <- lapply(prs[(nam.prs <- names(priors))], length)
    if (!is.list(priors) || !isTRUE(all.equal(lngths, lapply(priors, length)))) {
      warning("'priors' is not a list with elements numeric vectors of appropriate ",
              "length; default priors are used instead.\n")
    } else {
      prs[nam.prs] <- priors
    }
  }
  Data <- c(Data, prs)
######################################################################################

###############
#  STAN.code  #
###############

myt <- function(times = 1) {
  tb <- "    "
  paste(rep(tb, times), collapse = "")
}

data_part <- function (families, ncolsXs, ncolsZs, extraXs, n_RE, colmns_HC, colmns_nHC, optionHC, scaling, totalM) {
  outcomes <- seq_along(families)
  data_outcome <- function (outcome, family, ncolsZ, extraX = FALSE, colmns_HC, colmns_nHC, optionHC, scaling) {
    type_long_outcome <- switch(family$family,
                                "gaussian" = paste0("vector[N", outcome, "] y", outcome),
                                "binomial" = paste0("int<lower=0, upper=1> y", outcome, "[N", outcome, "]"),
                                "poisson" = paste0("int<lower=0> y", outcome, "[N", outcome, "]"))
    if (optionHC == "HC") {
      paste0(myt(), "int N", outcome, ";\n",
             myt(), "int ncx", outcome, ";\n",
             myt(), "int id", outcome, "[N", outcome, "];\n",
             myt(), "int RE_ind", outcome, if (ncolsZ > 1) paste0("[", ncolsZ, "]"), ";\n",
             myt(), type_long_outcome, ";\n",
             myt(), "matrix[N", outcome, ", ", ncolsZ, "] ", "Z_", outcome, ";\n",
             myt(), "matrix[n, ncx", outcome, "] ", "Xhc", outcome, ";\n",
             if (scaling == "Non" & extraX) paste0(myt(), "matrix[N", outcome, ", ", "ncx", outcome, "] ",
                                                  "X", outcome, ";\n"),
             if (scaling != "Non" & length(colmns_HC) > 1) paste0(myt(), "int ZrowsStart", outcome, "[n];\n",
                                         myt(), "int ZrowsEnd", outcome, "[n];\n",
                                         myt(), "matrix[", ncolsZ,", N", outcome, "] ", "Zv", outcome, ";\n"))
        }
    else {
      paste0(myt(), "int N", outcome, ";\n",
             myt(), "int ncx", outcome, ";\n",
             myt(), "int id", outcome, "[N", outcome, "];\n",
             myt(), "int RE_ind", outcome, if (ncolsZ > 1) paste0("[", ncolsZ, "]"), ";\n",
             myt(), type_long_outcome, ";\n",
             myt(), "matrix[N", outcome, ", ", ncolsZ, "] ", "Z_", outcome, ";\n",
             myt(), "matrix[N", outcome, ", ncx", outcome, "] ", "X", outcome, ";\n") 
    }
  }
  def_HC_center <- function (colmns_HC_i, outcome, ncolsZ, ncolsX) {
    ncol_xhc = length(colmns_HC_i[[1]]) - 1 
    paste0(
      if (ncolsX > 2) paste0( myt(), "matrix[n, ", ncolsX, "] ", "XhcC", outcome, ";\n")
      else  paste0( myt(), "vector[n] XhcC", outcome, ";\n"),
      if (ncol_xhc > 1) 
        paste0( myt(), "vector[", ncol_xhc, "] means_Xhc",outcome, ";\n")
      else
        paste0(myt(), "real means_Xhc", outcome, ";\n"),
        paste0( myt(), "matrix[N", outcome,",", ncolsZ, "] ", "Zc", outcome, ";\n",
                if (ncolsZ > 2) 
                        paste0(myt(), "vector[", (ncolsZ -1), "] means_Z",outcome, ";\n")     
                       else if (ncolsZ == 2)
                        paste0( myt(), "real means_Z", outcome, ";\n")))}
  def_HC_std <- function (colmns_HC_i, outcome, ncolsZ, ncolsX) {
    ncol_xhc = length(colmns_HC_i[[1]]) - 1 
    paste0(
      if (ncolsX > 2) paste0( myt(), "matrix[n, ", ncolsX, "] ", "XhcS", outcome, ";\n")
      else  paste0( myt(), "vector[n] XhcS", outcome, ";\n"),
      if (ncol_xhc > 1) 
        paste0( myt(), "vector[", ncol_xhc, "] SDs_Xhc",outcome, ";\n",
                myt(), "vector[", ncol_xhc, "] mean_sd_Xhc",outcome, ";\n")
      else
        paste0( myt(), "real SDs_Xhc", outcome, ";\n",
                myt(), "real mean_sd_Xhc", outcome, ";\n"),
        paste0( myt(), "matrix[N", outcome,",", ncolsZ, "] ", "Zs", outcome, ";\n",
                if (ncolsZ > 2)  
                  paste0(myt(), "vector[", (ncolsZ -1), "] SDs_Z",outcome, ";\n",
                         myt(), "vector[", (ncolsZ -1), "] mean_sd_Z",outcome, ";\n")
                else if (ncolsZ == 2)
                  paste0(myt(), "real SDs_Z", outcome, ";\n",
                         myt(), "real mean_sd_Z", outcome, ";\n")))}
  def_NHC_center <- function (outcome, ncolsX) {
    paste0(if (ncolsX > 1) 
             paste0( myt(), "matrix[N", outcome, ", ", ncolsX, "] ", "Xc", outcome, ";\n",
                     myt(), "vector[", (ncolsX -1), "] means_X",outcome, ";\n")
           else
             paste0( myt(), "vector[N", outcome,  "] Xc", outcome, ";\n",
                     myt(), "real means_X", outcome, ";\n"))
  }
  def_NHC_std <- function (outcome, ncolsX) {
    paste0(if (ncolsX > 1) paste0( myt(), "matrix[N", outcome, ", ",ncolsX, "] ", "Xs", outcome, ";\n",
                                   myt(), "vector[", (ncolsX -1), "] SDs_X",outcome, ";\n",
                                   myt(), "vector[", (ncolsX -1), "] mean_sd_X",outcome, ";\n")
           else paste0(myt(), "vector[N", outcome,  "] Xs", outcome, ";\n",
                       myt(), "real SDs_X",outcome, ";\n",
                       myt(), "real mean_sd_X",outcome, ";\n"))
  } 
  df_extraX <- function (outcome, ncolsX, scaling) {
     paste0(if (scaling == "center") def_NHC_center(outcome, ncolsX),
            if (scaling == "standardize") def_NHC_std(outcome, ncolsX))
  }
  df_trf <- function (outcome, family, ncolsX, ncolsZ, extraX, colmns_HC, optionHC, scaling) {
    if (optionHC == "HC") {
      paste0(if (scaling == "center") def_HC_center(colmns_HC, outcome, ncolsZ, ncolsX),
             if (scaling == "standardize") def_HC_std(colmns_HC, outcome, ncolsZ, ncolsX),
             if (extraX)  df_extraX(outcome, ncolsX, scaling))
    }
    else { 
      paste0(if (scaling == "center") def_NHC_center(outcome, ncolsX),
             if (scaling == "standardize") def_NHC_std(outcome, ncolsX))
    }
  }
  def_inv <- function (outcomes, ncolsZs){
    if ( length(outcomes) > 1 )
    paste0(myt(), "matrix[", sum(unlist(ncolsZs)), ", ",totalM , "] ZZinv", ";\n",
           myt(), "matrix[", totalM, ", ", sum(unlist(ncolsZs)), "] ZZtinv", ";\n",
           myt(), "matrix[", totalM, ", ", sum(unlist(ncolsZs)),"] ZZ",  
           if (scaling == "center") paste0("c") else paste0("s"), ";\n",
           myt(), "matrix[", sum(unlist(ncolsZs)), ", ", totalM,"] ZZ",  
           if (scaling == "center") paste0("c") else paste0("s"), "t;\n")
    else
      paste0(myt(), "matrix[N1,", ncolsZs, "] Ztinv1", ";\n",
             myt(), "matrix[", ncolsZs, ",N1] Zinv1",  ";\n")
  }
  paste0("data {\n", myt(), "int n;\n", myt(), "int n_RE;\n", 
         paste0(mapply(data_outcome, outcomes, families, ncolsZs, extraXs, colmns_HC, colmns_nHC, optionHC, scaling), 
                collapse = ""),
         paste0(sapply(outcomes, function (outcome) 
           paste0(myt(), "real<lower=0> scale_betas", outcome, ";\n")), collapse = ""),
         if (any(sapply(families, `[[`, 'family') == "gaussian")) 
           paste0(myt(), "real<lower=0> scale_sigmas;\n"),
         myt(), "real<lower=0> scale_diag_D;\n",
         if (n_RE > 1)
           paste0(myt(), "real<lower=0> lkj_shape;\n"),
         if (scaling != "Non")paste0(mapply(df_trf, outcomes, families, ncolsXs, ncolsZs, extraXs, colmns_HC, optionHC, scaling), collapse = ""),
         if (optionHC == "HC" & scaling != "Non") def_inv(outcomes, ncolsZs),
         "}\n")
}

##########################################################################################

parameters <- function (families, ncolsXs, n_RE, scaling) {
  outcomes <- seq_along(families)
  set_parms <- function (outcome, ncolsX, family, scaling) {
    set_betas <- function (outcome, ncolsX, family, scaling) {
    if (scaling == "center"){if (ncolsX > 2) 
      paste0(myt(), "vector[", ncolsX, "] temp_betas", outcome, ";\n") else
        paste0(myt(), "real temp_betas", outcome, ";\n")
      }
    else if (scaling == "standardize") { if (ncolsX > 2) 
      paste0(myt(), "vector[", ncolsX, "] temp_betas", outcome, ";\n") else
        paste0(myt(), "real temp_betas", outcome, ";\n")
    }
    else{ paste0(myt(), "vector[ncx", outcome, "] betas", outcome, ";\n")
    }
    }
    paste0( 
    paste0(set_betas(outcome, ncolsX, family, scaling)),
    if (any(family$family == "gaussian")) 
       paste0(myt(), "real<lower = 0> sigma", outcome, ";\n"))
  }
  paste0("\nparameters {\n",
         paste0(mapply(set_parms, outcomes, ncolsXs, families, scaling), collapse = ""),
         if (optionHC == "HC")
           paste0(myt(), "matrix[n, n_RE] u;\n")
         else
           paste0(myt(), "matrix[n, n_RE] b;\n"),
         if (n_RE > 1)
           paste0(myt(), "vector<lower = 0>[n_RE] L_var_D;\n",
                  myt(), "cholesky_factor_corr[n_RE] L_corr_D;\n")
         else 
           paste0(myt(), "real<lower = 0> D;\n"),
         "}\n")
}

##########################################################################################

transformed_parameters <- function (families, ncolsZs, colmns_HC, colmns_nHC, RE_inds, optionHC, scaling) {
  outcomes <- seq_along(families)
  def_etas <- function (outcome) {
    paste0(myt(), "vector[N", outcome, "] eta", outcome, ";\n")
  }
  def_etas_int <- function (outcome) {
    paste0(myt(), "eta", outcome, " = rep_vector(0,N", outcome, ")",";\n")
  }
  colmns_HC2 <- unlist(colmns_HC, recursive = FALSE)
  nams <- names(colmns_HC2)
  ncols <- sapply(colmns_HC2, length)
  HC_part <- function (outcome, columns, nam, ncol, i, scaling) {
    if (ncol) {
      paste0(myt(2), "mu_u[i, ", i,"] = ", 
             if (scaling == "center") paste0("XhcC", outcome, "[i, ", columns, "] * temp_betas",outcome, "[", columns, "]", collapse = " + ") else
             if (scaling == "standardize") paste0("XhcS", outcome, "[i, ", columns, "] * temp_betas",outcome, "[", columns, "]", collapse = " + ") else
               paste0("Xhc", outcome, "[i, ", columns, "] * betas",outcome, "[", columns, "]", collapse = " + "), ";\n")
    } else {
      paste0(myt(2), "mu_u[i, ", i,"] = 0.0;\n")
    }
  }
  linpred_part <- function (outcome, ncolsZ, colmns_HC, colmns_nHC, RE_ind, scaling) {
    j_i <- paste0("j", outcome)
    nZ = seq_len(ncolsZ)
    paste0(myt(), "for (", j_i, " in 1:N", outcome, ") {\n",
           myt(2), "eta", outcome, "[", j_i, "] = ", 
           paste0(if (scaling == "Non") paste0("Z_"), 
                  if (scaling == "center") paste0("Zc"),
                  if (scaling == "standardize") paste0("Zs"),
                  outcome, "[", j_i, ", ", nZ, "] * u[id", outcome, "[", j_i, "], ",
                  RE_ind, "]", collapse = " + "),
           if (length(clm <- colmns_nHC)) {
             paste0("\n", myt(4), " + ", 
                    paste0(if (scaling == "Non") paste0("X"),
                           if (scaling == "center") paste0("Xc"),
                           if (scaling == "standardize") paste0("Xs"),
                           outcome, "[", j_i, ", ", clm, "] * ",
                           if (scaling == "Non") paste0("betas") else paste0("temp_betas"), 
                           outcome, "[", clm, "]", collapse = " + "), ";\n")
           } else ";\n",
           myt(), "}\n")
  }
  ########################
  ### functions for nHC
  def_yhats_part1 <- function (outcome, scaling) {
    paste0(if (scaling == "center") paste0(myt(), "vector[N", outcome, "] y", outcome, 
                                           "_hat = Xc",outcome, " * temp_betas",outcome, ";\n") else
           if (scaling == "standardize") paste0(myt(), "vector[N", outcome, "] y", outcome, 
                                                "_hat = Xs",outcome, " * temp_betas",outcome, ";\n") else
           paste0(myt(), "vector[N", outcome, "] y", outcome, "_hat = X",outcome, " * betas",outcome,";\n"))
  }
  def_yhats_part2 <- function (outcome, ncolsZ, colmns_HC, colmns_nHC, RE_ind) {
    j_i <- paste0("j", outcome)
    q_i <- paste0("q", outcome)
    paste0(myt(), "for (", j_i, " in 1:N", outcome, ") {\n",
           if (length(colmns_HC) > 1)
             paste0(myt(2), "for (", q_i, " in 1:", ncolsZ , ") {\n",
                    myt(3), "y", outcome, "_hat[", j_i, "] += ", 
                    "Z_", outcome, "[", j_i, ",", q_i,
                    "] * b[id", outcome, "[", j_i, "], RE_ind", outcome, "[",q_i,"]];\n",
                    myt(2), "}\n")
           else
             paste0(myt(3), "y", outcome, "_hat[", j_i, "] += ", 
                    "Z_", outcome, "[", j_i, ", 1] * b[id", outcome, "[", j_i, "], ",
                    RE_ind, "];\n"),
           myt(), "}\n")
  }
  ########################
  if (optionHC == "HC") {
    paste0("\ntransformed parameters {\n", paste0(mapply(def_etas, outcomes), collapse = ""), 
           myt(), "matrix[n, n_RE] mu_u;\n",
           myt(), "for (i in 1:n) {\n",
           paste0(mapply(HC_part, rep(outcomes, sapply(colmns_HC, length)), 
                         colmns_HC2, nams, ncols, seq_along(colmns_HC2), scaling), collapse = ""),
           myt(), "}\n",
           paste0(mapply(linpred_part, outcomes, ncolsZs, colmns_HC, colmns_nHC, RE_inds, scaling), collapse = ""),
           "}\n")
  }
  else{
  paste0("\ntransformed parameters {\n", paste0(mapply(def_yhats_part1, outcomes, scaling), collapse = ""), 
         paste0(mapply(def_yhats_part2, outcomes, ncolsZs, colmns_HC, colmns_nHC, RE_inds), collapse = ""),
         "}\n") 
  }
}

##########################################################################################

model <- function (families, ncolsXs, n_RE, scaling) {
  outcomes <- seq_along(families)
  RE_part <- paste0("\nmodel {\n",
                    if (n_RE > 1)
                      paste0(myt(),  "matrix[n_RE, n_RE] L_D;\n",
                             if (optionHC == "NHC") paste0(myt(),  "vector[n_RE] mu;\n"),
                             myt(),  "L_D = diag_pre_multiply(L_var_D, L_corr_D);\n",
                             myt(),  "L_var_D ~ cauchy(0, scale_diag_D);\n",
                             myt(),  "L_corr_D ~ lkj_corr_cholesky(lkj_shape);\n",
                             if (optionHC == "NHC") paste0(myt(),  "mu = rep_vector(0, n_RE);\n"),
                             myt(),  "for (i in 1:n) {\n",
                             if (optionHC == "HC")
                               paste0(myt(2), "u[i, ] ~ multi_normal_cholesky(mu_u[i, ], L_D);\n")
                             else 
                               paste0(myt(2), "b[i, ] ~ multi_normal_cholesky(mu, L_D);\n"))
                    else
                      paste0(myt(), "D ~ cauchy(0, scale_diag_D);\n",
                             myt(), "for (i in 1:n) {\n",
                             if (optionHC == "HC")
                               paste0(myt(2), "u[i, ] ~ normal(mu_u[i, ], D);\n")
                             else
                               paste0(myt(2), "b[i, ] ~ normal(mu, D);\n")),
                      myt(), "}\n")
  priors_part <- function (outcome, ncolsX, family, scaling) {
    if (ncolsX > 2) 
    paste0(if (scaling == "center")   paste0(myt(), "for (k", outcome, " in 2:", ncolsX , ") {\n"),
           if (scaling == "standardize")   paste0(myt(), "for (k", outcome, " in 2:", ncolsX, ") {\n"),
           if (scaling == "Non")   paste0(myt(), "for (k", outcome, " in 1:ncx", outcome, ") {\n"),
           if (scaling != "Non")  
             paste0( myt(2), "temp_betas") 
           else
             paste0( myt(2), "betas"), outcome, "[k", outcome, "] ~ normal(0.0, scale_betas", outcome, ");\n",
           myt(), "}\n")
    else
    paste0( if (scaling != "Non")  paste0( myt(), "temp_betas") 
            else paste0( myt(), "betas"), 
            outcome, " ~ normal(0.0, scale_betas", outcome, ");\n")
  }
  dist_part <- function (outcome, family, optionHC) {
    if (optionHC == "HC"){
      if (family$family == "gaussian") {
        paste0(myt(), "y", outcome, " ~ normal(eta", outcome, ", sigma", outcome, ");\n")
      } else if (family$family == "binomial") {
        switch (family$link,
                "logit" = paste0(myt(), "y", outcome, " ~ bernoulli_logit(eta", outcome, ");\n"),
                "probit" = paste0(myt(), "y", outcome, " ~ bernoulli(Phi_approx(eta", outcome, "));\n"),
                "cloglog" = paste0(myt(), "y", outcome, " ~ bernoulli(inv_cloglog(eta", outcome, "));\n")
        )
      } else if (family$family == "poisson") {
        paste0(myt(), "y", outcome, " ~ poisson_log(eta", outcome, ");\n")
      }
    }
    else{
      if (family$family == "gaussian") {
          paste0(myt(), "sigma", outcome, " ~ cauchy(0, scale_sigmas);\n",
            myt(), "y", outcome, " ~ normal(y", outcome, "_hat",", sigma", outcome, ");\n")
      } else if (family$family == "binomial") {
        switch (family$link,
                "logit" = paste0(myt(), "y", outcome, " ~ bernoulli_logit(y", outcome,"_hat", ");\n"),
                "probit" = paste0(myt(), "y", outcome, " ~ bernoulli(Phi_approx(y", outcome,"_hat", "));\n"),
                "cloglog" = paste0(myt(), "y", outcome, " ~ bernoulli(inv_cloglog(y", outcome,"_hat", "));\n")
        )
      } else if (family$family == "poisson") {
        paste0(myt(), "y", outcome, " ~ poisson_log(y", outcome, "_hat", ");\n")
      }
    }
  }
  paste0(RE_part, paste0(mapply(priors_part, outcomes, ncolsXs, families, scaling), collapse = ""),
         paste0(mapply(dist_part, outcomes, families, optionHC), collapse = ""), "}\n")
}

##########################################################################################

generated_quantities <- function (families, extraXs, colmns_HC, colmns_nHC, ncolsXs, ncolsZs, n_RE, optionHC, scaling) {
    outcomes <- seq_along(families)
    transform_int_nhc <- function (outcome, ncolsX, scaling) {
      if (ncolsX > 2) {
        paste0(myt(), "real Intercept", outcome," = temp_betas", outcome, 
               "[1] - dot_product(", 
               if (scaling == "center") paste0("means_X", outcome, ", temp_betas", outcome, "[2:", ncolsX, "]);\n"),
               if (scaling == "standardize") paste0("mean_sd_X", outcome, ", temp_betas", outcome, "[2:", ncolsX, "]);\n",
                                                    myt(), "vector[", (ncolsX -1),"] betas", outcome,"_part;\n"))
      }
      else{
        paste0(myt(), "real Intercept", outcome," = temp_betas", outcome, "[1] - ", 
               if (scaling == "center") paste0("means_X", outcome, " * temp_betas", outcome, ";\n"), 
               if (scaling == "standardize") paste0("mean_sd_X", outcome, " * temp_betas", outcome, ";\n",
                                                    myt(), "real betas", outcome,"_part;\n"))
      }
    }
    transform_int_hc <- function (outcome, extraX, colmns_HC_i, colmns_nHC_i, ncolsX, ncolsZ, scaling) {
      k = colmns_HC_i[[1]][-1]
      m = seq_len(ncolsZ)[-length(seq_len(ncolsZ))]
      xhc_part <- function (outcome, k, scaling) {
        i = seq_len(length(k))
        if (length(k) > 1) paste0(" - ", if (scaling == "center") paste0("means_Xhc") else paste0("mean_sd_Xhc"), 
                                  outcome, "[", i,"]"," * temp_betas", outcome, "[", k, "]",collapse = "")
        else paste0(" - ", if (scaling == "center") paste0("means_Xhc") else paste0("mean_sd_Xhc"), 
                    outcome, " * temp_betas", outcome, "[", k, "]", collapse = "")
      }
      xhc_extra <- function (outcome, colmns_nHC_i) {
        paste0(" - ", if (scaling == "center") paste0("means_X") else paste0("mean_sd_X"), 
               outcome, "[", (colmns_nHC_i -1),"]"," * temp_betas", outcome, "[", colmns_nHC_i, "]",collapse = "")
      }
      z_part <- function (outcome, j, m, scaling) {
        if (length(m) > 1) paste0(" - ",if (scaling == "center") paste0("means_Z") else paste0("mean_sd_Z"), 
                                  outcome, "[", m, "] * temp_betas", outcome, "[", j, "]",collapse = "")
        else paste0(" - ", if (scaling == "center") paste0("means_Z") else paste0("mean_sd_Z"), 
                    outcome, " * temp_betas", outcome,"[", j, "]", collapse = "")
      }
      paste0(myt(), "real Intercept", outcome," = temp_betas", outcome, "[1]",
             xhc_part(outcome, k, scaling),
             if (length(colmns_HC_i) > 1) z_part(outcome, unlist(colmns_HC_i[[-1]]), m, scaling),
             if (extraX) xhc_extra(outcome, colmns_nHC_i),";\n",
             if (ncolsX > 1 & scaling == "standardize") paste0(myt(), "vector[", (ncolsX - 1),"] betas", outcome,"_part;\n") )
    }
    def_final_betas <- function (outcome) {
      paste0(myt(), "vector[ncx", outcome, "] betas", outcome, ";\n")
    }
    transform_betas_nhc <- function (outcome, ncolsX) {
      if (ncolsX > 2) 
        paste0(myt(), "for (i in 1:", (ncolsX-1), ") {\n",
               myt(2),"betas", outcome, "_part[i] = temp_betas", outcome, "[i+1] / SDs_X", outcome, "[i];\n", myt(), "}\n")
      else paste0(myt(),"betas", outcome, "_part = temp_betas[2]", outcome, " / SDs_X", outcome, ";\n")
    }
    write_z <- function (outcome,colmns_HC_i,r) {
      if (length(r) > 1) {
        out_r <- vector("character", length(r))
        for (i in seq_along(r)) {
          ind <- r[i]
          out_r[i] <- paste0(myt(),"betas", outcome, "_part[", (ind-1), "] = temp_betas", outcome, "[", ind, "] / SDs_Z", outcome,"[", i, "]", ";\n")
        }
        paste(out_r, collapse = "")
      } else {
        paste0(myt(),"betas", outcome, "_part[", (r-1), "]", " = temp_betas", outcome, "[", r, "] / SDs_Z", outcome, ";\n")
      }
    }
    write_x <- function(outcome, k, sd, ind_extraX){
      if (length(k) > 1) {
        out_k <- vector("character", length(k))
        for (i in seq_along(k)) {
          ind <- k[i]
          out_k[i] <- paste0(myt(),"betas", outcome, "_part[", (ind-1), "] = temp_betas", outcome, "[", ind, "] / ", sd, outcome,"[", i, "]", ";\n")
        }
        paste(out_k, collapse = "")
      } else{
        paste0(myt(),"betas", outcome, "_part[", (k-1), "]", " = temp_betas", outcome, "[", k, "] / ", sd, outcome, ind_extraX, ";\n")
      }
    }
    transform_betas_hc <- function (outcome,colmns_HC_i,colmns_nHC_i,extraX) {
      k1 = unlist(colmns_HC_i[1])[-1]
      k2 = colmns_nHC_i 
      paste0(write_x(outcome, k1, "SDs_Xhc", ""),
      if (extraX) write_x(outcome, k2, "SDs_X", paste0("[",(k2-1),"]")),
      if (length(colmns_HC_i) > 1) write_z(outcome,colmns_HC_i, unlist(colmns_HC_i[[-1]])))
    }
    merge_betas <- function (outcome, ncolsX, scaling) {
      if (ncolsX > 2) {
        paste0(myt(), "betas", outcome, "= append_row(Intercept", outcome, 
               if (scaling == "standardize") paste0(", betas", outcome, "_part);\n")
               else paste0(", temp_betas", outcome, "[2: ", ncolsX, "]);\n"))
      }
      else{
        paste0(myt(), "betas", outcome,"= [Intercept", outcome, 
               if (scaling == "standardize") paste0(", betas", outcome, "_part]' ;\n")
               else paste0(", temp_betas", outcome, "[2]);\n"))
      }
    }
    transform_D <- function (outcomes, scaling) {
      if ( length(outcomes) == 1 ) {
        paste0("D =  Zinv1 * Z", 
               if (scaling == "standardize") paste0("s1 * tmp_D * (Zs1)' * Ztinv1;\n")
               else paste0("c1 * tmp_D * (Zc1)' * Ztinv1;\n") )
      }else{ paste0("D =  ZZinv * ZZ", 
                    if (scaling == "standardize") paste0("s * tmp_D * ZZst * ZZtinv;\n")
                    else paste0("c * tmp_D * ZZct * ZZtinv;\n") )
      }
    }
    transform_b <- function (outcome, n_RE, scaling, colmns_HC_i) {
      if (length(colmns_HC_i) > 1){
        paste0(
        if (n_RE > 1) 
        paste0(myt(2), "b[i, RE_ind", outcome, "] = (Zv", outcome, "[, ZrowsStart", outcome, 
               "[i]:ZrowsEnd", outcome, "[i]] * Z", 
               if (scaling == "standardize") paste0("s") else paste0("c"), outcome, "[ZrowsStart", outcome,
               "[i]:ZrowsEnd", outcome, "[i], ] * tmp_b[i, RE_ind", outcome, "]')';\n")
        else
        paste0(myt(2), "b[i] = (Zv1[, ZrowsStart1[i]:ZrowsEnd1[i]] * Zs1[ZrowsStart1[i]:ZrowsEnd1[i], ] * tmp_b[i])';\n"))
      } else{
        paste0(myt(2), "b[i, RE_ind", outcome, "] = tmp_b[i, RE_ind", outcome, "];\n")
      }
    }
    paste0("\ngenerated quantities {\n",
           if (scaling != "Non")  paste0(paste0(mapply(def_final_betas, outcomes), collapse = ""),
                                         if (optionHC == "NHC") paste0(mapply(transform_int_nhc, outcomes, ncolsXs, scaling), collapse = ""),
                                         if (optionHC == "HC") paste0(mapply(transform_int_hc, outcomes, extraXs, colmns_HC, colmns_nHC, ncolsXs, ncolsZs, scaling), collapse = "")),
           if (n_RE > 1) paste0(myt(), "matrix[n_RE, n_RE] D;\n",
                                if (optionHC == "HC" & scaling != "Non") paste0(myt(), "matrix[n_RE, n_RE] tmp_D;\n")),
           if (optionHC == "HC")  paste0(myt(), "matrix[n, n_RE] b;\n", 
                                         if (scaling != "Non") paste0(myt(),"matrix[n, n_RE] tmp_b;\n")),
           if (scaling == "standardize"& optionHC == "HC") 
              paste0(mapply(transform_betas_hc, outcomes, colmns_HC,colmns_nHC, extraXs), collapse = ""),
           if (scaling == "standardize" & optionHC == "NHC") 
              paste0(mapply(transform_betas_nhc, outcomes, ncolsXs), collapse = ""),
           if (scaling != "Non") paste0(mapply(merge_betas, outcomes, ncolsXs, scaling), collapse = ""), 
           if (n_RE > 1) paste0(myt(),if (optionHC == "HC" & scaling != "Non") paste0("tmp_D") else paste0("D"),
                                " = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)';\n"),
           if (optionHC == "HC")  paste0(myt(), 
                                         if (scaling != "Non") paste0("tmp_"), "b = u - mu_u;\n",
                                         if (scaling != "Non") paste0(myt(), transform_D(outcomes, scaling),
                                                                      myt(), "for (i in 1:n) {\n",
                                                                      paste0(mapply(transform_b, outcomes, n_RE, scaling, colmns_HC), collapse = ""), myt(), "}\n")),
           "}\n")
  }
##########################################################################################

###################
# Write STAN code #
###################

model_name <- paste0("mvglmer", sample(1e06, 1), if (engine == "JAGS") ".txt" else ".stan")
if (engine == "JAGS") {
  cat(build_model(families, seq_along(families), colmns_HC, colmns_nHC, overdispersion,
                  Data$n_RE), file = file.path(con$working.directory, model_name))
} else {
  cat(data_part(families, ncol_X, lapply(colmns_HC, length), lapply(colmns_nHC, length), 
                Data$n_RE, colmns_HC, colmns_nHC, optionHC, scaling, totalM),
      parameters(families, ncol_X, Data$n_RE, scaling),
      transformed_parameters(families, lapply(colmns_HC, length), colmns_HC, colmns_nHC, RE_inds, optionHC, scaling),
      model(families, ncol_X, Data$n_RE, scaling), 
      generated_quantities(families, lapply(colmns_nHC, length), colmns_HC, colmns_nHC, ncol_X, lapply(colmns_HC, length), Data$n_RE, optionHC, scaling),
      file = file.path(con$working.directory, model_name))
}

##########################################################################################

########
# fit  #
########

# parameters to save

params <- paste0('betas', seq_len(n_outcomes))
if (any(ind_gs <- sapply(families, function (x) x$family == "gaussian"))) {
  params <- c(params, paste0("sigma", which(ind_gs)))
}
params <- if (engine == "JAGS") c(params, "inv_D", "b") else c(params, "D", "b")

inits <- function () {
  ints <- lapply(components[grep("ncx", names(components), fixed = TRUE)],
                 rnorm, sd = 0.1)
  names(ints) <- paste0('betas', seq_len(n_outcomes))
  ints$u <- drop(matrix(rnorm(Data$n * Data$n_RE), Data$n,
                        Data$n_RE))
  ints$inv_D <- ints$D <- if (Data$n_RE > 1) diag(Data$n_RE) else 1
  if (any(ind_gs)) {
    nms <- which(ind_gs)
    taus <- rep(list(1), length(nms))
    names(taus) <- paste0("tau", nms)
    ints <- c(ints, taus)
  }
  ints
}
fit <- if (engine == "JAGS") {
  jags(data = Data, inits = inits, parameters.to.save = params,
       model.file = file.path(con$working.directory, model_name),
       parallel = con$n.processors > 1, n.chains = con$n.chains,
       n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
       n.thin = con$n.thin, seed = con$seed, verbose = con$verbose)
} else {
  options(mc.cores = con$n.chains)
  if (con$optimize_only) {
    model <- rstan::stan_model(file = file.path(con$working.directory, model_name))
    if (is.null(init)) {
      init <- "random"
    } else {
      if (!is.list(init)) {
        stop("'init' should be a list with appropriate names; run once ",
             "mvglmer() and see the components of object$par.")
      }
    }
    out <- rstan::optimizing(model, data = Data, hessian = TRUE, as_vector = FALSE,
                             seed = con$n.iter, init = init)
    if (con$clear.model) {
      file.remove(file.path(con$working.directory, model_name))
    }
    return(out)
  }
  out <- rstan::stan(file = file.path(con$working.directory, model_name), data = Data, 
                     pars = params, iter = con$n.iter, chains = con$n.chains, 
                     thin = con$n.thin, seed = con$seed, 
                     control = list('adapt_delta' = con$adapt_delta))
  sims.list <- lapply(lapply(params, extract, object = out, permuted = FALSE), bind_chains)
  sims.list[] <- lapply(sims.list, function (x) 
    if (length(dim(x)) == 1) as.matrix(x) else x)
  names(sims.list) <- params
  sims.list[['D']] <- fix_D(sims.list[['D']])
  sims.list[['b']] <- fix_b(sims.list[['b']], Data$n_RE)
  splts <- rep(seq_along(sims.list), sapply(sims.list, function (x) prod(dim(x)[-1])))
  rhats <- head(rstan::summary(out)$summary[, "Rhat"], -1)
  Rhat <- split(rhats, splts)
  names(Rhat) <- params
  list(sims.list = sims.list, Rhat = Rhat,
       mcmc.info = list(n.chains = con$n.chains, n.thin = con$n.thin,
                        n.warmup = con$n.warmup,
                        n.samples = nrow(sims.list[["betas1"]]),
                        elapsed.mins = sum(get_elapsed_time(out)) / 60), 
       DIC = NULL, pD = NULL)
}
if (con$clear.model) {
  file.remove(file.path(con$working.directory, model_name))
}
out <- list(mcmc = fit$sims.list, components = components, data = data,
            families = families, control = con, mcmc.info = fit$mcmc.info,
            DIC = fit$DIC, pD = fit$pD, Rhat = fit$Rhat,
            priors = prs, engine = engine)
if (Data$n_RE == 1) {
  if (engine == "JAGS") 
    out$mcmc$inv_D <- array(out$mcmc$inv_D, c(length(out$mcmc$inv_D), 1, 1))
  else 
    out$mcmc$D <- array(out$mcmc$D, c(length(out$mcmc$D), 1, 1))
  out$mcmc$b <- array(out$mcmc$b, c(nrow(out$mcmc$b), ncol(out$mcmc$b), 1))
}
if (engine == "JAGS") {
  out$mcmc$D <- out$mcmc$inv_D
  for (i in seq_len(nrow(out$mcmc$betas1))) {
    out$mcmc$D[i, , ] <- solve(out$mcmc$D[i, , ])
  }
} else {
  out$mcmc$inv_D <- out$mcmc$D
  for (i in seq_len(nrow(out$mcmc$betas1))) {
    out$mcmc$inv_D[i, , ] <- solve(out$mcmc$D[i, , ])
  }
}
# fix names
Xnams <- lapply(components[grep("^X[0-9]", names(components))], colnames)
for (i in seq_along(Xnams)) {
  colnames(out$mcmc[[paste0("betas", i)]]) <- Xnams[[i]]
}
Znams <- lapply(components[grep("^Z_[0-9]", names(components))], colnames)
Znams <- unlist(mapply(paste0, Znams, seq_len(n_outcomes), SIMPLIFY = FALSE))
dimnames(out$mcmc$inv_D) <- list(NULL, Znams, Znams)
dimnames(out$mcmc$D) <- list(NULL, Znams, Znams) 
dimnames(out$mcmc$b) <- list(NULL, NULL, Znams)
# calculate statistics
summary_fun <- function (FUN, ...) {
  out <- lapply(out$mcmc, function (x) {
    if (!is.null(dim(x)) && length(dim(x)) > 1) {
      d <- if (is.matrix(x)) 2L else c(2L, 3L)
      apply(x, d, FUN, ...)
    } else {
      FUN(x, ...)
    }
  })
  out[!sapply(out, is.null)]
}
out$postMeans <- summary_fun(mean, na.rm = TRUE)
out$postModes <- summary_fun(modes)
out$EffectiveSize <- summary_fun(effectiveSize)
out$StErr <- summary_fun(stdErr)
out$StDev <- summary_fun(sd, na.rm = TRUE)
out$CIs <- summary_fun(quantile, probs = c(0.025, 0.975))
out$Pvalues <- summary_fun(computeP)
out$call <- cl
class(out) <- "mvglmer"
out
}



