#bootRegress: a function with input arguments cleanDat.unreg, traits=numericMeta, regvars (a vector of string values that specify colnames in traits), protectVars (an optional argument that specifies regvars to be included in the regression model but not subtracted in the normExpr.reg output matrix, and therefore protected), and varType (a required input vector specifying either 'factor' or 'numeric', for each regvars value).
#             The function handles regvars of varType 'factor' (e.g., condition, Batch, and Sex) as hard-coded with a reference level first alphabetical (or 'Control' or 'CT' is set to reference if found), whereas 'numeric' varType regvars are handled as, e.g., Age.
#             If a variable  in regvars is also in protectVars, the coefficient times variable is not subtracted in the foreach call that generates normExpr.reg, whereas if it is an unprotected regression variable, it is subtracted, either as a single term for regvars of varType 'numeric', or as each of the non-reference level dummy (binary) variables for each factor level times the coefficient for that dummy variable.
#
#              The function returns regressed expression as a matrix, but first outputs to the console the quantiles check for before and after regression on a random, single sample (column) of the data in cleanDat.unreg and normExpr.reg, respectively.
#
#              A progress bar for coefficient matrix foreach parallel calculation is displayed (optional, default). Each parallel worker of this foreach loop, sets the RNG seed to 8675309 for reproducibility.
#              The function checks for the existence of a doFuture or doParallel backend, and:
#                (1) if doParallel, resets the progress option to FALSE and prints a message about this to the console;
#                (2) if doFuture backend is active, sets the future.globals.maxSize option to 1.5*1e9 and print a message about this to the console;
#                (3) if no parallel backend for multithreading is enabled, sets up a doFuture backend using parallel::detectCores()-1 threads, and prints a message about this to the console.
#
# - by Eric Dammer and Vivek Swarup for use in the Seyfried Lab Multiomics Pipeline

bootRegress <- function(
  cleanDat.unreg,
  traits,  # numericMeta
  regvars, # colnames of numericMeta to include in regression model
  varType, # "numeric" or "factor" for each regvar above, in same order
  protectVars = character(0), # if any regvars' covariance should not be subtracted (protect them), include them here as a vector.
  numboot = 1000,
  progress = TRUE,            # show progress bar? requires doFuture package and backend setup outside or within function.
  progressIncrement = 30      # how often to update progress bar, 1 is after every worker completion.
) {
  # ---- Validation ----
  stopifnot(is.matrix(cleanDat.unreg) || is.data.frame(cleanDat.unreg))
  cleanDat.unreg <- as.matrix(cleanDat.unreg)

  if (missing(traits)) stop("Argument 'traits' is required (e.g., numericMeta).")
  if (length(regvars) == 0) stop("'regvars' must be a non-empty character vector calling out colname(s) of traits.")
  if (missing(varType)) stop("'varType' must be provided and match 'regvars' in length & order. Valid values are 'numeric' or 'factor'.")
  if (length(varType) != length(regvars)) stop("'varType' length must equal 'regvars' length.")
  if (!all(varType %in% c("factor", "numeric")))
    stop("Each 'varType' entry must be either 'factor' or 'numeric'.")

  # ---- Packages (lazy-load) ----
  req <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is required.", p))
  req("foreach"); req("boot")
  library("foreach")
  has_doParallel <- requireNamespace("doParallel", quietly = TRUE)
  has_doFuture   <- requireNamespace("doFuture", quietly = TRUE)
  has_future     <- requireNamespace("future", quietly = TRUE)
  has_progressr  <- requireNamespace("progressr", quietly = TRUE)

  # ---- Backend detection / setup (SAFE) ----
  req <- function(p) if (!requireNamespace(p, quietly = TRUE)) FALSE else TRUE
  has_doParallel <- req("doParallel")
  has_doFuture   <- req("doFuture")
  has_future     <- req("future")
  has_progressr  <- req("progressr")
  
  safe_get_name <- function() {
    nm <- tryCatch(foreach::getDoParName(), error = function(e) NA_character_)
    if (is.null(nm)) NA_character_ else nm
  }
  safe_get_workers <- function() {
    w <- tryCatch(foreach::getDoParWorkers(), error = function(e) NA_integer_)
    if (is.null(w) || is.na(w)) 0L else as.integer(w)
  }
  
  backend_name <- safe_get_name()
  workers      <- safe_get_workers()
  
  # No backend (or sequential) -> set up doFuture multisession
  if (is.na(backend_name) || backend_name %in% c("", "doSEQ")) {
    if (!has_doFuture || !has_future) {
      stop(
        "No parallel backend configured. ",
        "Please install 'doFuture' and 'future' (e.g., install.packages(c('doFuture','future'))), ",
        "or register a backend manually (doParallel/doFuture)."
      )
    }
    doFuture::registerDoFuture()
    n_cores <- max(1L, parallel::detectCores() - 1L)
    future::plan(future::multisession, workers = n_cores)
    backend_name <- safe_get_name()
    workers      <- safe_get_workers()
    message(sprintf("[backend] No backend detected; configured doFuture multisession with %d workers.", workers))
    options(future.globals.maxSize = 1.5 * 1e9)
    message("[backend] Set options(future.globals.maxSize = 1.5e9).")
  
  # doFuture already active
  } else if (identical(backend_name, "doFuture")) {
    if (workers <= 0L && has_future) {
      n_cores <- max(1L, parallel::detectCores() - 1L)
      future::plan(future::multisession, workers = n_cores)
      workers <- safe_get_workers()
    }
    message(sprintf("[backend] doFuture backend detected with %d workers.", workers))
    options(future.globals.maxSize = 1.5 * 1e9)
    message("[backend] Set options(future.globals.maxSize = 1.5e9).")
  
  # doParallel active -> disable progressr
  } else if (identical(backend_name, "doParallel")) {
    message(sprintf("[backend] doParallel backend detected with %d workers.", workers))
    if (isTRUE(progress)) {
      progress <<- FALSE
      message("[progress] Disabled progress bar because 'progressr' does not receive updates from doParallel workers.")
    }
  
  # Some other backend
  } else {
    message(sprintf("[backend] '%s' backend detected with %d workers.", backend_name, workers))
  }

  # ---- Build covariate data with requested types ----
  if (!all(regvars %in% colnames(traits))) {
    miss <- setdiff(regvars, colnames(traits))
    stop("The following 'regvars' are missing in 'traits': ", paste(miss, collapse = ", "))
  }

  covar_df <- setNames(
    lapply(seq_along(regvars), function(j) {
      v <- regvars[j]
      if (varType[j] == "factor") {
        factor(traits[[v]])
      } else {
        as.numeric(traits[[v]])
      }
    }),
    regvars
  )
  covar_df <- as.data.frame(covar_df, check.names = FALSE)

  # ---- Relevel protected factor variables to Control/CT (if present) ----
  protectVars <- intersect(protectVars, regvars)
  for (v in intersect(protectVars, regvars[varType == "factor"])) {
    f <- factor(covar_df[[v]])
    lev <- levels(f)
    target <- if ("Control" %in% lev) "Control" else if ("CT" %in% lev) "CT" else NA_character_
    if (!is.na(target)) {
      covar_df[[v]] <- stats::relevel(f, ref = target)
      message(sprintf("[relevel] Protected factor '%s' re-leveled to reference '%s' automatically.", v, target))
    }
  }

  # ---- Print reference levels for all factor variables ----
  for (v in regvars[varType == "factor"]) {
    ref_level <- levels(covar_df[[v]])[1]
    message(sprintf("[factor] Variable '%s' uses reference level '%s'.", v, ref_level))
  }

  # ---- Prototype design to lock coefficient layout ----
  proto_thisexp <- as.numeric(cleanDat.unreg[1, ])
  # formula: ~ thisexp + var1 + var2 + ...
  rhs <- paste(regvars, collapse = " + ")
  fml <- stats::as.formula(paste0("~ thisexp + ", rhs))
  proto_df <- stats::model.matrix(fml, data = cbind(data.frame(thisexp = proto_thisexp, check.names = FALSE), covar_df))
  design_colnames <- colnames(proto_df)

  # Map numeric and factor columns for later subtraction
  numeric_vars <- regvars[varType == "numeric"]
  factor_vars  <- regvars[varType == "factor"]

  idx_numeric <- setNames(
    vapply(numeric_vars, function(v) match(v, design_colnames, nomatch = 0), integer(1L)),
    numeric_vars
  )

  factor_dummy_cols <- lapply(factor_vars, function(v) {
    # columns begin with 'VarLevel' pattern from model.matrix
    grep(paste0("^", gsub("([\\W])", "\\\\\\1", v)), design_colnames, value = TRUE)
  })
  names(factor_dummy_cols) <- factor_vars

  # Prepare full dummy matrices for factors (for subtraction across all samples)
  factor_dummy_mats <- lapply(factor_vars, function(v) {
    fv <- factor(covar_df[[v]])
    mm <- stats::model.matrix(~ 0 + f, data = data.frame(f = fv))
    colnames(mm) <- paste0(v, levels(fv))
    mm
  })
  names(factor_dummy_mats) <- factor_vars

  # Numeric covariate vectors
  numeric_vectors <- lapply(numeric_vars, function(v) as.numeric(covar_df[[v]]))
  names(numeric_vectors) <- numeric_vars

  # ---- Bootstrap helper ----
  bs_fun <- function(formula, data, indices) {
    d <- data[indices, ]
    fit <- stats::lm(formula, data = d)
    stats::coef(fit)
  }

  # ---- ETA banner based on your benchmark (2736x89 -> 120 min) ----
  n_rows <- nrow(cleanDat.unreg); n_cols <- ncol(cleanDat.unreg)
  if (workers <= 0) workers <- 1L
  baseline_minutes <- 120
  est_minutes <- baseline_minutes * (n_rows / 2736) * (n_cols / 89) * (1 / workers)
  message(sprintf("[estimate] Working on %d rows x %d columns using %d workers.", n_rows, n_cols, workers))
  message(sprintf("[estimate] Estimated time to complete: ~%.1f minutes (scaled from 2736x89 / 120 min).", est_minutes))

  # ---- Progress setup ----
  use_progress <- isTRUE(progress) && has_progressr && (backend_name == "doFuture")
  if (use_progress) {
    # Choose a good default handler (shows ETA in RStudio/console)
    progressr::handlers(progressr::handler_progress(
      format = "[:bar] :percent | elapsed: :elapsed | ETA: :eta"
    ))
    options(progressr.enable = TRUE, progressr.force_handlers = TRUE)
  }

  # ---- Run row-wise bootstrap to produce coefmat (WITH OPTIONAL PROGRESS) ----
  start_time <- Sys.time()
  message(sprintf("Regression Start: %s", format(start_time, "%H:%M:%S")))

  coefmat <- NULL

  if (use_progress) {
    coefmat <- suppressWarnings(progressr::with_progress({
      p <- progressr::progressor(steps = ceiling(n_rows / max(1L, progressIncrement)))
      foreach::foreach(i = 1:n_rows, .combine = rbind,
                       .options.future = list(seed = TRUE)) %dopar% {
        # Worker RNG as requested
        set.seed(8675309)
        thisexp <- as.numeric(cleanDat.unreg[i, ])
        df_i <- stats::model.matrix(fml, data = cbind(data.frame(thisexp = thisexp, check.names = FALSE), covar_df))

        # Column alignment (if levels drop in this row)
        if (!identical(colnames(df_i), design_colnames)) {
          missing_cols <- setdiff(design_colnames, colnames(df_i))
          if (length(missing_cols)) {
            add <- matrix(0, nrow = nrow(df_i), ncol = length(missing_cols),
                          dimnames = list(NULL, missing_cols))
            df_i <- cbind(df_i, add)
          }
          df_i <- df_i[, design_colnames, drop = FALSE]
        }

        bs.res <- boot::boot(
          data = as.data.frame(df_i),
          statistic = bs_fun,
          R = numboot,
          formula = thisexp ~ .
        )
        # tick progress sparingly
        if (i %% progressIncrement == 0L || i == n_rows || i == 1L) p()

        apply(bs.res$t, 2, function(x) stats::median(x, na.rm = TRUE))
      }
    }) )
  } else {
    coefmat <- suppressWarnings( foreach::foreach(i = 1:n_rows, .combine = rbind,
                                .options.future = list(seed = TRUE)) ) %dopar% {
      set.seed(8675309)
      thisexp <- as.numeric(cleanDat.unreg[i, ])
      df_i <- stats::model.matrix(fml, data = cbind(data.frame(thisexp = thisexp, check.names = FALSE), covar_df))

      if (!identical(colnames(df_i), design_colnames)) {
        missing_cols <- setdiff(design_colnames, colnames(df_i))
        if (length(missing_cols)) {
          add <- matrix(0, nrow = nrow(df_i), ncol = length(missing_cols),
                        dimnames = list(NULL, missing_cols))
          df_i <- cbind(df_i, add)
        }
        df_i <- df_i[, design_colnames, drop = FALSE]
      }

      bs.res <- boot::boot(
        data = as.data.frame(df_i),
        statistic = bs_fun,
        R = numboot,
        formula = thisexp ~ .
      )
      apply(bs.res$t, 2, function(x) stats::median(x, na.rm = TRUE))
    }
  }

  colnames(coefmat) <- design_colnames
  rownames(coefmat) <- rownames(cleanDat.unreg)
  coefmat[is.na(coefmat)] <- 0
  coefmat <<- coefmat  # assign to global environment variable of same name, outside of function

  message(sprintf("[step 1] Finished calculating variable median coefficients (%d boostraps per each of %d input rows)\n         - Coefficient matrix saved to global environment as variable 'coefmat'.",numboot,n_rows))

  # ---- Create normalized matrix by subtracting UNPROTECTED effects ----
  message(sprintf("[step 2] Creating normalized matrix by subtracting UNPROTECTED effects. | ETA: <%.2f min", 2/11.5*round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)))

  normExpr.reg <- foreach::foreach(i = 1:n_rows, .combine = rbind) %dopar% {
    row_vec <- as.numeric(cleanDat.unreg[i, ])
    coefs_i <- coefmat[i, ]

    # numeric vars: single term (skip protected)
    if (length(numeric_vars)) {
      for (v in setdiff(numeric_vars, protectVars)) {
        idx <- idx_numeric[[v]]
        if (!is.na(idx) && idx > 0) {
          row_vec <- row_vec - (coefs_i[idx] * numeric_vectors[[v]])
        }
      }
    }

    # factor vars: sum over non-reference dummy columns (skip protected)
    if (length(factor_vars)) {
      for (v in setdiff(factor_vars, protectVars)) {
        dnames <- factor_dummy_cols[[v]]
        if (length(dnames)) {
          common <- intersect(dnames, colnames(factor_dummy_mats[[v]]))
          if (length(common)) {
            adj <- as.numeric(factor_dummy_mats[[v]][, common, drop = FALSE] %*% coefs_i[common])
            row_vec <- row_vec - adj
          }
        }
      }
    }

    row_vec
  }
  rownames(normExpr.reg) <- rownames(cleanDat.unreg)
  colnames(normExpr.reg) <- colnames(cleanDat.unreg)

  # ---- Sanity check on a random sample (column) ----
  set.seed(8675309)
  randomSample <- sample.int(ncol(cleanDat.unreg), 1)
  cat("# Quantiles BEFORE regression (sample ", colnames(cleanDat.unreg)[randomSample], "):\n", sep = "")
  print(stats::quantile(cleanDat.unreg[, randomSample],
                        c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE))
  cat("# Quantiles AFTER regression  (same sample):\n")
  print(stats::quantile(normExpr.reg[, randomSample],
                        c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE))

  # ---- End banner ----
  end_time <- Sys.time()
  elapsed_m <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  message(sprintf("Regression End: %s (%.2f minutes)", format(end_time, "%H:%M:%S"), elapsed_m))

  return(normExpr.reg)
}
