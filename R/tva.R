#'@importFrom rstan extract stan_model sampling optimizing gqs extract sflist2stanfit rstan_options read_stan_csv
#'@importFrom dplyr summarize mutate group_by %>% across if_else select bind_cols bind_rows rename
#'@importFrom tidyr pivot_longer pivot_wider crossing
#'@importFrom readr read_table write_tsv
#'@importFrom methods formalArgs new as callNextMethod show
#'@importFrom stats na.omit as.formula model.matrix pnorm runif terms
#'@importFrom cli col_cyan col_magenta col_grey col_blue ansi_strwrap style_underline style_bold
#'@importFrom tibble tibble
#'@importFrom utils citation str combn
#'@importFrom lme4 findbars subbars fixef ranef nobars
#'@importFrom gtools inv.logit
#'@importFrom cmdstanr cmdstan_model write_stan_file
#'


nolhs <- function(f) {
  stopifnot(inherits(f, "formula") || (is.call(f) && f[[1]] == "~"))
  if(length(f) == 2L) f else f[-2L]
}

rhs <- function(f) {
  stopifnot(inherits(f, "formula") || (is.call(f) && f[[1]] == "~"))
  if(length(f) == 2L) f[[2L]] else f[[3L]]
}

barnames <- function(exprs) {
  vapply(exprs, function(e) {
    stopifnot(is.call(e) && e[[1L]] == "|")
    as.character(as.expression(e[[3L]]))
  }, character(1L))
}




#'@export
stantva_path <- function() {
  file.path(find.package("RStanTVA"), "StanTVA")
}

#'@export
read_tva_data <- function(file, ...) {
  if(inherits(file, "connection")) f <- file
  else f <- base::file(file, "rb")
  n <- as.integer(readLines(f, n = 1))
  dat <- read_table(f, col_names = c("condition", "exposure", "targets", "distractors", "report"), col_types = c("inccc"), n_max = n, na = character(), ...)
  if(n != nrow(dat)) warning("Expected ",n," trials but only read ",nrow(dat)," lines!")
  close(f)
  dat %>% mutate(
    across(c(targets, distractors), ~do.call(rbind, strsplit(.x, "", TRUE) %>% lapply(function(x) if_else(x=="0",NA_character_,x)))),
    report = strsplit(if_else(report == "-", "", report), "", TRUE),
    S = (!is.na(targets) | !is.na(distractors))+0L,
    D = (!is.na(distractors) )+ 0L,
    R = t(vapply(seq_len(n()), function(i) targets[i,] %in% report[[i]] | distractors[i,] %in% report[[i]], logical(ncol(targets)))) + 0L,
    items = t(vapply(seq_len(n()), function(i) if_else(is.na(targets[i,]), distractors[i,], targets[i,]), character(ncol(targets))))
  ) %>% select(condition, S, D, items, T = exposure, R) %>% as("tvadata")
}


#'@export
write_tva_data <- function(data, file, ...) {
  if(inherits(file, "connection")) f <- file
  else f <- base::file(file, "w")
  stopifnot(!is.null(data$S), !is.null(data$R))
  if(is.null(data$items)) {
    data$items <- t(vapply(seq_len(nrow(data)), function(i) sample(LETTERS, ncol(data$S)), character(ncol(data$S))))
  }
  if(is.null(data$condition)) {
    data$condition <- rep(1L, nrow(data))
  }
  if(is.null(data$D)) {
    data$D <- matrix(0L, nrow = nrow(data), ncol = ncol(data$S))
  }
  writeLines(as.character(nrow(data)), f)
  bind_cols(
    condition = data$condition,
    exposure = data$T,
    targets = vapply(seq_len(nrow(data)), function(i) paste0(if_else(data$S[i,] == 1L & data$D[i,] == 0L, data$items[i,], "0"), collapse = ""), character(1)),
    distractors = vapply(seq_len(nrow(data)), function(i) paste0(if_else(data$S[i,] == 1L & data$D[i,] == 1L, data$items[i,], "0"), collapse = ""), character(1)),
    report = vapply(seq_len(nrow(data)), function(i) if(sum(data$R[i,])>0) paste0(if_else(data$S[i,] == 1L & data$D[i,] == 0L & data$R[i,] == 1L, data$items[i,], ""), collapse = "") else "-", character(1))
  ) %>% write_tsv(f, col_names = FALSE, ...)
}

setGeneric("model_code", function(object, ...) {})

#'@export
setMethod("model_code", "stanmodel", function(object, type = c("stan","stan2","cpp")) {
  type <- match.arg(type)
  if(type == "stan") object@code
  else if(type == "stan2") object@model_code
  else if(type == "cpp") object@model_cpp
  else stop("Unknown type '", type,"'!")
})

#'@importClassesFrom rstan stanfit
#'@export
setMethod("model_code", "stanfit", function(object, type) {
  model_code(object@stanmodel, type)
})


## add hierarchical stuff here

prepare_data <- function(trials, model, require_outcome = TRUE) {
  fs <- if(inherits(model, "stantvamodel")) model@code@config$formula else if(inherits(model, "stantvacode")) model@config$formula else stop("`model` must be a StanTVA model or StanTVA model code object!")

  required_columns <- if(isFALSE(require_outcome)) c() else if(model@code@config$task == "wr") c("R","S","T") else if(model@code@config$task == "pr") c("R","S","D","T") else stop("Unsupported task ",sQuote(model@code@config$task),"!")

  if(any(!required_columns %in% colnames(trials))) {
    stop("The supplied data must contain columns ",paste(sQuote(required_columns), collapse=", ")," but at least one is missing: ",paste(sQuote(setdiff(required_columns, colnames(trials))), collapse=", "))
  }

  ltrials <- as.list(trials[,required_columns,drop=FALSE])

  pf <- bind_rows(lapply(fs, parse_formula))

  ltrials$N <- nrow(trials)


  for(i in seq_len(nrow(pf))) {
    for(j in seq_len(nrow(pf$random[[i]]))) {
      x <- eval(pf$random[[i]]$factor[[j]], trials)
      ltrials[[pf$random[[i]]$factor_txt[[j]]]] <- if(is.character(x)) {
        rfc <- as.factor(x)
        rfcl <- levels(rfc)
        rfc <- as.integer(rfc)
        attr(rfc, "levels") <- rfcl
        rfc
      } else if(is.factor(x)) {
        rfc <- as.integer(x)
        attr(rfc, "levels") <- levels(x)
        rfc
      } else if(is.integer(x)) {
        rfc <- x
        attr(rfc, "levels") <- seq_len(max(x))
        rfc
      } else {
        stop("Random factor ",sQuote(deparse1(pf$random[[i]]$factor[[j]]))," must be integer, factor, or character!")
      }
      ltrials[[paste0("N_",pf$random[[i]]$factor_txt[[j]])]] <- length(attr(rfc, "levels"))
    }

    f_var <- pf$param[i]

    Cmatf <- model.matrix(pf$fixed_formula[[i]], trials)
    colnames(Cmatf) <- clean_name(colnames(Cmatf))

    ltrials[[paste0("M_",f_var)]] <- ncol(Cmatf)
    int_C <- which(attr(Cmatf, "assign") == 0L)
    ltrials[[paste0("int_",f_var)]] <- if(length(int_C) != 1L) 0L else int_C
    for(j in seq_len(model@code@df[f_var])) {
      ltrials$X <- cbind(ltrials$X, Cmatf)
      ltrials[[if(model@code@df[f_var] == 1L) paste0("map_",f_var) else paste0("map_",f_var,"_",j)]] <- array((ncol(ltrials$X)-ncol(Cmatf)+1L):ncol(ltrials$X), dim = ncol(Cmatf))
    }

    for(k in seq_len(nrow(pf$random[[i]]))) {
      rf <- pf$random[[i]]$factor_txt[k]
      Cmatr <- model.matrix(pf$random[[i]]$formula[[k]], trials)
      colnames(Cmatr) <- clean_name(colnames(Cmatr))
      for(j in seq_len(model@code@df[f_var])) {
        g <- if(model@code@df[f_var] == 1L || pf$random[[i]]$custom[k]) pf$random[[i]]$group[k] else paste0(pf$random[[i]]$group[k], "_", j)
        ltrials[[paste0("M_",f_var,"_",g)]] <- ncol(Cmatr)
        ltrials[[paste0("Z_",g)]] <- cbind(ltrials[[paste0("Z_",g)]], Cmatr)
        ltrials[[if(model@code@df[f_var] == 1L) paste0("map_",f_var,"_",g) else paste0("map_",f_var,"_",j,"_",g)]] <- array((ncol(ltrials[[paste0("Z_",g)]])-ncol(Cmatr)+1L):ncol(ltrials[[paste0("Z_",g)]]), dim = ncol(Cmatr))
      }
    }
  }

  ltrials
}

clean_name <- function(str) gsub("(^_+)|(_+$)$","",gsub("[^a-zA-Z0-9]+","_",str))


nested_parameter <- function(param, transform, type = "real", dim = 1L, prior_fixed_intercept = NULL) {
  ret <- list(
    target = param,
    is_simplex = grepl("^simplex\\b", type),
    is_vector = grepl("^(vector|simplex)\\b", type),
    dim = dim,
    type_constraint = gsub("^[^<]*(<[^>]+>)?.*$", "\\1", type)
  )
  ret$fdim <- if(ret$is_simplex) dim - 1L else dim

  if(!is.null(prior_fixed_intercept)) {
    ret$prior <- c(
      sprintf("if(int_%1$s) {", param),
      sprintf("\t%2$s ~ %1$s;", as.character(as.expression(rhs(prior_fixed_intercept))), sprintf(transform, if(dim == 1L) sprintf("b[map_%1$s[int_%1$s]]", param) else sprintf("b[map_%1$s_%2$d[int_%1$s]]", param, seq_len(ret$fdim)))),
      "}"
    )
  }

  ret
}


parse_formula <- function(f) {
  stopifnot(is.call(f) && f[[1L]] == "~" && length(f) == 3L)
  lhs <- f[[2L]]
  if(is.symbol(lhs)) {
    parname <- as.character(lhs)
    link_name <- "identity"
    link <- function(.) .
    stan_link <- "%s"
    inverse_link <- function(.) .
    stan_inverse_link <- "%s"
    args <- list()
  } else if(is.call(lhs) && length(lhs) >= 2L && is.symbol(lhs[[2L]])) {
    parname <- as.character(lhs[[2L]])
    link_name <- as.character(lhs[[1L]])
    args <- as.list(lhs[c(-1,-2)])
    if(link_name == "log") {
      link <- log
      stan_link <- "log(%s)"
      inverse_link <- exp
      stan_inverse_link <- "exp(%s)"
    } else if(link_name == "logit") {
      link <- logit
      stan_link <- "logit(%s)"
      inverse_link <- inv.logit
      stan_inverse_link <- "inv_logit(%s)"
    } else if(link_name == "probit") {
      link <- qnorm
      stan_link <- "inv_Phi(%s)"
      inverse_link <- pnorm
      stan_inverse_link <- "Phi(%s)"
    } else if(link_name == "scaled_logit" && length(args) == 2L) {
      link <- function(x) logit((x-args[[1]])/(args[[2]]-args[[1]]))
      stan_link <- sprintf("logit((%%s-%1$g)/(%2$g-%1$g))", args[[1]], args[[2]])
      inverse_link <- function(x) args[[1]]+(args[[2]]-args[[1]])*inv.logit(x)
      stan_inverse_link <- sprintf("%1$g+(%2$g-%1$g)*inv_logit(%%s)", args[[1]], args[[2]])
    } else if(link_name == "scaled_logit" && length(args) == 1L) {
      link <- function(x) logit(x/args[[2]])
      stan_link <- sprintf("logit(%%s/%g)", args[[1]])
      inverse_link <- function(x) args[[2]]*inv.logit(x)
      stan_inverse_link <- sprintf("%g*inv_logit(%%s)", args[[1]])
    } else if(link_name == "scaled_probit" && length(args) == 2L) {
      link <- function(x) qnorm((x-args[[1]])/(args[[2]]-args[[1]]))
      stan_link <- sprintf("inv_Phi((%%s-%1$g)/(%2$g-%1$g))", args[[1]], args[[2]])
      inverse_link <- function(x) args[[1]]+(args[[2]]-args[[1]])*pnorm(x)
      stan_inverse_link <- sprintf("%1$g+(%2$g-%1$g)*Phi(%%s)", args[[1]], args[[2]])
    } else if(link_name == "scaled_probit" && length(args) == 1L) {
      link <- function(x) qnorm(x/args[[2]])
      stan_link <- sprintf("inv_Phi(%%s/%g)", args[[1]])
      inverse_link <- function(x) args[[2]]*pnorm(x)
      stan_inverse_link <- sprintf("%g*Phi(%%s)", args[[1]])
    } else stop("Unknown link ",sQuote(link_name),"!")
  } else stop("lhs must be a call or symbol!")
  fxs <- findbars(f[[3]])
  tibble(
    param = parname,
    link_name = link_name,
    link = list(link),
    stan_link = stan_link,
    inverse_link = list(inverse_link),
    stan_inverse_link = stan_inverse_link,
    fixed_formula = list(formula(call("~",nobars(f[[3]])))),
    random = lapply(seq_along(fxs), function(i) {
      fx <- fxs[[i]]
      if(is.call(fx[[2L]]) && fx[[2L]][[1L]] == "|") {
        tibble(formula = list(formula(call("~",fx[[2L]][[2L]]))), custom = TRUE, group = clean_name(paste0(deparse1(fx[[2L]][[3L]]),"_",deparse1(fx[[3L]]))), factor = list(fx[[3L]]))
      } else {
        tibble(formula = list(formula(call("~",fx[[2L]]))), custom = FALSE, group = clean_name(paste0(parname,"_",i,"_",deparse1(fx[[3L]]))), factor = list(fx[[3L]]))
      }
    }) %>% bind_rows() %>% {
      if(nrow(.) > 0) {
        .$param <- parname
        .$factor_txt <- vapply(.$factor, deparse1, character(1))
      }
      .
    } %>% list()
  )

}


#'@export
stantva_code <- function(formula = NULL, locations, task = c("wr","pr"), regions = list(), C_mode = c("equal","locations","regions"), w_mode = c("locations","regions","equal"), t0_mode = c("constant", "gaussian", "gamma", "exponential"), K_mode = c("bernoulli", "free", "binomial", "hypergeometric"), parallel = isTRUE(rstan_options("threads_per_chain") > 1L), save_log_lik = FALSE, priors = TRUE, sanity_checks = TRUE) {

  task <- match.arg(task)
  C_mode <- match.arg(C_mode)
  t0_mode <- match.arg(t0_mode)
  K_mode <- match.arg(K_mode)
  w_mode <- match.arg(w_mode)

  hierarchical_config <- bind_rows(lapply(formula, parse_formula))

  prior_random_corr <- ~lkj_corr(1)
  prior_random_sd <- ~std_normal()

  includeFile <- function(f) sprintf("#include %s", f)

  call_args <- character(0)
  call_args_list <- list()

  for(x in formalArgs(sys.function())) {
    if(x %in% c("data","type")) next
    val <- get(x)
    call_args[x] <- if(is.numeric(val) || is.character(val) || is.logical(val)) as.character(val) else deparse1(val)
    call_args_list[[x]] <- val
  }


  code_blocks <- list(
    `functions` = character(0),
    `data` = character(0),
    `transformed data` = character(0),
    `parameters` = character(0),
    `transformed parameters` = character(0),
    `model` = character(0),
    `generated quantities` = character(0)
  )


  parameters <- list()
  mydata <- list()

  if(!is.null(formula)) {
    if(!is.list(formula)) {
      stop("`formula` must be a list of regression-style formulas for model parameters!")
    }
    formula_lhs <- vapply(formula, function(f) {
      if(!inherits(f, "formula")) stop("All elements of `formula` must be formulas!")
      if(length(f) != 3L) stop("All formulas in `formula` must have a left-hand side!")
      lhs <- f[[2]]
      if(is.call(lhs) && length(lhs) >= 2L && is.symbol(lhs[[2]])) {
        as.character(lhs[1:2])
      } else if(is.symbol(lhs)) {
        c("identity", as.character(lhs))
      } else stop("Left-hand sides of formulas in `formula` must be a single parameter name (e.g., `C`) or a transformed parameter (e.g., `log(C)`)!")
    }, character(2L))

    has_formulas <- TRUE
  } else {
    has_formulas <- FALSE
  }

  add_code <- function(name, ..., prepend = FALSE) code_blocks[[name]] <<- if(isTRUE(prepend))  c(..., code_blocks[[name]]) else c(code_blocks[[name]], ...)
  add_param <- function(name, ...) {
    if(has_formulas && name %in% formula_lhs[2,]) {
      parameters[[name]] <<- list(...)
      parameters[[name]]$class <<- c(setdiff(parameters[[name]]$class, "phi"), "theta")
      hc <- nested_parameter(
        name,
        hierarchical_config$stan_inverse_link[match(name, hierarchical_config$param)],
        parameters[[name]]$type,
        if(is.null(parameters[[name]]$dim)) 1L else parameters[[name]]$dim,
        prior_fixed_intercept = parameters[[name]]$prior
      )
      parameters[[name]]$hierarchical <<- hc
    } else {
      parameters[[name]] <<- list(...)
      add_code(if(isTRUE(parameters[[name]]$transformed)) "transformed parameters" else "parameters", sprintf("%s %s;", parameters[[name]]$type, name))
    }
  }
  add_data <- function(name, ...) {
    mydata[[name]] <<- list(...)
    add_code(if(isTRUE(mydata[[name]]$transformed)) "transformed data" else "data", sprintf("%s %s;", mydata[[name]]$type, name))
  }
  add_region <- function(name, locations) regions[[name]] <<- as.integer(locations)



  if(length(regions) == 1) {
    stop("If you define regions, please define at least two!")
  } else if(length(regions) > 1) {
    all_covered_locations <- unlist(regions)
    duplicate_locations <- unique(all_covered_locations[duplicated(all_covered_locations)])
    if(length(duplicate_locations) > 0) stop("You may not assign location(s) ", paste0(duplicate_locations, collapse=", "), " to multiple regions!")
    if(!setequal(all_covered_locations, seq_len(locations))) stop("You did not cover all locations with regions!")
  }

  datmap <- function(names = base::names(mydata)) {
    names_x_i <- Filter(function(name) !is.null(mydata[[name]]$class) && mydata[[name]]$class == "x_i", names)
    names_x_r <- Filter(function(name) !is.null(mydata[[name]]$class) && mydata[[name]]$class == "x_r", names)
    par_size_x_i <- vapply(mydata[names_x_i], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_i <- c(1L,1L+cumsum(par_size_x_i[-length(par_size_x_i)]))
    par_size_x_r <- vapply(mydata[names_x_r], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_r <- c(1L,1L+cumsum(par_size_x_r[-length(par_size_x_r)]))
    c(
      sprintf("array[N,%d] int x_i;", sum(par_size_x_i)),
      vapply(seq_along(names_x_i), function(i) if(par_size_x_i[i] == 1) sprintf("x_i[:,%d] = %s;", par_offset_x_i[i], names_x_i[i]) else sprintf("x_i[:,%d:%d] = %s;", par_offset_x_i[i], par_offset_x_i[i]+par_size_x_i[i]-1L, names_x_i[i]), character(1)),
      sprintf("array[N,%d] real x_r;", sum(par_size_x_r)),
      vapply(seq_along(names_x_r), function(i) if(par_size_x_r[i] == 1) sprintf("x_r[:,%d] = %s;", par_offset_x_r[i], names_x_r[i]) else sprintf("x_r[:,%d:%d] = %s;", par_offset_x_r[i], par_offset_x_r[i]+par_size_x_r[i]-1L, names_x_r[i]), character(1))
    )
  }

  datremap <- function(names = base::names(mydata), names_back = names) {
    names_x_i <- Filter(function(name) !is.null(mydata[[name]]$class) && "x_i" %in% mydata[[name]]$class, names)
    names_x_r <- Filter(function(name) !is.null(mydata[[name]]$class) && "x_r" %in% mydata[[name]]$class, names)
    which_names_back_x_i <- na.omit(match(names_back, names_x_i))
    which_names_back_x_r <- na.omit(match(names_back, names_x_r))
    par_size_x_i <- vapply(mydata[names_x_i], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_i <- c(1L,1L+cumsum(par_size_x_i[-length(par_size_x_i)]))
    par_size_x_r <- vapply(mydata[names_x_r], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_r <- c(1L,1L+cumsum(par_size_x_r[-length(par_size_x_r)]))
    c(
      vapply(which_names_back_x_i, function(i) if(par_size_x_i[i] == 1) sprintf("x_i[%1$d]", par_offset_x_i[i]) else sprintf("x_i[%1$d:%2$d]", par_offset_x_i[i], par_offset_x_i[i]+par_size_x_i[i]-1L), character(1)),
      vapply(which_names_back_x_r, function(i) if(par_size_x_r[i] == 1) sprintf("x_r[%1$d]", par_offset_x_r[i]) else sprintf("to_vector(x_r[%1$d:%2$d])", par_offset_x_r[i], par_offset_x_r[i]+par_size_x_r[i]-1L), character(1))
    )
  }

  parmap <- function(names = base::names(parameters)) {
    names_theta <- Filter(function(name) !is.null(parameters[[name]]$class) && "theta" %in% parameters[[name]]$class, names)
    names_phi <- Filter(function(name) !is.null(parameters[[name]]) && (is.null(parameters[[name]]$class) || "phi" %in% parameters[[name]]$class), names)
    par_size_theta <- vapply(parameters[names_theta], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_theta <- c(1L,1L+cumsum(par_size_theta[-length(par_size_theta)]))
    par_size_phi <- vapply(parameters[names_phi], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_phi <- c(1L,1L+cumsum(par_size_phi[-length(par_size_phi)]))
    c(
      sprintf("vector[%d] phi;", sum(par_size_phi)),
      vapply(seq_along(names_phi), function(i) if(par_size_phi[i] == 1) sprintf("phi[%d] = %s;", par_offset_phi[i], names_phi[i]) else sprintf("phi[%d:%d] = %s;", par_offset_phi[i], par_offset_phi[i]+par_size_phi[i]-1L, names_phi[i]), character(1)),
      sprintf("array[N] vector[%d] theta;", sum(par_size_theta)),
      unlist(lapply(seq_along(names_theta), function(i) if(par_size_theta[i] == 1) sprintf("theta[:,%d] = to_array_1d(%s);", par_offset_theta[i], names_theta[i]) else sprintf("theta[:,%d] = to_array_1d(%s[,%d]);", seq.int(from = par_offset_theta[i], length.out = par_size_theta[i]), names_theta[i], seq_len(par_size_theta[i]))))
      #vapply(seq_along(names_theta), function(i) if(par_size_theta[i] == 1) sprintf("theta[:,%d] = to_array_1d(%s);", par_offset_theta[i], names_theta[i]) else sprintf("theta[:,%d:%d] = to_vector_array(%s);", par_offset_theta[i], par_offset_theta[i]+par_size_theta[i]-1L, names_theta[i]), character(1))
    )
  }

  parremap <- function(names = base::names(parameters), names_back = names) {
    names_theta <- Filter(function(name) !is.null(parameters[[name]]$class) && "theta" %in% parameters[[name]]$class, names)
    names_phi <- Filter(function(name) !is.null(parameters[[name]]) && (is.null(parameters[[name]]$class) || "phi" %in% parameters[[name]]$class), names)
    which_names_back_phi <- na.omit(match(names_back, names_phi))
    which_names_back_theta <- na.omit(match(names_back, names_theta))
    par_size_theta <- vapply(parameters[names_theta], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_theta <- c(1L,1L+cumsum(par_size_theta[-length(par_size_theta)]))
    par_size_phi <- vapply(parameters[names_phi], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_phi <- c(1L,1L+cumsum(par_size_phi[-length(par_size_phi)]))
    c(
      vapply(which_names_back_phi, function(i) if(par_size_phi[i] == 1) sprintf("phi[%1$d]", par_offset_phi[i]) else sprintf("phi[%1$d:%2$d]", par_offset_phi[i], par_offset_phi[i]+par_size_phi[i]-1L), character(1)),
      vapply(which_names_back_theta, function(i) if(par_size_theta[i] == 1) sprintf("theta[%1$d]", par_offset_theta[i]) else sprintf("theta[%1$d:%2$d]", par_offset_theta[i], par_offset_theta[i]+par_size_theta[i]-1L), character(1))
    )
  }



  parsig <- function(names = base::names(parameters), names_back = names, types = TRUE, index = NULL) {
    names_theta <- Filter(function(name) !is.null(parameters[[name]]$class) && "theta" %in% parameters[[name]]$class, names)
    names_phi <- Filter(function(name) !is.null(parameters[[name]]) && (is.null(parameters[[name]]$class) || "phi" %in% parameters[[name]]$class), names)
    which_names_back_phi <- na.omit(match(names_back, names_phi))
    which_names_back_theta <- na.omit(match(names_back, names_theta))
    par_size_theta <- vapply(parameters[names_theta], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_theta <- c(1L,1L+cumsum(par_size_theta[-length(par_size_theta)]))
    par_size_phi <- vapply(parameters[names_phi], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_phi <- c(1L,1L+cumsum(par_size_phi[-length(par_size_phi)]))
    c(
      vapply(which_names_back_phi, function(i) {
        x <- if(isTRUE(types)) sprintf("%s %s", parameters[[names_phi[i]]]$rtype, names_phi[i]) else names_phi[i]
        if(!isTRUE(types) && !is.null(parameters[[names_phi[i]]]$dim) && parameters[[names_phi[i]]]$dim > 1) x <- sprintf("to_vector(%s)", x)
        x
      }, character(1)),
      vapply(which_names_back_theta, function(i) {
        x <- if(isTRUE(types)) sprintf("%s %s", parameters[[names_theta[i]]]$rtype, names_theta[i]) else if(!is.null(index)) sprintf("%s[%s]", names_theta[i], index) else names_theta[i]
        if(!isTRUE(types) && !is.null(parameters[[names_theta[i]]]$dim) && parameters[[names_theta[i]]]$dim > 1) x <- sprintf("to_vector(%s)", x)
        x
      }, character(1))
    )
  }

  datsig <- function(names = base::names(mydata), names_back = names, types = TRUE, index = NULL) {
    names_x_i <- Filter(function(name) !is.null(mydata[[name]]$class) && "x_i" %in% mydata[[name]]$class, names)
    names_x_r <- Filter(function(name) !is.null(mydata[[name]]$class) && "x_r" %in% mydata[[name]]$class, names)
    which_names_back_x_i <- na.omit(match(names_back, names_x_i))
    which_names_back_x_r <- na.omit(match(names_back, names_x_r))
    par_size_x_i <- vapply(mydata[names_x_i], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_i <- c(1L,1L+cumsum(par_size_x_i[-length(par_size_x_i)]))
    par_size_x_r <- vapply(mydata[names_x_r], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_r <- c(1L,1L+cumsum(par_size_x_r[-length(par_size_x_r)]))
    c(
      vapply(which_names_back_x_i, function(i) if(isTRUE(types)) sprintf("data %s %s", mydata[[names_x_i[i]]]$rtype, names_x_i[i]) else if(!is.null(index)) sprintf("%s[%s]", names_x_i[i], index) else names_x_i[i], character(1)),
      vapply(which_names_back_x_r, function(i) if(isTRUE(types)) sprintf("data %s %s", mydata[[names_x_r[i]]]$rtype, names_x_r[i]) else if(!is.null(index)) sprintf("%s[%s]", names_x_r[i], index) else names_x_r[i], character(1))
    )
  }


  add_code("functions", includeFile("tva.stan"))

  add_data(name = "N", type = "int<lower=1>")

  if(C_mode == "equal") {
    add_param(name = "C", type = "real<lower=machine_precision()>", ctype = "real", rtype = "real", prior = ~gamma(3.5,0.035))
    s_pars <- "C"
    s_body <- c(
      sprintf("vector[%1$d] s = rep_vector(C, %1$d);", locations)
    )
  } else if(C_mode == "regions") {
    if(length(regions) == 0) stop("You must define regions if C_mode == `regions`!")
    s_pars <- "C"
    s_body <- c(
      sprintf("vector[%d] s;", locations),
      unlist(
        lapply(seq_along(regions), function(i) {
          sprintf("s[%d] = C[%d];", regions[[i]], i)
        })
      )
    )
    add_param(name = "C", type = sprintf("vector<lower=machine_precision()>[%d]", length(regions)), ctype=sprintf("vector[%d]", length(regions)), rtype = "vector", dim = length(regions))
    for(i in seq_along(regions)) {
      #add_code(
      #  "generated quantities",
      #  paste0("real C_", names(regions)[i]," = C[",i,"];"),
      #  paste0("real A_", names(regions)[i]," = C_", names(regions)[i],"/sum(C);")
      #)
    }
  } else if(C_mode == "locations") {
    s_body <- NULL
    s_pars <- "s"
    add_param(name = "s", type = sprintf("vector<lower=machine_precision()>[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype="vector", dim = locations, prior = ~gamma(3.5,0.035))
    for(i in seq_along(regions)) {
      #add_code(
      #  "generated quantities",
      #  paste0("real A_", names(regions)[i]," = (", paste0("s[", regions[[i]], "]", collapse=" + "),")/sum(s);")
      #)
    }
  }

  if(w_mode == "equal") {
    w_pars <- c()
    w_body <- sprintf("vector[%1$d] w = rep_vector(1.0/%1$d.0, %1$d);", locations)
  } else if(w_mode == "regions") {
    if(length(regions) == 0) stop("You must define regions if w_mode = 'regions'!")
    add_param(name = "b", type = sprintf("simplex[%d]", length(regions)), ctype=sprintf("vector[%d]", length(regions)), rtype = "vector", dim = length(regions))
    w_pars <- "b"
    w_body <- c(
      sprintf("vector[%1$d] w;", locations),
      vapply(seq_along(regions), function(i) {
        sprintf("w[%d] = b[%d]/%d.0;", regions[[i]], i, length(regions[[i]]))
      }, character(1))
    )
    #add_code(
    #  "generated quantities",
    #  paste0("real w_", names(regions)," = b[", seq_along(regions),"];")
    #)
  } else if(w_mode == "locations") {
    w_pars <- "w"
    w_body <- NULL
    add_param(name = "w", type = sprintf("simplex[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype ="vector", dim = locations)
    for(i in seq_along(regions)) {
      #add_code(
      #  "generated quantities",
      #  paste0("real w_", names(regions)[i]," = ", paste0("w[", regions[[i]], "]", collapse=" + "),";")
      #)
    }
  }

  K_args <- "[]'"
  if(K_mode == "bernoulli") {
    add_code("functions", includeFile("bernoulliK.stan"))
    add_param(name = "K", class = c("phi", "K"), type = sprintf("real<lower=0,upper=%d>", locations), ctype="real", rtype="real", prior = substitute(~uniform(0,nS), list(nS = locations)))
    K_args <- "[K]'"
  } else if(K_mode == "free") {
    add_code("functions", includeFile("freeK.stan"))
    add_param(name = "pK", class = c("phi","K"), type = sprintf("simplex[%d]", locations+1L), ctype=sprintf("vector[%d]", locations+1L), rtype="vector", dim = locations+1L)
    #add_code("generated quantities", paste0("real mK = ",paste(sprintf("%d * pK[%d]", seq_len(locations), seq_len(locations)+1L), collapse=" + "),";"));
    K_args <- "pK"
  } else if(K_mode == "binomial") {
    add_code("functions", includeFile("binomialK.stan"))
    # TODO add prior!
    add_param(name = "nK", class = c("phi", "K"), type = "real<lower=machine_precision()>", ctype="real", rtype="real")
    add_param(name = "pK", class = c("phi", "K"), type = "real<lower=machine_precision(),upper=1.0-machine_precision()>", ctype="real", rtype="real", prior = ~uniform(0,1))
    #add_code("generated quantities", "real mK = nK * pK;")
    K_args <- "[nK, pK]'"
  } else if(K_mode == "hypergeometric") {
    add_code("functions", includeFile("hypergeometricK.stan"))
    add_code(
      "transformed data",
      "int min_gK = 0;",
      "for(i in 1:N) {",
      "\tint nR = sum(R[i,]);",
      "\tif(nR > min_gK) min_gK = nR;",
      "}",
      sprintf("int min_bK = %d - min_gK;", locations)
    )
    # TODO add priors!
    add_param(name = "gK", class = c("phi", "K"), type = "real<lower=min_gK>", ctype="real", rtype="real")
    add_param(name = "bK", class = c("phi", "K"), type = "real<lower=min_bK>", ctype="real", rtype="real")
    #add_code("generated quantities", sprintf("real mK = %d * gK / (gK + bK);", locations))
    K_args <- "[gK, bK]'"
  }

  t0_args <- "[]'"
  if(t0_mode == "constant") {
    add_code("functions", includeFile("constantt0.stan"))
    add_data(name = "max_mu", type = "real", ctype = "real", rtype="real", transformed = TRUE)
    add_code("transformed data", "max_mu = max(T);", "for(i in 1:N) if(sum(R[i,]) && T[i] < max_mu) max_mu = T[i];")
    add_param(name = "t0", class = c("phi"), type = "real<upper=max_mu>", ctype="real", rtype="real", prior = ~normal(20, 30))
  } else if(t0_mode == "gaussian") {
    add_code("functions", includeFile("gaussiant0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real", ctype="real", rtype="real", prior = ~normal(20, 30))
    add_param(name = "sigma0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real", prior = ~gamma(2,0.08))
    t0_args <- "[mu0, sigma0]'"
  } else if(t0_mode == "exponential") {
    # TODO implement default priors!
    add_code("functions", includeFile("exponentialt0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real")
    #add_param(name = "t0", class = c("phi"), type = "real", ctype="real", rtype="real")
    t0_args <- "[1/mu0]'"
  } else if(t0_mode == "gamma") {
    # TODO implement default priors!
    add_code("functions", includeFile("gammat0.stan"))
    add_param(name = "a0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real")
    add_param(name = "b0", class = c("phi","t0"), type = "real<lower=C/1000.0>", ctype="real", rtype="real")
    #add_param(name = "t0", class = c("phi"), type = "real", ctype="real", rtype="real")
    t0_args <- "[a0,b0]'"
  }

  add_data(name = "nS", type = "array[N] int", ctype = "int", rtype="int", class="x_i", transformed = TRUE)
  add_code(
    "transformed data",
    "for(i in 1:N) nS[i] = sum(S[i,]);",
    "int total_nS = sum(nS);"
  )



  add_data(name = "T", class="x_r", type = "array[N] real<lower=0>", ctype = "real", rtype="real")
  add_data(name = "S", class="x_i", type = sprintf("array[N,%d] int<lower=0,upper=1>", locations), ctype = sprintf("array[%d] int", locations), rtype="array[] int", dim = locations)
  add_data(name = "R", class="x_i", type = sprintf("array[N,%d] int<lower=0,upper=1>", locations), ctype = sprintf("array[%d] int", locations), rtype="array[] int", dim = locations)


  if(task == "wr") {
    v_data <- c("nS","S")
    v_pars <- c(w_pars,s_pars)
    v_body <- c(
      w_body,
      s_body,
      "array[nS] int Ss = get_matches(S);",
      "vector[nS] v = s[Ss] .* w[Ss] / sum(w[Ss]);",
      "for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();",
      "return v/1000.0;"
    )
    l_data <- union(c("S","R","T","nS"), v_data)
    l_pars <- c(v_pars, if(!is.null(parameters$t0))"t0",Filter(function(p) any(c("t0","K") %in% parameters[[p]]$class), names(parameters)))
    l_body <- c(
      sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("log_lik = tva_wr_log(R, S, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
    )
    p_data <- setdiff(l_data, "R")
    p_pars <- l_pars
    p_body <- c(
      sprintf("vector[%d] p;", locations+1),
      sprintf("for(i in 0:%d) {", locations-1L),
      sprintf("\tvector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("\tp[i+1] = exp(tva_wr_score_log(i, S, %s, %s, %s, v));", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args),
      "}",
      sprintf("p[%d] = 1.0 - sum(p[:%d]);", locations+1, locations),
      "return p;"
    )
    s_data <- p_data
    s_pars <- p_pars
    s_body <- c(
      sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("return tva_wr_rng(S, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
    )
  } else if(task == "pr") {
    v_data <- c("nS","S","D")
    v_pars <- c(w_pars,s_pars,"alpha")
    v_body <- c(
      w_body,
      s_body,
      "array[nS] int Ss = get_matches(S);",
      sprintf("vector[%d] w_alpha = w;", locations),
      sprintf("for(i in 1:%d) if(D[i]) w_alpha[i] *= alpha;", locations),
      "vector[nS] v = s[Ss] .* w_alpha[Ss] / sum(w_alpha[Ss]);",
      "for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();",
      "return v/1000.0;"
    )
    add_data(name = "D", class="x_i", type = sprintf("array[N,%d] int<lower=0,upper=1>", locations), ctype=sprintf("array[%d] int", locations), rtype="array[] int", dim = locations)
    add_param(name = "alpha", type = "real<lower=machine_precision()>", ctype = "real", rtype="real", prior = ~lognormal(-0.4,0.6))
    l_data <- union(c("S","D","R","T"), v_data)
    l_pars <- c(v_pars,if(!is.null(parameters$t0))"t0",Filter(function(p) any(c("t0","K") %in% parameters[[p]]$class), names(parameters)))
    l_body <- c(
      sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("log_lik = tva_pr_log(R, S, D, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
    )
    p_data <- setdiff(l_data, "R")
    p_pars <- l_pars
    p_body <- c(
      sprintf("vector[%d] p;", locations+1),
      sprintf("for(i in 0:%d) {", locations-1),
      sprintf("\tvector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("\tp[i+1] = exp(tva_pr_score_log(i, S, D, %s, %s, %s, v));", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args),
      "}",
      sprintf("p[%d] = 1.0 - sum(p[:%d]);", locations+1, locations),
      "return p;"
    )
    s_data <- p_data
    s_pars <- p_pars
    s_body <- c(
      sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("return tva_pr_rng(S, D, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
    )
  }

  add_code(
    "functions",
    paste0("vector calculate_v(",paste(c(datsig(names_back = v_data), parsig(v_pars)), collapse = ", "),") {"),
    paste0("\t", v_body),
    "}"
  )

  # hierarchical stuff




  if(nrow(hierarchical_config) > 0L) {
    all_params <- hierarchical_config %>% rename(name = param) %>% mutate(dim = vapply(name, function(x) if(!is.null(parameters[[x]]$dim)) as.integer(parameters[[x]]$dim) else 1L, integer(1)), fdim = vapply(name, function(x) if(grepl("^simplex", parameters[[x]]$type)) as.integer(parameters[[x]]$dim)-1L else if(!is.null(parameters[[x]]$dim)) as.integer(parameters[[x]]$dim) else 1L, integer(1)))
    for(i in seq_len(nrow(all_params))) {
      if(all_params$dim[i] > 1L) {
        is_custom <- all_params$random[[i]]$custom
        all_params$random[[i]] <- crossing(all_params$random[[i]], index = seq_len(all_params$fdim[i]))
        if(!is_custom) {
          all_params$random[[i]]$group <- sprintf("%s_%d", all_params$random[[i]]$group, all_params$random[[i]]$index)
        }
      } else {
        all_params$random[[i]]$index <- NA_integer_
      }
    }
    all_random_effects <- bind_rows(tibble(param = character(), group=character(),factor_txt=character()), all_params$random)
    all_random_factors <- unique(all_random_effects$factor_txt)
    all_random_params <- all_random_effects %>% group_by(group, factor_txt) %>% summarize(dim = sum(all_params$dim[match(param, all_params$name)]), fdim = sum(all_params$fdim[match(param, all_params$name)]), M_var = paste(sprintf("M_%s_%s", param, group), collapse = "+"))
    M_var <- paste(if_else(all_params$fdim!=1L,sprintf("%d*M_%s", all_params$fdim, all_params$name),sprintf("M_%s", all_params$name)), collapse = "+")
    add_code(
      "data",
      sprintf("int<lower=1,upper=N> N_%s;", clean_name(all_random_factors)),
      sprintf("array[N] int<lower=1,upper=N_%1$s> %1$s;", clean_name(all_random_factors)),
      sprintf("int<lower=0> M_%s;", all_params$name),
      unique(sprintf("int<lower=0> M_%s_%s;", all_random_effects$param, all_random_effects$group)),
      sprintf("int<lower=0,upper=M_%1$s> int_%1$s;", all_params$name),
      sprintf("matrix[N,%s] X;", M_var),
      unlist(lapply(seq_len(nrow(all_params)), function(i) if(all_params$dim[i] == 1L) sprintf("array[M_%1$s] int map_%1$s;", all_params$name[i]) else sprintf("array[M_%1$s] int map_%1$s_%2$d;", all_params$name[i], seq_len(all_params$fdim[i])))),
      sprintf("matrix[N,%2$s] Z_%1$s;", all_random_params$group, all_random_params$M_var),
      if_else(is.na(all_random_effects$index), sprintf("array[M_%1$s_%2$s] int map_%1$s_%2$s;", all_random_effects$param, all_random_effects$group), sprintf("array[M_%1$s_%2$s] int map_%1$s_%3$d_%2$s;", all_random_effects$param, all_random_effects$group, all_random_effects$index))
    )
    add_code(
      "transformed data",
      sprintf("vector[%1$d] ones_vector_%1$d = rep_vector(1.0, %1$d);", unique(all_params$fdim[all_params$fdim < all_params$dim])),
      sprintf("int M = %s;", M_var),
      sprintf("int M_%s = %s;", all_random_params$group, all_random_params$M_var),
      sprintf("array[N_%1$s] int<lower=0> Ntrials_by_%1$s = rep_array(0, N_%1$s);", clean_name(all_random_factors)),
      sprintf("array[N_%1$s,N] int<lower=0> trials_by_%1$s = rep_array(0, N_%1$s, N);", clean_name(all_random_factors)),
      if(length(all_random_factors) > 0) "for(i in 1:N) {",
      sprintf("\tNtrials_by_%1$s[%1$s[i]] += 1;", clean_name(all_random_factors)),
      sprintf("\ttrials_by_%1$s[%1$s[i],Ntrials_by_%1$s[%1$s[i]]] = i;", clean_name(all_random_factors)),
      if(length(all_random_factors) > 0) "}"
    )
    add_code(
      "parameters",
      "vector[M] b;",
      sprintf("corr_matrix[M_%1$s] r_%1$s;", all_random_params$group),
      sprintf("vector<lower=machine_precision()>[M_%1$s] s_%1$s;", all_random_params$group),
      sprintf("array[N_%2$s] vector[M_%1$s] w_%1$s;", all_random_params$group, all_random_params$factor_txt)
    )
    add_code(
      "transformed parameters",
      unlist(lapply(seq_len(nrow(all_params)), function(i) {
        if(all_params$dim[i] == 1) {
          sprintf("vector[N] %s;", all_params$name[i])
        } else {
          sprintf("matrix[N,%2$d] %1$s;", all_params$name[i], all_params$dim[i])
        }
      })),
      "{",
      paste0("\t",
        c(
          unlist(lapply(seq_len(nrow(all_params)), function(i) {
            if(all_params$dim[i] == 1) {
              sprintf("%1$s = X[,map_%1$s] * b[map_%1$s];", all_params$name[i])
            } else {
              sprintf("%2$s[,%1$d] = X[,map_%2$s_%1$d] * b[map_%2$s_%1$d];", seq_len(all_params$fdim[i]), all_params$name[i])
            }
          })),
          unlist(lapply(all_random_factors, function(rf) {
            j <- which(all_random_effects$factor_txt == rf)
            c(
              sprintf("for(i in 1:N_%s) {", clean_name(rf)),
              sprintf("\tarray[Ntrials_by_%1$s[i]] int j = trials_by_%1$s[i,:Ntrials_by_%1$s[i]];", clean_name(rf)),
              #sprintf("\tmatrix[Ntrials_by_%1$s[i],M_%2$s] Z_%2$s_i = Z_%2$s[j,:];", clean_name(rf), unique(all_random_effects$group[j])),
              #sprintf("\tvector[M_%1$s] z_%1$s_i = w_%1$s[i,:];", unique(all_random_effects$group[j])),
              unlist(
                lapply(j, function(i) {
                  k <- match(all_random_effects$param[i], all_params$name)
                  if(all_params$dim[k] == 1L) {
                    sprintf("\t%1$s[j] += Z_%2$s[j,map_%1$s_%2$s] * w_%2$s[i,map_%1$s_%2$s];", all_random_effects$param[i], all_random_effects$group[i])
                  } else {
                    sprintf("\t%1$s[j,%3$d] += Z_%2$s[j,map_%1$s_%3$d_%2$s] * w_%2$s[i,map_%1$s_%3$d_%2$s];", all_random_effects$param[i], all_random_effects$group[i], all_random_effects$index[i])
                  }
                })
              ),
              "}"
            )
          })),
          unlist(lapply(which(all_params$link_name != "identity"), function(i) if(all_params$fdim[i] < all_params$dim[i]) sprintf("%1$s[,:%2$d] = %3$s;", all_params$name[i], all_params$fdim[i], sprintf(all_params$stan_inverse_link[i], sprintf("%s[,:%d]", all_params$name[i], all_params$fdim[i]))) else sprintf("%1$s = %2$s;", all_params$name[i], sprintf(all_params$stan_inverse_link[i], all_params$name[i])))),
          unlist(lapply(which(all_params$fdim < all_params$dim), function(i) c(sprintf("%1$s[,%2$d] = 1.0 / (1.0 + %3$s);", all_params$name[i], all_params$dim[i], paste(sprintf("%s[,%d]", all_params$name[i], seq_len(all_params$fdim[i])), collapse = " + ")), sprintf("%1$s[,%3$d] .*= %1$s[,%2$d];", all_params$name[i], all_params$dim[i], seq_len(all_params$fdim[i])))))
        )
      ),
      "}"
    )
    add_code(
      "transformed data",
      sprintf("vector[M_%1$s] mu_w_%1$s = rep_vector(0.0, M_%1$s);", all_random_params$group)
    )
    add_code(
      "model",
      unlist(lapply(seq_len(nrow(all_random_params)), function(i) {
        c(
          if(isTRUE(priors)) sprintf("s_%1$s ~ %2$s;", all_random_params$group[i], as.character(as.expression(rhs(prior_random_sd)))),
          sprintf("if(M_%s > 1) {", all_random_params$group[i]),
          if(isTRUE(priors)) sprintf("\tr_%1$s ~ %2$s;", all_random_params$group[i], as.character(as.expression(rhs(prior_random_corr)))),
          sprintf("\tw_%1$s ~ multi_normal(mu_w_%1$s, quad_form_diag(r_%1$s, s_%1$s));", all_random_params$group[i]),
          "} else {",
          sprintf("\tw_%1$s[,1] ~ normal(0.0, s_%1$s[1]);", all_random_params$group[i]),
          "}"
        )
      }))
    )

  }


  ##



  add_code(
    "functions",
    paste0("real log_lik_single(",paste(c(datsig(names_back = l_data), parsig(l_pars)), collapse = ", "),") {"),
    paste0("\t", c("real log_lik;", l_body, sprintf("if(0 && log_lik == negative_infinity()) print(\"logLik(%s|%s) = -inf !\");", paste0(sprintf("%1$s=\",%1$s,\"",datsig(names_back = l_data, types = FALSE)),collapse=","), paste0(sprintf("%1$s=\",%1$s,\"",parsig(l_pars, types = FALSE)),collapse=",")), "return log_lik;")),
    "}"
  )


  if(isTRUE(parallel)) {
    add_code(
      "functions",
      "vector log_lik_rect(vector phi, vector theta, data array[] real x_r, data array[] int x_i) {",
      paste0("\treturn [log_lik_single(",paste(c(datremap(names_back = l_data), parremap(l_pars)),collapse=", "),")]';"),
      "}"
    )
  }


  add_code(
    "model",
    if(isTRUE(priors)) {
      c(
        "// default priors",
        unlist(lapply(names(parameters), function(name) {
          p <- NULL
          if(!is.null(parameters[[name]]$hierarchical)) {
            p <- parameters[[name]]$hierarchical$prior
          } else if(!is.null(parameters[[name]]$prior)) {
            items <- parameters[[name]]$prior
            if(!is.list(items)) items <- list(items)
            p <- unlist(lapply(items, function(prior) if(is.null(prior)) NULL else if(is.character(prior)) prior else if(length(prior) == 2) sprintf("%s ~ %s;", name, as.character(as.expression(prior[[2]]))) else if(length(prior) == 3L) as.character(as.expression(prior))))
          }
          if(length(p) == 0L) sprintf("// no prior for %s", name)
          else p
        }))
      )
    } else if(is.list(priors)) {
      c(
        "// user-defined priors",
        vapply(priors, function(prior) {
          sprintf("%s;", as.character(as.expression(prior)))
        }, character(1))
      )
    }
  )

  if(isTRUE(parallel)) {
    add_code(
      "transformed data",
      datmap()
    )
  }

  if(isTRUE(parallel)) {
    add_code(
      "transformed parameters",
      parmap(l_pars)
    )
  }


  add_code(
    "model",
    "// likelihood (only if prior != 0)",
    "if(target() != negative_infinity()) {",
    if(isTRUE(parallel)) {
      paste0("\ttarget += map_rect(log_lik_rect, phi, theta, x_r, x_i);")
    } else {
      paste0("\tfor(i in 1:N) target += log_lik_single(",paste(c(datsig(names_back = l_data, types = FALSE, index = "i"), parsig(l_pars, types = FALSE, index = "i")),collapse=", "),");")
    },
    "}"
  )

  if(isTRUE(save_log_lik)) {
    add_code(
      "generated quantities",
      "// likelihood",
      "vector[N] log_lik;",
      if(isTRUE(parallel)) {
        paste0("\tlog_lik = map_rect(log_lik_rect, phi, theta, x_r, x_i);")
      } else {
        paste0("for(i in 1:N) log_lik[i] = log_lik_single(",paste(c(datsig(names_back = l_data, types = FALSE, index = "i"), parsig(l_pars, types = FALSE, index = "i")),collapse=", "),");")
      }
    )
  }

  if(isTRUE(sanity_checks)) {
    add_code(
      "transformed data",
      "for(i in 1:N) {",
      "\tif(nS[i] < 1) reject(\"Inconsistency detected: According to the data, trial \",i,\" did not display any items!\");",
      sprintf("\tfor(j in 1:%d) {", locations),
      "\t\tif(R[i,j] && !S[i,j]) reject(\"Inconsistency detected: According to the data, the item at location #\",j,\", of trial \",i,\", was reported but there was no item displayed at that location!\");",
      if(!is.null(mydata$D)) "\t\tif(D[i,j] && !S[i,j]) reject(\"Inconsistency detected: According to the data, there should be a distractor at location #\",j,\", of trial \",i,\", but there was no item displayed at that location!\");",
      if(!is.null(mydata$D)) "\t\tif(R[i,j] && D[i,j]) reject(\"Inconsistency detected: According to the data, the item at location #\",j,\", of trial \",i,\", was reported as a target but it was a distractor!\");",
      "\t}",
      "}"
    )
  }

  header <- c(
    "StanTVA",
    "=======",
    "This is a StanTVA program, generated with RStanTVA. Please cite as:",
    "",
    strsplit(format(citation("RStanTVA"), style="text"),"\n")[[1]],
    "",
    "Configuration",
    "=============",
    sprintf(" - %s = %s", names(call_args), call_args),
    "",
    "License",
    "=======",
    "StanTVA and RStanTVA are licensed under the GNU General Public License 3. For a copy of the license agreement, see: https://www.gnu.org/licenses/gpl-3.0.html"
  ) %>% ansi_strwrap(width = 80L)


  ret <- paste(
    c(
      sprintf("/*%s*", strrep("*",4+max(nchar(header)))),
      #sprintf(" *  %s%s  *", header, strrep(" ",max(nchar(header))-nchar(header))),
      paste0(" *  ", header),
      sprintf(" *%s*/", strrep("*",4+max(nchar(header)))),
      "",
      vapply(names(code_blocks)[vapply(code_blocks, length, integer(1)) > 0L], function(block_name) paste(sprintf("%s {", block_name), paste("\t", code_blocks[[block_name]], collapse ="\n", sep=""), "}", sep = "\n"), character(1)),
      ""
    ),
    collapse="\n"
  )


  df <- vapply(names(parameters), function(pn) if(grepl("^simplex\\b", parameters[[pn]]$type)) as.integer(parameters[[pn]]$dim-1L) else if(is.null(parameters[[pn]]$dim)) 1L else as.integer(parameters[[pn]]$dim), integer(1), USE.NAMES = TRUE)
  if(has_formulas) {
    attr(df, "formula_lhs") <- formula_lhs
    attr(df, "random_factors") <- all_random_factors
  }
  new("stantvacode", code = ret, config = call_args_list, include_path = stantva_path(), df = df)
}

#'@export
stantvacode <- setClass("stantvacode", slots = c("code" = "character", "config" = "list", "include_path" = "character", "df" = "integer"))


#'@export
setMethod("show", "stantvacode", function(object) {
  cat(col_grey("// Include path(s): ", paste0(object@include_path, collapse="; ")),"\n")
  cat(object@code)
})


#'@export
stantva_model <- function(..., stan_options = list()) {
  args <- list(...)
  mc <- if(length(args) == 1 && inherits(args[[1]], "stantvacode")) args[[1]] else do.call(stantva_code, args)
  stan_options$model_code <- mc@code
  stan_options$isystem <- c(mc@include_path, stan_options$isystem)
  if(isTRUE(mc@config$parallel) && rstan_options("threads_per_chain") <= 1L) {
    stop("You requested a parallel model but RStan has not been configured to use multithreading! Try `rstan_options(threads_per_chain = ...)` to set the appropriate number of parallel threads within each chain before compiling the model code! To use all available CPUs, try `rstan_options(threads_per_chain = parallel::detectCores())`. If you do not wish to use multithreading, regenerate the model with `parallel = FALSE`!")
  }
  m <- do.call(stan_model, stan_options) %>% as("stantvamodel")
  m@code <- mc
  m
}


#'@export
write_stantva_model <- function(model, file = stdout()) {
  code <- if(inherits(model, "stantvamodel") || inherits(model, "stantvafit")) {
    model_code(model)
  } else if(inherits(model, "stantvacode")) {
    model
  } else {
    stop("`model` must be of type stantvamodel, stantvafit or stantvacode!")
  }
  writeLines(code@code, file)
}




#'@export
tvadata <- setClass("tvadata", contains = "tbl_df")

#'@export
setMethod("show", "tvadata", function(object) {
  if(is.null(object$D)) {
    cat(col_cyan("TVA"), "data containing",nrow(object),"whole-report trial(s)\n")
  } else {
    cat(col_cyan("TVA"), "data containing",nrow(object),"whole- and/or partial-report trial(s)\n")
  }
  callNextMethod()
})





#'@importClassesFrom rstan stanmodel
#'@export
stantvamodel <- setClass("stantvamodel", contains = "stanmodel", slots = c("code" = "stantvacode"))

#'@importClassesFrom rstan stanfit
#'@export
stantvafit <- setClass("stantvafit", contains = "stanfit", slots = c("stanmodel" = "stantvamodel", "data" = "list"))



#'@export
setMethod("show", c(object="stantvamodel"), function(object) {
  cat(col_cyan("StanTVA"), "model with", length(object@code@df),"free parameter(s) and the following configuration:\n")
  for(cname in names(object@code@config)) {
    cat(ansi_strwrap(paste0("- ",col_magenta(cname)," = ",deparse1(object@code@config[[cname]])), indent = 2, exdent = 6),sep="\n")
  }
})


setGeneric("generate", function(x, ...) {})

#'@export
setMethod("generate", c(x="stantvamodel"), function(x, data, params, vars, seed = NULL) {
  s <- gqs(x, data, params, seed = if(is.null(seed)) sample.int(.Machine$integer.max, size = 1L) else seed)
  if(missing(vars)) extract(s)
  else extract(s, vars)
})

#'@export
setMethod("generate", "stantvafit", function(x, newdata, vars, seed = NULL) {
  if(missing(newdata) || is.null(newdata)) {
    generate(x@stanmodel, x@data, as.matrix(x), vars, seed)
  } else {
    generate(x@stanmodel, newdata, as.matrix(x), vars, seed)
  }
})


setGeneric("fit", function(object, ...) {})

init_sampler <- function(model, pdata) {
  function(num_chain = 1) {
    if(is.null(pdata$X)) {
      list()
    } else {
      ret <- list(b = double(ncol(pdata$X)))
      ret$b[colnames(pdata$X) == "Intercept"] <- as.array(runif(sum(colnames(pdata$X) == "Intercept"), -1, 1))
      rfs <- bind_rows(bind_rows(lapply(model@code@config$formula, parse_formula))$random)
      for(i in seq_len(nrow(rfs))) {
        for(j in seq_len(model@code@df[rfs$param[i]])) {
          rf <- if(model@code@df[rfs$param[i]] > 1L && !rfs$custom[i]) paste0(rfs$group[i], "_", j) else rfs$group[i]
          ret[[paste0("r_",rf)]] <- diag(ncol(pdata[[paste0("Z_",rf)]]))
          ret[[paste0("s_",rf)]] <- as.array(if_else(colnames(pdata[[paste0("Z_",rf)]]) == "Intercept", 0.1, 0.01))
          ret[[paste0("w_",rf)]] <- matrix(0, nrow = pdata[[paste0("N_", rfs$factor_txt[i])]], ncol = ncol(pdata[[paste0("Z_",rf)]]))
        }
      }
      ret
    }
  }
}

#'@export
setMethod("sampling", c(object = "stantvamodel"), function(object, data, pars = NULL, include = TRUE, chains = 4, init, ...) {
  if(object@code@config$locations != ncol(data$S)) stop("Cannot fit a StanTVA model compiled for ",object@code@config$locations," location(s) to a data set with ",ncol(data$S)," location(s)!")
  pdata <- prepare_data(data, object)
  formula_lhs <- attr(object@code@df, "formula_lhs")
  pars_to_exclude <- character()
  if(isTRUE(object@code@config$parallel)) {
    pars_to_exclude <- union(pars_to_exclude, c("theta","phi"))
  }
  if(!is.null(formula_lhs) && ncol(formula_lhs) > 0) {
    pars_to_exclude <- union(pars_to_exclude, formula_lhs[2,])
  }
  if(length(pars) == 0 || (length(pars) == 1 && is.na(pars) && isTRUE(include))) {
    include <- FALSE
    pars <- pars_to_exclude
  } else if(isFALSE(include)) {
    pars <- union(pars_to_exclude, na.omit(pars))
  } else {
    pars <- setdiff(na.omit(pars), pars_to_exclude)
  }
  if(missing(init) && !is.null(pdata$X)) {
    init <- init_sampler(object, pdata)
  } else if(missing(init)) {
    init <- "random"
  }
  if(length(pars) == 0 && isFALSE(include)) {
    f <- callNextMethod(object, pdata, pars = NA, include = TRUE, init = init, chains = chains, ...)
  } else {
    f <- callNextMethod(object, pdata, pars = pars, include = include, init = init, chains = chains, ...)
  }
  f@stanmodel <- object
  f <- as(f, "stantvafit")
  f@data <- pdata
  f
})

#'@export
setMethod("optimizing", c(object = "stantvamodel"), function(object, data, pars = NULL, include = TRUE, init, ...) {
  if(object@code@config$locations != ncol(data$S)) stop("Cannot fit a StanTVA model compiled for ",object@code@config$locations," location(s) to a data set with ",ncol(data$S)," location(s)!")
  formula_lhs <- attr(object@code@df, "formula_lhs")
  pdata <- prepare_data(data, object)
  r <- callNextMethod(object, pdata, init = if(missing(init)) init_sampler(object, pdata) else init)
  if(!is.null(r$par)) {
    original_names <- names(r$par)
    keep <- !grepl("^(theta|phi)\\[|^r_", original_names)
    r$par <- r$par[keep]
    r$theta_tilde <- r$theta_tilde[,keep,drop=FALSE]
    p <- translate_names(object, pdata, original_names[keep])
    names(r$par) <- attr(p, "alias")
    colnames(r$theta_tilde) <- attr(p, "alias")
  }
  r
})


#'@export
setMethod("fit", c(object="stantvamodel"), function(object, data, method = c("optimizing","sampling"), ...) {
  method <- match.arg(method)
  do.call(method, list(object = object, data = data, ...))
})

#'@export
setMethod("logLik", "stantvamodel", function(object, data, params) {
  if(!isTRUE(object@code@config$save_log_lik)) stop("StanTVA model must be compiled with `save_log_lik` = TRUE in order to use logLik()!")
  generate(object, data, params, "log_lik")$log_lik
})

#'@export
setMethod("logLik", "stantvafit", function(object, newdata) {
  if(!isTRUE(object@stanmodel@code@config$save_log_lik)) stop("StanTVA model must be compiled with `save_log_lik` = TRUE in order to use logLik()!")
  if(missing(newdata)) extract(object, "log_lik")$log_lik
  else {
    logLik(object@stanmodel, newdata, extract(object))
  }
})


#' @export
list2stantvafit <- function(fits) {
  r <- sflist2stanfit(fits) %>% as("stantvafit")
  r@data <- fits[[1]]@data
  r
}

alias.stantvafit <- function(object) {
  attr(names(object), "alias")
}

#'@export
setMethod("alias", "stantvafit", alias.stantvafit)

fixef.stantvafit <- function(object) {
  formula_lhs <- attr(object@stanmodel@code@df, "formula_lhs")
  r <- NULL
  if(!is.null(formula_lhs)) {
    b <- extract(object, "b")$b
    colnames(b) <- alias(object)[match(sprintf("b[%d]", seq_len(ncol(b))), names(object))]
    r <- rbind(r, b)
  }
  for(p in setdiff(names(object@stanmodel@code@df), formula_lhs[2,])) {
    b <- extract(object, p)[[1]]
    colnames(b) <- p
    r <- rbind(r, b)
  }
  r
}

#'@export
setMethod("fixef", "stantvafit", fixef.stantvafit)


ranef.stantvafit <- function(object) {
  all_params <- bind_rows(lapply(object@stanmodel@code@config$formula, parse_formula)) %>% rename(name = param) %>% mutate(dim = vapply(name, function(x) if(!is.null(parameters[[x]]$dim)) as.integer(parameters[[x]]$dim) else 1L, integer(1)), fdim = vapply(name, function(x) if(grepl("^simplex", parameters[[x]]$type)) as.integer(parameters[[x]]$dim)-1L else if(!is.null(parameters[[x]]$dim)) as.integer(parameters[[x]]$dim) else 1L, integer(1)))
  for(i in seq_len(nrow(all_params))) {
    if(all_params$dim[i] > 1L) {
      is_custom <- all_params$random[[i]]$custom
      all_params$random[[i]] <- crossing(all_params$random[[i]], index = seq_len(all_params$fdim[i]))
      if(!is_custom) {
        all_params$random[[i]]$group <- sprintf("%s_%d", all_params$random[[i]]$group, all_params$random[[i]]$index)
      }
    } else {
      all_params$random[[i]]$index <- NA_integer_
    }
  }
  g <- bind_rows(all_params$random)
  if(nrow(g) == 0L) return(list())
  sapply(unique(g$factor_txt), function(rf) {
    gs <- unique(filter(g, factor_txt == rf)$group)
    bos <- rstan::extract(object, paste0("w_", gs))
    b <- array(NA_real_, dim = c(dim(bos[[1]])[1:2], sum(vapply(bos, function(x) dim(x)[[3]], integer(1)))))
    dimnames(b)[[3]] <- character(dim(b)[[3]])
    y <- 1L
    for(o in gs) {
      bo <- bos[[paste0("w_",o)]]
      dimnames(bo)[[3]] <- character(dim(bo)[[3]])
      for(i in which(g$group == o)) {
        if(is.na(g$index[i])) {
          m <- object@data[[sprintf("map_%s_%s", g$param[i], o)]]
          dimnames(bo)[[3]][m] <- sprintf("%s_%s", g$param[i], colnames(object@data[[sprintf("Z_%s", o)]])[m])
        } else  {
          m <- object@data[[sprintf("map_%s_%d_%s", g$param[i], g$index[i], o)]]
          dimnames(bo)[[3]][m[,j]] <- sprintf("%s_%s[%d]", g$param[i], colnames(object@data[[sprintf("Z_%s", o)]])[m[,j]], g$index[i])
        }
      }
      b[,,y:(y+dim(bo)[[3]]-1L)] <- bo
      dimnames(b)[[3]][y:(y+dim(bo)[[3]]-1L)] <- dimnames(bo)[[3]]
      y <- y + dim(bo)[[3]]
    }
    b
  }, simplify = FALSE)
}

#'@export
setMethod("ranef", "stantvafit", ranef.stantvafit)


coef.stanfit <- function(object) {
  fixefs <- fixef(object)
  ranefs <- ranef(object)
  ret <- list()
  for(i in names(ranefs)) {
    ranefs_with_fixefs <- intersect(colnames(fixefs), dimnames(ranefs[[i]])[[3]])
    fixefs_without_ranefs <- setdiff(colnames(fixefs), dimnames(ranefs[[i]])[[3]])
    ret[[i]] <- array(NA_real_, dim = c(dim(ranefs[[i]])[1:2], length(ranefs_with_fixefs)+length(fixefs_without_ranefs)), dimnames = list(NULL,NULL,c(fixefs_without_ranefs,ranefs_with_fixefs)))
    for(j in dim(ranefs[[i]])[[2]]) {
      ret[[i]][,j,ranefs_with_fixefs] <- ranefs[[i]][,j,ranefs_with_fixefs] + fixefs[,ranefs_with_fixefs]
      if(length(fixefs_without_ranefs) > 0L) {
        ret[[i]][,j,fixefs_without_ranefs] <- fixefs[,fixefs_without_ranefs]
      }
    }
  }
  ranefs
}

#'@export
setMethod("coef", "stantvafit", coef.stanfit)

predict.stantvafit <- function(object, newdata, variables = names(object@stanmodel@code@df)) {
  p <- extract(object)

  newdata <- if(missing(newdata)) object@data else prepare_data(newdata, object@stanmodel, FALSE)
  fx <- bind_rows(lapply(object@stanmodel@code@config$formula, parse_formula))
  sapply(variables, function(parname) {
    which_formula <- match(parname, fx$param)
    if(is.na(which_formula)) {
      p[[parname]]
    } else {

      par_dim <- if(length(object@par_dims[[parname]]) > 1L) object@par_dims[[parname]][2L] else 1L
      par_df <- object@stanmodel@code@df[parname]
      r <- vapply(seq_len(par_dim), function(i) {
        m <- newdata[[if(par_dim > 1) paste0("map_",parname,"_",i) else paste0("map_",parname)]]
        if(i > par_df) return(matrix(1, length(p$lp__), newdata$N))
        y <- tcrossprod(newdata$X[,m], p$b[,m])
        bt <- fx$inverse_link[[which_formula]]
        rfs <- fx$random[[which_formula]]
        for(k in seq_len(nrow(rfs))) {
          mrf <- newdata[[if(par_dim > 1) paste0("map_",parname,"_",i,"_",rfs$group[k],"_",i) else paste0("map_",parname,"_",rfs$group[k])]]
          for(j in seq_len(newdata[[paste0("N_",rfs$factor_txt[k])]])) {
            J <- newdata[[rfs$factor_txt[k]]] == j
            #print(which(J))
            #print(if(par_dim > 1) paste0("w_",rfs$group[k],"_",i) else paste0("w_",rfs$group[k]))
            Z <- newdata[[if(par_dim > 1) paste0("Z_",rfs$group[k],"_",i) else paste0("Z_",rfs$group[k])]][J,mrf,drop=FALSE]
            w <- t(as.matrix(p[[if(par_dim > 1) paste0("w_",rfs$group[k],"_",i) else paste0("w_",rfs$group[k])]][,j,mrf]))
            W <- Z %*% w
            #print(if(par_dim > 1) paste0("Z_",rfs$group[k],"_",i) else paste0("Z_",rfs$group[k]))
            #print(names(newdata))
            #print(Z)
            y[J,] <- y[J,] + W
          }
        }
        t(fx$inverse_link[[which_formula]](y))
      }, matrix(NA_real_, length(p$lp__), newdata$N))
      if(par_dim > 1L && par_dim > par_df) {
        rs <- do.call(cbind, apply(r, 2, rowSums, simplify = FALSE))
        for(i in seq_len(par_dim)) {
          r[,,i] <- r[,,i] / rs
        }
      }
      if(object@stanmodel@code@df[parname] > 1) r else r[,,1]
    }
  }, simplify = FALSE)
}

#'@export
setMethod("predict", "stantvafit", predict.stantvafit)


fitted.stantvafit <- function(object) {
  predict(object)
}

#'@export
setMethod("fitted", "stantvafit", fitted.stantvafit)



#'@export
tva_report <- function(data) {
  tibble(
    condition = data$condition,
    exposure = data$T,
    score = as.integer(if(is.null(data$D)) rowSums(data$R == 1L & data$S == 1L) else rowSums(data$R == 1L & data$S == 1L & data$D == 0L)),
    n_items = as.integer(rowSums(data$S == 1L)),
    n_distractors = if(is.null(data$D)) integer(nrow(data)) else as.integer(rowSums(data$D))
  ) %>% mutate(n_targets = n_items - n_distractors)
}


#'@export
setMethod("show", "stantvafit", function(object) print(object))



#'@export
setMethod("print", "stantvafit", function(x, digits_summary = 2, ...) {

  par_names <- names(x)

  sampler <- attr(x@sim$samples[[1]], "args")$sampler_t

  cat(col_cyan("StanTVA"), "model with", length(x@stanmodel@code@df), "free parameter(s), fitted with ")
  cat(x@sim$chains," ",sampler, " chains, each with iter=", x@sim$iter,
      "; warmup=", x@sim$warmup, "; thin=", x@sim$thin,"\n", sep = "")


  heading <- function(txt) cat("\n", style_underline(style_bold(txt)), "\n", sep = "")

  heading("Model configuration:")
  cat(sprintf("%s = %s\n", names(x@stanmodel@code@config), vapply(x@stanmodel@code@config, function(x) deparse1(x), character(1))), sep = "")


  fx <- bind_rows(lapply(x@stanmodel@code@config$formula, parse_formula))


  global_pars <- setdiff(names(x@stanmodel@code@df), fx$param)

  not_converged <- character()

  if(length(global_pars) > 0) {

    heading("Global parameters:")

    g_summary <- rstan::summary(x, global_pars)$summary

    print(round(g_summary, digits_summary), max = prod(dim(g_summary)))

    not_converged <- c(not_converged, rownames(g_summary)[g_summary[,"Rhat"] >= 1.05])

  }


  if(nrow(fx) > 1) {

    rfs <- fx$random %>% lapply(function(fxi) if(x@stanmodel@code@df[fxi$param] > 1L) crossing(fxi, index = seq_len(x@stanmodel@code@df[fxi$param])) else bind_cols(fxi, index = NA_integer_)) %>% bind_rows() %>% mutate(group = if_else(custom | is.na(index), group, paste0(group,"_",index)))

    random_factors_txt <- if(nrow(rfs) > 0L) unique(rfs$factor_txt) else character()

    b_summary <- rstan::summary(x, "b")$summary

    rownames(b_summary) <- attr(par_names, "alias")[match(rownames(b_summary), par_names)]

    heading("Fixed effects:")

    print(round(b_summary, digits_summary), max = prod(dim(b_summary)))

    not_converged <- c(not_converged, rownames(b_summary)[b_summary[,"Rhat"] >= 1.05])


    for(rf in random_factors_txt) {

      heading(paste0("Hyperparameters on random effects (", col_blue(deparse1(rfs$factor[[match(rf, rfs$factor_txt)]])), " level, N = ", x@data[[paste0("N_",rf)]] ,"):"))



      gs <- unique(rfs$group[rfs$factor_txt == rf])



      s_summary <- rstan::summary(x, paste0("s_", gs))$summary

      rownames(s_summary) <- attr(par_names, "alias")[match(rownames(s_summary), par_names)]

      for(g in gs) {

        M_rf <- x@par_dims[[sprintf("r_%s", g)]][1]
        if(M_rf >= 2L) {
          c_pars <- combn(M_rf, 2L)
          r_summary <- rstan::summary(x, paste0("r_", g, "[", c_pars[1,], ",", c_pars[2,],"]"))$summary
          rownames(r_summary) <- attr(par_names, "alias")[match(rownames(r_summary), par_names)]
          s_summary <- rbind(s_summary, r_summary)
        }
      }

      print(round(s_summary, digits_summary), max = prod(dim(s_summary)))



      not_converged <- c(not_converged, rownames(s_summary)[s_summary[,"Rhat"] >= 1.05])

    }


    w_summary <- rstan::summary(f, sprintf("w_%s", gs), probs = double(0))$summary

    not_converged <- c(not_converged, attr(par_names, "alias")[match(rownames(w_summary)[w_summary[,"Rhat"] >= 1.05], par_names)])

  }



  if(length(not_converged)) {
    warning("Model did not converge (Rhat >= 1.05) for parameter(s) ", paste(not_converged, collapse=", "))
  }

})

translate_names <- function(model, data, names) {

  ret <- names
  new_ret <- names

  fx <- bind_rows(lapply(model@code@config$formula, parse_formula))

  if(nrow(fx) > 0L) {
    for(i in seq_len(nrow(fx))) {
      param <- fx$param[i]
      if(model@code@df[param] == 1L) {
        m <- data[[sprintf("map_%s", param)]]
        new_ret[match(sprintf("b[%d]", m), ret)] <- paste0(fx$param[i], "_", colnames(data$X)[m])
      } else for(j in seq_len(model@code@df[param])) {
        m <- data[[sprintf("map_%s_%d", param, j)]]
        new_ret[match(sprintf("b[%d]", m), ret)] <- paste0(fx$param[i], "_", colnames(data$X)[m], "[",j,"]")
      }
    }
    rfs <- if(is.null(fx$random) || nrow(bind_rows(fx$random)) == 0) tibble() else fx$random %>% lapply(function(fxi) if(!is.null(fxi$param) && model@code@df[fxi$param] > 1L) crossing(fxi, index = seq_len(model@code@df[fxi$param])) else bind_cols(fxi, index = NA_integer_)) %>% bind_rows() %>% mutate(group = if_else(custom | is.na(index), group, paste0(group,"_",index)))
    for(i in seq_len(nrow(rfs))) {
      param <- rfs$param[i]
      if(model@code@df[param] == 1L) {
        m <- data[[sprintf("map_%s_%s", rfs$param[i], rfs$group[i])]]
        new_ret[match(sprintf("s_%s[%d]", rfs$group[i], m), ret)] <- sprintf("sd(%s_%s_%s)", rfs$param[i], rfs$factor_txt[i], colnames(data[[paste0("Z_", rfs$group[i])]])[m])
        for(s in seq_len(data[[paste0("N_",rfs$factor_txt[i])]])) new_ret[match(sprintf("w_%s[%d,%d]", rfs$group[i], s, m), ret)] <- sprintf("%s_%s_%s[%d]", rfs$param[i], rfs$factor_txt[i], colnames(data[[paste0("Z_", rfs$group[i])]])[m], s)
      } else for(j in seq_len(model@code@df[param])) {
        m <- data[[sprintf("map_%s_%d_%s", rfs$param[i], j, rfs$group[i])]]
        new_ret[match(sprintf("s_%s[%d]", rfs$group[i], m), ret)] <- sprintf("sd(%s_%s_%s[%d])", rfs$param[i], rfs$factor_txt[i], colnames(data[[paste0("Z_", rfs$group[i])]])[m], j)
        for(s in seq_len(data[[paste0("N_",rfs$factor_txt[i])]])) new_ret[match(sprintf("w_%s[%d,%d]", rfs$group[i], s, m), ret)] <- sprintf("%s_%s_%s[%d,%d]", rfs$param[i], rfs$factor_txt[i], colnames(data[[paste0("Z_", rfs$group[i])]])[m], s, j)
      }
    }
    if(nrow(rfs) > 0L) for(g in unique(rfs$group)) {
      rf <- rfs$factor_txt[match(g, rfs$group)]
      pars <- unique(rfs$param[rfs$group==g])
      cnames <- colnames(data[[paste0("Z_", g)]])
      pnames <- character(length(cnames))
      for(i in seq_along(pars)) {
        param <- pars[i]
        if(model@code@df[param] == 1L) {
          m <- data[[sprintf("map_%s_%s", pars[i], g)]]
          pnames[as.vector(m)] <- param
        } else for(j in seq_len(model@code@df[param])) {
          m <- data[[sprintf("map_%s_%d_%s", pars[i], j, g)]]
          pnames[as.vector(m)] <- param
          cnames[m] <- sprintf("%s[%d]", cnames[m], j)
        }
      }
      for(a in seq_along(cnames)) {
        for(b in seq_along(cnames)) {
          new_ret[match(sprintf("r_%s[%d,%d]", g, a, b), ret)] <- sprintf("cor(%s_%s_%s,%s_%s_%s)", pnames[a], rf, cnames[a], pnames[b], rf, cnames[b])
        }
      }
    }

  }

  attr(ret, "alias") <- new_ret

  ret

}


#'@export
setMethod("names", "stantvafit", function(x) translate_names(x@stanmodel, x@data, x@sim$fnames_oi))


