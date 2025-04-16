#'@importFrom rstan extract stan_model sampling optimizing gqs sflist2stanfit rstan_options read_stan_csv
#'@importFrom dplyr summarize mutate group_by %>% across if_else select bind_cols bind_rows rename last filter transmute n
#'@importFrom tidyr pivot_longer pivot_wider crossing
#'@importFrom readr read_table write_tsv
#'@importFrom methods formalArgs new as callNextMethod show
#'@importFrom stats na.omit as.formula model.matrix pnorm runif terms contrasts qlogis plogis qnorm formula
#'@importFrom cli col_cyan col_magenta col_grey col_blue ansi_strwrap style_underline style_bold
#'@importFrom tibble tibble
#'@importFrom utils citation str combn packageName packageVersion
#'@importFrom lme4 findbars subbars fixef ranef nobars
#'@importFrom brms prior set_prior
#'@importFrom rlang .data .env
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

get_prior <- function(priors, class, dpar = NA_character_, group = NA_character_, coef = NA_character_) {
  if(is.null(priors)) return(NULL)
  m <- priors$class == class &
    (priors$dpar == "" | is.na(dpar) | dpar == priors$dpar) &
    (priors$group == "" | is.na(group) | group == priors$group) &
    (priors$coef == "" & coef != "Intercept" | is.na(coef) & priors$coef == "" | coef == priors$coef)
  if(any(m)) priors$prior[last(which(m))]
  else NULL
}

deparse_prior <- function(prior) {
  if(is.null(prior) || nrow(prior) == 0) {
    return(NULL)
  } else if(nrow(prior) == 1)  {
    args <- list(str2lang(prior$prior))
    if(prior$class != "b") args$class <- str2lang(prior$class)
    for(arg in c("coef","group","dpar")) if(prior[[arg]] != "") args[[arg]] <- str2lang(prior[[arg]])
    for(arg in c("lb","ub")) if(!is.na(prior[[arg]])) args[[arg]] <- prior[[arg]]
    return(as.call(c(list(as.name("prior")), args)))
  } else {
    return(call("+", sys.function()(prior[seq_len(nrow(prior)-1),]),sys.function()(prior[nrow(prior),])))
  }
}

#'@importFrom brms prior
#'@export
brms::prior


#' @title StanTVA path
#' @description Returns the path to the StanTVA directory.
#' @return The path to the StanTVA directory.
#' @details This function is used internally by the \code{\link[RStanTVA:stantva_model]{stantva_model()}} method.
#' @examples
#' path <- stantva_path()
#' path
#'@export
stantva_path <- function() {
  file.path(find.package(packageName()), "StanTVA")
}

#' @title Read TVA data
#' @description Reads TVA data from a file.
#' @param file The file name.
#' @param set The set of items.
#' @param ... Additional arguments passed to \code{\link[readr:read_table]{read_table()}}.
#' @return A TVA data object, which inherits from \code{data.frame}.
#' @examples
#' \donttest{
#' data <- read_tva_data("data.dat")
#' data
#' }
#'@export
read_tva_data <- function(file, set = LETTERS, ...) {
  if(inherits(file, "connection")) f <- file
  else f <- base::file(file, "rb")
  n <- as.integer(readLines(f, n = 1))
  dat <- read_table(f, col_names = c("condition", "exposure", "targets", "distractors", "report"), col_types = c("inccc"), n_max = n, na = character(), ...)
  if(n != nrow(dat)) warning("Expected ",n," trials but only read ",nrow(dat)," lines!")
  close(f)
  dat %>% mutate(
    across(c(.data$targets, .data$distractors), ~do.call(rbind, strsplit(.x, "", TRUE) %>% lapply(function(x) if_else(x=="0",NA_character_,x)))),
    report = strsplit(if_else(.data$report == "-", "", .data$report), "", TRUE),
    S = (!is.na(.data$targets) | !is.na(.data$distractors))+0L,
    D = (!is.na(.data$distractors) )+ 0L,
    R = t(vapply(seq_len(n()), function(i) .data$targets[i,] %in% .data$report[[i]] | .data$distractors[i,] %in% .data$report[[i]], logical(ncol(.data$targets)))) + 0L,
    items = t(vapply(seq_len(n()), function(i) if_else(is.na(.data$targets[i,]), .data$distractors[i,], .data$targets[i,]), character(ncol(.data$targets)))),
    E = vapply(seq_len(n()), function(i) sum(!.data$report[[i]] %in% .data$items[i,]), integer(1)),
    I = length(set)
  ) %>% select(.data$condition, .data$items, .data$report, .data$S, .data$D, T = .data$exposure, .data$R, .data$E, .data$I)
}

#' @title Write TVA data
#' @description Writes TVA data to a file.
#' @param data The TVA data object.
#' @param file The file name.
#' @param ... Additional arguments passed to \code{\link[readr:write_tsv]{write_tsv()}}.
#' @return No return value, called for side effects.
#' @examples
#' \donttest{
#' data <- read_tva_data("data.dat")
#' write_tva_data(data, "data.dat")
#' }
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

#' Extract Stan code
#' @export
setGeneric("model_code", function(object, type) UseMethod("model_code"))

#' @title Extract Stan code
#' @describeIn model_code method
#' @description Extracts the Stan code from a StanTVA model object.
#' @param object A StanTVA model object or fit.
#' @param type The type of code to return (\code{stan}: formatted StanTVA, \code{stan2}: ready-to-compile Stan code, \code{cpp}: generated C++ code).
#' @return A RStanTVA model code object (\code{stan}), or a string containing the code (\code{stan2} or \code{cpp}).
#' @examples
#' \donttest{
#' model <- stantva_model(locations = 2)
#' model_code(model)
#' }
#' @export
setMethod("model_code", "stanmodel", function(object, type = c("stan","stan2","cpp")) {
  type <- match.arg(type)
  if(type == "stan") object@code
  else if(type == "stan2") object@model_code
  else if(type == "cpp") object@model_cpp
  else stop("Unknown type '", type,"'!")
})

#'@describeIn model_code Extract code from a model fit
#'@importClassesFrom rstan stanfit
#'@export
setMethod("model_code", "stanfit", function(object, type) {
  model_code(object@stanmodel, type)
})


## add hierarchical stuff here

prepare_data <- function(trials, model, require_outcome = TRUE, contrasts = list()) {

  mc <- if(inherits(model, "stantvamodel")) model@code else if(inherits(model, "stantvacode")) model else stop("`model` must be a StanTVA model or StanTVA model code object!")


  for(cn in names(contrasts)) {
    if(cn %in% colnames(trials)) {
      if(is.numeric(trials[[cn]])) next
      if(!is.factor(trials[[cn]])) {
        trials[[cn]] <- as.factor(trials[[cn]])
      }
      x <- setdiff(levels(trials[[cn]]), rownames(contrasts[[cn]]))
      if(length(x) > 0) {
        stop("Originally fitted dataset did not contain level(s) ",paste0(x,collapse=", ")," for factor ",cn,"!")
      }
      trials[[cn]] <- factor(as.character(trials[[cn]]), levels = rownames(contrasts[[cn]]))
      stats::contrasts(trials[[cn]]) <- contrasts[[cn]]
    }
  }


  fs <- mc@config$formula

  required_columns <- if(isFALSE(require_outcome)) c() else if(mc@config$task == "wr") c("R","S","T") else if(mc@config$task == "pr") c("R","S","D","T") else stop("Unsupported task ",sQuote(mc@config$task),"!")
  if(isTRUE(mc@config$allow_guessing)) required_columns <- c(required_columns, "E", "I")


  if(any(!required_columns %in% colnames(trials))) {
    stop("The supplied data must contain columns ",paste(sQuote(required_columns), collapse=", ")," but at least one is missing: ",paste(sQuote(setdiff(required_columns, colnames(trials))), collapse=", "))
  }

  ltrials <- as.list(trials[,required_columns,drop=FALSE])

  pf <- bind_rows(lapply(fs, parse_formula))

  ltrials$N <- nrow(trials)


  for(i in seq_len(nrow(pf))) {
    for(j in seq_len(nrow(pf$random[[i]]))) {
      x <- eval(pf$random[[i]]$factor[[j]], trials)
      ltrials[[pf$random[[i]]$factor_txt[[j]]]] <- if(is.factor(x)) {
        rfc <- as.integer(x)
        attr(rfc, "levels") <- levels(x)
        rfc
      } else if(is.numeric(x)) {
        rfc <- x
        attr(rfc, "levels") <- seq_len(max(x))
        rfc
      } else {
        stop("Random factor ",sQuote(deparse1(pf$random[[i]]$factor[[j]]))," must be integer or factor!")
      }
      ltrials[[paste0("N_",pf$random[[i]]$factor_txt[[j]])]] <- length(attr(rfc, "levels"))
    }

    f_var <- pf$param[i]

    Cmatf <- model.matrix(pf$fixed_formula[[i]], trials)
    colnames(Cmatf) <- clean_name(colnames(Cmatf))

    ltrials[[paste0("M_",f_var)]] <- ncol(Cmatf)
    int_C <- which(attr(Cmatf, "assign") == 0L)
    ltrials[[paste0("int_",f_var)]] <- if(length(int_C) != 1L) 0L else int_C
    for(j in seq_len(mc@df[f_var])) {
      ltrials$X <- cbind(ltrials$X, Cmatf)
      ltrials[[if(mc@dim[f_var] == 1L) paste0("map_",f_var) else paste0("map_",f_var,"_",j)]] <- array((ncol(ltrials$X)-ncol(Cmatf)+1L):ncol(ltrials$X), dim = ncol(Cmatf))
    }

    for(k in seq_len(nrow(pf$random[[i]]))) {
      rf <- pf$random[[i]]$factor_txt[k]
      Cmatr <- model.matrix(pf$random[[i]]$formula[[k]], trials)
      colnames(Cmatr) <- clean_name(colnames(Cmatr))
      for(j in seq_len(mc@df[f_var])) {
        g <- if(mc@dim[f_var] == 1L || pf$random[[i]]$custom[k]) pf$random[[i]]$group[k] else paste0(pf$random[[i]]$group[k], "_", j)
        ltrials[[paste0("M_",f_var,"_",g)]] <- ncol(Cmatr)
        int_C <- which(attr(Cmatr, "assign") == 0L)
        ltrials[[paste0("int_",f_var,"_",g)]] <- if(length(int_C) != 1L) 0L else int_C
        ltrials[[paste0("Z_",g)]] <- cbind(ltrials[[paste0("Z_",g)]], Cmatr)
        ltrials[[if(mc@dim[f_var] == 1L) paste0("map_",f_var,"_",g) else paste0("map_",f_var,"_",j,"_",g)]] <- array((ncol(ltrials[[paste0("Z_",g)]])-ncol(Cmatr)+1L):ncol(ltrials[[paste0("Z_",g)]]), dim = ncol(Cmatr))
      }
    }
  }

  attr(ltrials, "contrasts") <- list()
  for(cn in colnames(trials)) if(is.factor(trials[[cn]]) && nlevels(trials[[cn]]) > 1) {
    attr(ltrials, "contrasts")[[cn]] <- stats::contrasts(trials[[cn]])
  }

  ltrials
}

clean_name <- function(str) gsub("(^_+)|(_+$)$","",gsub("[^a-zA-Z0-9]+","_",str))


nested_parameter <- function(param, transform, type = "real", dim = 1L, prior_fixed_intercept = NULL, prior_fixed_slope = ~normal(0, 1)) {
  ret <- list(
    target = param,
    is_simplex = grepl("^simplex\\b", type),
    is_vector = grepl("^(vector|simplex)\\b", type),
    dim = dim,
    type_constraint = gsub("^[^<]*(<[^>]+>)?.*$", "\\1", type),
    transform = transform
  )
  ret$fdim <- if(ret$is_simplex) dim - 1L else dim

  ret
}


parse_formula <- function(f) {
  stopifnot(is.call(f) && f[[1L]] == "~" && length(f) == 3L)
  lhs <- f[[2L]]
  if(is.symbol(lhs)) {
    parname <- as.character(lhs)
    link_name <- "identity"
    inverse_link <- link <- quote(function(.) .)
    stan_link <- "%s"
    stan_inverse_link <- "%s"
    args <- list()
  } else if(is.call(lhs) && length(lhs) >= 2L && is.symbol(lhs[[2L]])) {
    parname <- as.character(lhs[[2L]])
    link_name <- as.character(lhs[[1L]])
    args <- as.list(lhs[c(-1,-2)])
    if(link_name == "log") {
      link <- quote(log)
      stan_link <- "log(%s)"
      inverse_link <- quote(exp)
      stan_inverse_link <- "exp(%s)"
    } else if(link_name == "logit") {
      link <- quote(qlogis)
      stan_link <- "logit(%s)"
      inverse_link <- quote(plogis)
      stan_inverse_link <- "inv_logit(%s)"
    } else if(link_name == "probit") {
      link <- quote(qnorm)
      stan_link <- "inv_Phi(%s)"
      inverse_link <- quote(pnorm)
      stan_inverse_link <- "Phi(%s)"
    } else if(link_name == "scaled_logit" && length(args) == 2L) {
      link <- substitute(function(x) qlogis((x-a)/(b-a)), list(a = args[[1]], b = args[[2]]))
      stan_link <- sprintf("logit((%%s-%1$g)/(%2$g-%1$g))", args[[1]], args[[2]])
      inverse_link <- substitute(function(x) a+(b-a)*plogis(x), list(a = args[[1]], b = args[[2]]))
      stan_inverse_link <- sprintf("%1$g+(%2$g-%1$g)*inv_logit(%%s)", args[[1]], args[[2]])
    } else if(link_name == "scaled_logit" && length(args) == 1L) {
      link <- substitute(function(x) qlogis(x/a), list(a = args[[1]]))
      stan_link <- sprintf("logit(%%s/%g)", args[[1]])
      inverse_link <- substitute(function(x) a*plogis(x), list(a = args[[1]]))
      stan_inverse_link <- sprintf("%g*inv_logit(%%s)", args[[1]])
    } else if(link_name == "scaled_probit" && length(args) == 2L) {
      link <- substitute(function(x) qnorm((x-a)/(b-a)), list(a = args[[1]], b = args[[2]]))
      stan_link <- sprintf("inv_Phi((%%s-%1$g)/(%2$g-%1$g))", args[[1]], args[[2]])
      inverse_link <- substitute(function(x) a+(b-a)*pnorm(x), list(a = args[[1]], b = args[[2]]))
      stan_inverse_link <- sprintf("%1$g+(%2$g-%1$g)*Phi(%%s)", args[[1]], args[[2]])
    } else if(link_name == "scaled_probit" && length(args) == 1L) {
      link <- substitute(function(x) qnorm(x/a), list(a = args[[1]]))
      stan_link <- sprintf("inv_Phi(%%s/%g)", args[[1]])
      inverse_link <- substitute(function(x) a*pnorm(x), list(a = args[[1]]))
      stan_inverse_link <- sprintf("%g*Phi(%%s)", args[[1]])
    } else stop("Unknown link ",sQuote(link_name),"!")
  } else stop("lhs must be a call or symbol!")
  fxs <- findbars(f[[3]])
  tibble(
    param = parname,
    link_name = link_name,
    link = list(eval(link)),
    stan_link = stan_link,
    inverse_link = list(eval(inverse_link)),
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

#' @title Generate StanTVA code
#' @description Creates a StanTVA model code object.
#' @param formula Optional formulas for nested and hierarchical model parameters.
#' @param locations The number of display locations (items).
#' @param task The task. Currently implemented: \code{wr} (whole report) and \code{pr} (partial report)
#' @param regions An optional list of groups of display locations (regions).
#' @param C_mode The mode/family for the $C$ parameter.
#' @param w_mode The mode/family for the $w$ parameter.
#' @param t0_mode The mode/family for the $t_0$ parameter.
#' @param K_mode The mode for the $K$ parameter.
#' @param max_K The upper bound of $K$.
#' @param allow_guessing (logical) Whether to allow guessing.
#' @param parallel (logical) Whether to use parallel chains.
#' @param save_log_lik (logical) Whether to save the log likelihood (needed for likelihood-based model comparison such as loo).
#' @param priors The priors.
#' @param sanity_checks (logical) Whether to perform sanity checks.
#' @param debug_neginf_loglik (logical) Whether to debug negative infinity log likelihood.
#' @return The StanTVA model code object.
#' @examples
#' model <- stantva_code(locations = 4, task = "pr")
#' model
#'@export
stantva_code <- function(formula = NULL, locations, task = c("wr","pr"), regions = list(), C_mode = c("equal","locations","regions"), w_mode = c("locations","regions","equal"), t0_mode = c("constant", "gaussian", "exponential", "shifted_exponential"), K_mode = c("bernoulli", "free", "binomial", "hypergeometric"), max_K = locations, allow_guessing = FALSE, parallel = isTRUE(rstan_options("threads_per_chain") > 1L), save_log_lik = FALSE, priors = NULL, sanity_checks = TRUE, debug_neginf_loglik = FALSE) {

  task <- match.arg(task)
  C_mode <- match.arg(C_mode)
  t0_mode <- match.arg(t0_mode)
  K_mode <- match.arg(K_mode)
  w_mode <- match.arg(w_mode)

  hierarchical_config <- bind_rows(lapply(formula, parse_formula))

  includeFile <- function(f) sprintf("#include %s", f)

  call_args <- character(0)
  call_args_list <- list()

  for(x in formalArgs(sys.function())) {
    if(x %in% c("data","type")) next
    val <- get(x)
    if(x == "priors") val <- deparse_prior(val)
    call_args[x] <- if(is.numeric(val) || is.character(val) || is.logical(val)) as.character(val) else deparse1(val)
    call_args_list[[x]] <- val
  }


  code_blocks <- list(
    `functions` = character(0),
    `data` = character(0),
    `transformed data` = character(0),
    `parameters` = character(0),
    `transformed parameters` = character(0),
    `transformed parameters 2` = character(0),
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
        hierarchical_config$stan_link[match(name, hierarchical_config$param)],
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
    add_param(name = "C", type = sprintf("vector<lower=machine_precision()>[%d]", length(regions)), ctype=sprintf("vector[%d]", length(regions)), rtype = "vector", dim = length(regions), prior = substitute(~gamma(a, 0.035), list(a=3.5/length(regions))))
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
    add_param(name = "r", type = sprintf("simplex[%d]", length(regions)), ctype=sprintf("vector[%d]", length(regions)), rtype = "vector", dim = length(regions))
    w_pars <- "r"
    w_body <- c(
      sprintf("vector[%1$d] w;", locations),
      unlist(lapply(seq_along(regions), function(i) {
        sprintf("w[%d] = r[%d]/%d.0;", regions[[i]], i, length(regions[[i]]))
      }))
    )
    #add_code(
    #  "generated quantities",
    #  paste0("real w_", names(regions)," = b[", seq_along(regions),"];")
    #)
  } else if(w_mode == "locations") {
    w_pars <- "w"
    w_body <- NULL
    add_param(name = "w", type = sprintf("simplex[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype ="vector", dim = locations, prior = ~lognormal(0,0.5))
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
    add_param(name = "K", class = c("phi", "K"), type = sprintf("real<lower=0,upper=%d>", max_K), ctype="real", rtype="real", prior = ~lognormal(1.1,0.3))
    K_args <- "[K]'"
  } else if(K_mode == "free") {
    add_code("functions", includeFile("freeK.stan"))
    add_param(name = "pK", class = c("phi","K"), type = sprintf("simplex[%d]", max_K+1L), ctype=sprintf("vector[%d]", max_K+1L), rtype="vector", dim = max_K+1L, prior = ~lognormal(0,0.5))
    #add_code("generated quantities", paste0("real mK = ",paste(sprintf("%d * pK[%d]", seq_len(locations), seq_len(locations)+1L), collapse=" + "),";"));
    K_args <- "pK"
  } else if(K_mode == "binomial") {
    add_code("functions", includeFile("binomialK.stan"))
    # TODO add prior!
    add_param(name = "nK", class = c("phi", "K"), type = "real<lower=machine_precision()>", ctype="real", rtype="real")
    add_param(name = "pK", class = c("phi", "K"), type = "real<lower=machine_precision(),upper=1.0-machine_precision()>", ctype="real", rtype="real", prior = ~beta(2,2))
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
    add_param(name = "t0", class = c("phi"), type = "real<upper=max_mu>", ctype="real", rtype="real", prior = ~normal(20, 15))
  } else if(t0_mode == "gaussian") {
    add_code("functions", includeFile("gaussiant0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real", ctype="real", rtype="real", prior = ~normal(20, 15))
    add_param(name = "sigma0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real", prior = ~gamma(2,0.2))
    t0_args <- "[mu0, sigma0]'"
  } else if(t0_mode == "exponential") {
    # TODO implement default priors!
    add_code("functions", includeFile("exponentialt0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real", prior = ~normal(20,15))
    #add_param(name = "t0", class = c("phi"), type = "real", ctype="real", rtype="real")
    t0_args <- "[1/mu0]'"
  }  else if(t0_mode == "shifted_exponential") {
    # TODO implement default priors!
    add_code("functions", includeFile("exponentialt0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real", prior = ~normal(10,10))
    add_param(name = "t0", class = c("phi"), type = "real", ctype="real", rtype="real", prior = ~normal(10,10))
    t0_args <- "[1/mu0]'"
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

  if(allow_guessing) {
    add_data(name = "E", class="x_i", type = "array[N] int<lower=0>", ctype = "int", rtype="int")
    add_data(name = "I", class="x_i", type = "array[N] int<lower=0>", ctype = "int", rtype="int")
    add_param(name = "g", class = "phi", type = "real<lower=machine_precision(),upper=1.0-machine_precision()>", ctype="real", rtype="real", prior = ~beta(2,20))
  }


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
    if(allow_guessing) {
      l_pars <- c(l_pars, "g")
      l_data <- c(l_data, "E", "I")
      l_body <- c(
        sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
        sprintf("log_lik = tva_wrg_log(R, S, %s, %s, %s, v, g, E, I);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
      )
    } else {
      l_body <- c(
        sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
        sprintf("log_lik = tva_wr_log(R, S, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
      )
    }
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
    if(allow_guessing) {
      l_pars <- c(l_pars, "g")
      l_data <- c(l_data, "E", "I")
      l_body <- c(
        sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
        sprintf("log_lik = tva_prg_log(R, S, D, %s, %s, %s, v, g, E, I);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
      )
    } else {
      l_body <- c(
        sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
        sprintf("log_lik = tva_pr_log(R, S, D, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
      )
    }
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

  if(isTRUE(debug_neginf_loglik)) {
    l_data <- c(l_data, "trial_no")
    add_data(name = "trial_no", class="x_i", type = "array[N] int", ctype="int", rtype="int", dim = 1, transformed = TRUE)
    add_code(
      "transformed data",
      "trial_no = linspaced_int_array(N, 1, N);"
    )
  }

  add_code(
    "functions",
    paste0("vector calculate_v(",paste(c(datsig(names_back = v_data), parsig(v_pars)), collapse = ", "),") {"),
    paste0("\t", v_body),
    "}"
  )

  # default parameter-unspecific priors
  default_priors <- prior("normal(0.0,0.05)", "sd") +
    prior("normal(0.0,0.1)", "sd", coef = "Intercept") +
    prior("lkj_corr(3.0)", "cor") +
    prior("normal(0.0,5.0)") +
    prior("normal(0.0,10.0)", coef = "Intercept")

  # default parameter-specific priors
  for(name in names(parameters)) {
    if(!is.null(parameters[[name]]$prior)) default_priors <- default_priors + set_prior(deparse1(rhs(parameters[[name]]$prior)), class = "global", dpar = name)
  }

  priors <- default_priors + priors

  # hierarchical stuff




  if(nrow(hierarchical_config) > 0L) {
    all_params <- hierarchical_config %>% rename(name = .data$param) %>% mutate(dim = vapply(name, function(x) if(!is.null(parameters[[x]]$dim)) as.integer(parameters[[x]]$dim) else 1L, integer(1)), fdim = vapply(name, function(x) if(grepl("^simplex", parameters[[x]]$type)) as.integer(parameters[[x]]$dim)-1L else if(!is.null(parameters[[x]]$dim)) as.integer(parameters[[x]]$dim) else 1L, integer(1)))
    for(i in seq_len(nrow(all_params))) {
      if(all_params$dim[i] > 1L) {
        is_custom <- !is.null(all_params$random[[i]]$custom) && all_params$random[[i]]$custom
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
    if(isTRUE(debug_neginf_loglik)) {
      l_data <- c(l_data, clean_name(all_random_factors))
    }
    all_random_params <- all_random_effects %>% group_by(.data$group, .data$factor_txt) %>% summarize(dim = sum(all_params$dim[match(.data$param, all_params$name)]), fdim = sum(all_params$fdim[match(.data$param, all_params$name)]), M_var = paste(sprintf("M_%s_%s", .data$param, .data$group), collapse = "+"))
    M_var <- paste(if_else(all_params$fdim!=1L,sprintf("%d*M_%s", all_params$fdim, all_params$name),sprintf("M_%s", all_params$name)), collapse = "+")
    for(rf in clean_name(all_random_factors)) {
      add_data(name = sprintf("N_%s", rf), type = "int<lower=1,upper=N>", ctype="int", rtype="int", dim = 1)
      add_data(name = rf, class="x_i", type = sprintf("array[N] int<lower=1,upper=N_%1$s>", rf), ctype="int", rtype="int", dim = 1)
    }
    add_code(
      "data",
      sprintf("int<lower=0> M_%s;", all_params$name),
      unique(sprintf("int<lower=0> M_%s_%s;", all_random_effects$param, all_random_effects$group)),
      unique(sprintf("int<lower=0,upper=M_%1$s_%2$s> int_%1$s_%2$s;", all_random_effects$param, all_random_effects$group)),
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
          sprintf("vector%2$s[N] %1$s;", all_params$name[i], gsub("^[^<>]*(<.*>)?[^<>]*$","\\1",parameters[[all_params$name[i]]]$type))
        } else {
          sprintf("matrix%3$s[N,%2$d] %1$s;", all_params$name[i], all_params$dim[i], gsub("^[^<>]*(<.*>)?[^<>]*$","\\1",parameters[[all_params$name[i]]]$type))
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
      "}",
      unlist(lapply(seq_len(nrow(all_params)), function(i) {
        if(all_params$dim[i] == 1) {
          sprintf("for(i in 1:N) if(is_nan(%1$s[i]) || is_inf(%1$s[i])) reject(\"Rejecting proposal because %1$s[\",i,\"] = \",%1$s[i],\" !\");", all_params$name[i])
        } else {
          sprintf("for(i in 1:N) for(j in 1:cols(%1$s)) if(is_nan(%1$s[i,j]) || is_inf(%1$s[i,j])) reject(\"Rejecting proposal because %1$s[\",i,\",\",j,\"] = \",%1$s[i,j],\" !\");", all_params$name[i])
        }
      }))
    )
    add_code(
      "transformed data",
      sprintf("vector[M_%1$s] mu_w_%1$s = rep_vector(0.0, M_%1$s);", all_random_params$group)
    )
    add_code(
      "model",
      unlist(lapply(seq_len(nrow(all_random_effects)), function(i) {
        prior_slope <- get_prior(priors, "sd", group = all_random_effects$group[i], dpar = all_random_effects$param[i])
        prior_intercept <- get_prior(priors, "sd", group = all_random_effects$group[i], dpar = all_random_effects$param[i], coef = "Intercept")
        name <- all_random_effects$param[i]
        gr <- all_random_effects$group[i]
        index <- if(is.null(all_random_effects$index)) 1L else all_random_effects$index[i]
        if(!is.null(prior_intercept) && !is.null(prior_slope)) {
          if(prior_intercept == prior_slope) {
            sprintf("%s ~ %s;", if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s]", name, gr, seq_len(parameters[[name]]$hierarchical$fdim)), prior_intercept)
          } else {
            c(
              sprintf("if(int_%1$s_%2$s) {", name, gr),
              sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s[:(int_%1$s_%2$s-1)]]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s[:(int_%1$s_%2$s-1)]]", name, gr, index)),
              sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s[(int_%1$s_%2$s+1):]]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s[(int_%1$s_%2$s+1):]]", name, gr, index)),
              sprintf("\t%2$s ~ %1$s;", prior_intercept, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s[int_%1$s_%2$s]]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s[int_%1$s_%2$s]]", name, gr, index)),
              "} else {",
              sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s]", name, gr, index)),
              "}"
            )
          }
        } else if(!is.null(prior_intercept)) {
          c(
            sprintf("// no prior for %s random slopes (group: %s)", name, gr),
            sprintf("if(int_%1$s_%2$s) {", name, gr),
            sprintf("\t%2$s ~ %1$s;", prior_intercept, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s[int_%1$s_%2$s]]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s[int_%1$s_%2$s]]", name, gr, index)),
            "}"
          )
        } else if(!is.null(prior_slope)) {
          c(
            sprintf("// no prior for %s random intercepts (group: %s)", name, gr),
            sprintf("if(int_%1$s_%2$s) {", name, gr),
            sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s[:(int_%1$s_%2$s-1)]]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s[:(int_%1$s_%2$s-1)]]", name, gr, index)),
            sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s[(int_%1$s_%2$s+1):]]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s[(int_%1$s_%2$s+1):]]", name, gr, index)),
            "} else {",
            sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("s_%2$s[map_%1$s_%2$s]", name, gr) else sprintf("s_%2$s[map_%1$s_%3$d_%2$s]", name, gr, index)),
            "}"
          )
        } else {
          c(
            sprintf("// no prior for %s random intercepts (group: %s)", name, gr),
            sprintf("// no prior for %s random slopes (group: %s)", name, gr)
          )
        }
      })),
      unlist(lapply(seq_len(nrow(all_random_params)), function(i) {
        p_r <- get_prior(priors, "cor", group = all_random_params$group[i])
        c(
          sprintf("if(M_%s > 1) {", all_random_params$group[i]),
          if(is.null(p_r)) sprintf("\t// no prior for %s random effects correlations", all_random_params$group[i]) else  sprintf("\tr_%1$s ~ %2$s;", all_random_params$group[i], p_r),
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
    paste0("\t", c("real log_lik;", l_body, if(isTRUE(debug_neginf_loglik)) sprintf("if(log_lik == negative_infinity()) print(\"logLik(%s|%s) = -inf !\");", paste0(sprintf("%1$s=\",%1$s,\"",datsig(names_back = l_data, types = FALSE)),collapse=","), paste0(sprintf("%1$s=\",%1$s,\"",parsig(l_pars, types = FALSE)),collapse=",")), "return log_lik;")),
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
    unlist(lapply(names(parameters), function(name) {
      if(!is.null(parameters[[name]]$hierarchical)) {
        prior_intercept <- get_prior(priors, "b", dpar=name, coef = "Intercept")
        prior_slope <- get_prior(priors, "b", dpar=name)
        if(!is.null(prior_intercept) && !is.null(prior_slope)) {
          if(prior_intercept == prior_slope) {
            sprintf("%s ~ %s;", if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s]", name) else sprintf("b[map_%1$s_%2$d]", name, seq_len(parameters[[name]]$hierarchical$fdim)), prior_intercept)
          } else {
            c(
              sprintf("if(int_%1$s) {", name),
              sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s[:(int_%1$s-1)]]", name) else sprintf("b[map_%1$s_%2$d[:(int_%1$s-1)]]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
              sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s[(int_%1$s+1):]]", name) else sprintf("b[map_%1$s_%2$d[(int_%1$s+1):]]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
              sprintf("\t%2$s ~ %1$s;", prior_intercept, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s[int_%1$s]]", name) else sprintf("b[map_%1$s_%2$d[int_%1$s]]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
              "} else {",
              sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s]", name) else sprintf("b[map_%1$s_%2$d]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
              "}"
            )
          }
        } else if(!is.null(prior_intercept)) {
          c(
            sprintf("// no prior for %s fixed slopes", name),
            sprintf("if(int_%1$s) {", name),
            sprintf("\t%2$s ~ %1$s;", prior_intercept, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s[int_%1$s]]", name) else sprintf("b[map_%1$s_%2$d[int_%1$s]]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
            "}"
          )
        } else if(!is.null(prior_slope)) {
          c(
            sprintf("// no prior for %s fixed intercepts", name),
            sprintf("if(int_%1$s) {", name),
            sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s[:(int_%1$s-1)]]", name) else sprintf("b[map_%1$s_%2$d[:(int_%1$s-1)]]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
            sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s[(int_%1$s+1):]]", name) else sprintf("b[map_%1$s_%2$d[(int_%1$s+1):]]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
            "} else {",
            sprintf("\t%2$s ~ %1$s;", prior_slope, if(parameters[[name]]$hierarchical$dim == 1L) sprintf("b[map_%1$s]", name) else sprintf("b[map_%1$s_%2$d]", name, seq_len(parameters[[name]]$hierarchical$fdim))),
            "}"
          )
        } else {
          c(
            sprintf("// no prior for %s fixed intercepts", name),
            sprintf("// no prior for %s fixed slopes", name)
          )
        }
      } else {
        p <- get_prior(priors, "global", dpar=name)
        if(is.null(p)) sprintf("// no prior for global %s", name) else sprintf("%s ~ %s;", name, p)
      }
    }))
  )

  if(isTRUE(parallel)) {
    add_code(
      "transformed data",
      datmap()
    )
    add_code(
      "transformed parameters 2",
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
      "{",
      paste0("\t", code_blocks$`transformed parameters 2`),
      if(isTRUE(parallel)) {
        paste0("\tlog_lik = map_rect(log_lik_rect, phi, theta, x_r, x_i);")
      } else {
        paste0("\tfor(i in 1:N) log_lik[i] = log_lik_single(",paste(c(datsig(names_back = l_data, types = FALSE, index = "i"), parsig(l_pars, types = FALSE, index = "i")),collapse=", "),");")
      },
      "}"
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

  add_code("model", code_blocks$`transformed parameters 2`, prepend = TRUE)

  header <- c(
    "StanTVA",
    "=======",
    sprintf("This is a StanTVA program, generated with %s (v%s) Please cite as:", packageName(), packageVersion(packageName())),
    "",
    strsplit(format(citation(packageName()), style="text"),"\n")[[1]],
    "",
    "Configuration",
    "=============",
    sprintf(" - %s = %s", names(call_args), call_args),
    "",
    "License",
    "=======",
    "This program is licensed under the GNU General Public License 3. For a copy of the license agreement, see: https://www.gnu.org/licenses/gpl-3.0.html"
  ) %>% ansi_strwrap(width = 80L)

  code_blocks$`transformed parameters 2` <- NULL

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
  dfdim <- vapply(names(parameters), function(pn) if(is.null(parameters[[pn]]$dim)) 1L else as.integer(parameters[[pn]]$dim), integer(1), USE.NAMES = TRUE)
  if(has_formulas) {
    attr(df, "formula_lhs") <- formula_lhs
    attr(df, "random_factors") <- all_random_factors
  }
  new("stantvacode", code = ret, config = call_args_list, include_path = stantva_path(), df = df, dim = dfdim, version = packageVersion(packageName()), priors = priors)
}


#'StanTVA code class
#'@slot code The generated Stan code.
#'@slot config A list of model configuration parameters, as passed to \code{stantva_code()} or \code{stantva_model()}.
#'@slot include_path The path to the StanTVA includes (usually identical to \code{stantva_path()}).
#'@slot df The degrees of freedom of the model parameters.
#'@slot dim The dimensions of the model parameters.
#'@slot version The RStanTVA package version that was used to generate this model fit.
#'@slot priors Priors for the model parameters.
#'@export
setClass("stantvacode", slots = c("code" = "character", "config" = "list", "include_path" = "character", "df" = "integer", "dim" = "integer", "version" = "ANY", "priors" = "ANY"))


#'Show StanTVA code
#'
#'Display the content of the StanTVA code object in the console.
#'@param object The StanTVA code object.
#'@returns Returns \code{object} invisibly but the function is usually only called for its side effects.
#'@export
setMethod("show", "stantvacode", function(object) {
  cat(col_grey("// Include path(s): ", paste0(object@include_path, collapse="; ")),"\n")
  cat(object@code)
  invisible(object)
})


#'Read StanTVA fit from CSV
#'@description This function may be used to read an RStan or CmdStan fit from CSV files. Note that you also need to provide the fitted model.
#'@param csv_file The CSV file to be read.
#'@param data The data to which the model was fitted.
#'@param model The fitted model as an StanTVA model or StanTVA code object.
#'@param contrasts Any contrasts specified to factors in the data set.
#'@return The StanTVA fit object.
#'@examples
#'\donttest{
#'data <- read_tva_data("data.dat")
#'model <- stantva_code(locations = 6)
#'fit <- stancsv2stantvafit("chain1.csv", data, model)
#'fit
#'}
#'@export
stancsv2stantvafit <- function(csv_file, data, model, contrasts = list()) {
  mm <- if(inherits(model, "stantvamodel")) model else if(inherits(model, "stantvacode")) {
    mm <- stan_model(model_code = model@code, isystem = stantva_path())
    mm <- as(mm, "stantvamodel")
    mm@code <- model
    mm
  } else stop("model must be stantvamodel or stantvacode!")
  fx <- read_stan_csv(csv_file)
  fx@stanmodel <- mm
  fx <- as(fx, "stantvafit")
  fx@data <- prepare_data(data, mm, contrasts=contrasts)
  fx
}

#' @title StanTVA model
#' @description Creates a StanTVA model object.
#' @param ... Additional arguments passed to \code{\link[RStanTVA:stantva_code]{stantva_code()}}.
#' @param stan_options The Stan options, passed to \code{\link[rstan:stan_model]{stan_model()}}
#' @return The StanTVA model object.
#' @examples
#' \donttest{
#' model <- stantva_model(locations = 2, task = "pr")
#' model
#' }
#' @export
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

#' @title Write StanTVA model
#' @description Writes a StanTVA model to a file.
#' @param model The StanTVA model object.
#' @param file The file name.
#' @return No return value, called for side effects.
#' @examples
#' \donttest{write_stantva_model(model, "model.stan")}
#'@export
write_stantva_model <- function(model, file = stdout()) {
  code <- if(inherits(model, "stantvamodel") || inherits(model, "stantvafit")) {
    model@code
  } else if(inherits(model, "stantvacode")) {
    model
  } else {
    stop("`model` must be of type stantvamodel, stantvafit or stantvacode!")
  }
  writeLines(code@code, file)
}


#'StanTVA model class
#'@slot code The StanTVA code object that was used to compile this model.
#'@importClassesFrom rstan stanmodel
#'@export
setClass("stantvamodel", contains = "stanmodel", slots = c("code" = "stantvacode"))

#'StanTVA fit class
#'@slot stanmodel The StanTVA model object that was fitted to the data.
#'@slot data The data to which the StanTVA model was fitted.
#'@importClassesFrom rstan stanfit
#'@export
setClass("stantvafit", contains = "stanfit", slots = c("stanmodel" = "stantvamodel", "data" = "list"))


#' @title Show StanTVA model
#' @description Prints a StanTVA model object.
#' @param object The StanTVA model object.
#' @return The printed object.
#' @examples
#' \donttest{
#' model <- stantva_model(locations = 4)
#' show(model)
#' }
#'@export
setMethod("show", c(object="stantvamodel"), function(object) {
  cat(col_cyan("StanTVA"), "model with", length(object@code@df),"free parameter(s) and the following configuration:\n")
  for(cname in names(object@code@config)) {
    cat(ansi_strwrap(paste0("- ",col_magenta(cname)," = ",deparse1(if(cname == "priors") deparse_prior(object@code@config$priors) else object@code@config[[cname]])), indent = 2, exdent = 6),sep="\n")
  }
  invisible(object)
})

init_sampler <- function(model, pdata) {
  mc <- if(inherits(model, "stantvamodel")) model@code else if(inherits(model, "stantvacode")) model else stop("`model` must be stantvamodel or stantvacode!")
  function(chain_id = 1) {
    if(is.null(pdata$X)) {
      list(.p = 0L)
    } else {
      ret <- list(.p = 0L)
      #ret$b <- double(ncol(pdata$X))
      #ret$b[colnames(pdata$X) == "Intercept"] <- as.array(runif(sum(colnames(pdata$X) == "Intercept"), -0.2, 0.2))
      rfs <- bind_rows(bind_rows(lapply(mc@config$formula, parse_formula))$random)
      for(i in seq_len(nrow(rfs))) {
        for(j in seq_len(mc@df[rfs$param[i]])) {
          rf <- if(mc@dim[rfs$param[i]] > 1L && !rfs$custom[i]) paste0(rfs$group[i], "_", j) else rfs$group[i]
          #ret[[paste0("r_",rf)]] <- diag(ncol(pdata[[paste0("Z_",rf)]]))
          #ret[[paste0("s_",rf)]] <- as.array(if_else(colnames(pdata[[paste0("Z_",rf)]]) == "Intercept", 0.1, 0.01))
          #message(paste0("N_", rfs$factor_txt[i]),",",paste0("Z_",rf))
          ret[[paste0("w_",rf)]] <- matrix(0, nrow = pdata[[paste0("N_", rfs$factor_txt[i])]], ncol = ncol(pdata[[paste0("Z_",rf)]]))
        }
      }
      ret
    }
  }
}

fix_cmdstanr_output <- function(fp) {
  n <- 0L
  lines <- readLines(fp)
  save_warmup_lines <- grep("^\\s?#.*save_warmup", lines)
  for(i in save_warmup_lines) {
    lines[i] <- gsub("false|False|FALSE","0",lines[i])
  }
  writeLines(lines, fp)
}

#'Draw posterior samples from an RStanTVA model
#'
#'Draw samples from the model defined by \code{object}.
#'@export
setGeneric("sampling")

#'Maximum-likelihood estimation
#'
#'Obtain a point estimate by maximizing the joint posterior from the StanTVA model.
#'@export
setGeneric("optimizing")

#'@describeIn sampling method
#'@param object The StanTVA model object.
#'@param data The data to which the model should be fitted, usually a \code{data.frame}.
#'@param init How to initialize the individual chains, see \code{\link[rstan:sampling]{rstan::sampling()}}. Note that for \code{random}, any lower-level hierarchical (e.g., subject-level) parameters are initialized to zero.
#'@param backend Which backend to use for fitting (default: \code{rstan})
#'@param cpp_options Which options to pass to \code{stan_model()} for compiling the C++ code.
#'@param ... Further arguments passed to the sampling handler of the specified backend.
#'@return Returns a \code{stantva_fit} object, which inherits from \code{\link[rstan:stanfit]{stanfit}}, representing the fit of \code{object} to \code{data}.
#'@export
setMethod("sampling", c(object = "stantvamodel"), function(object, data,init = "random", ..., backend = c("rstan","cmdstanr","cmdstanr_mpi"), cpp_options = if(match.arg(backend) == "cmdstanr") list(stan_threads = object@code@config$parallel) else if(match.arg(backend) == "cmdstanr_mpi") list(CXX = "mpicxx", TBB_CXX_TYPE = "gcc", STAN_MPI = TRUE)) {
  if(object@code@config$locations != ncol(data$S)) stop("Cannot fit a StanTVA model compiled for ",object@code@config$locations," location(s) to a data set with ",ncol(data$S)," location(s)!")
  pdata <- prepare_data(data, object)
  formula_lhs <- attr(object@code@df, "formula_lhs")
  pars_to_exclude <- character()
  backend <- match.arg(backend)
  pars <- NULL
  include <- TRUE
  if(length(pars) == 0 || (length(pars) == 1 && is.na(pars) && isTRUE(include))) {
    include <- FALSE
    pars <- pars_to_exclude
  } else if(isFALSE(include)) {
    pars <- union(pars_to_exclude, na.omit(pars))
  } else {
    pars <- setdiff(na.omit(pars), pars_to_exclude)
  }
  if((missing(init) || identical(init, "random")) && !is.null(pdata$X)) {
    init <- init_sampler(object, pdata)
  } else if(missing(init)) {
    init <- "random"
  }
  if(backend == "rstan") {
    if(length(pars) == 0 && isFALSE(include)) {
      f <- callNextMethod(object, pdata, pars = NA, include = TRUE, init = init, ...)
    } else {
      f <- callNextMethod(object, pdata, pars = pars, include = include, init = init, ...)
    }
    f@stanmodel <- object
    f <- as(f, "stantvafit")
    f@data <- pdata
  } else if(backend == "cmdstanr") {
    m <- cmdstanr::cmdstan_model(stan_file = cmdstanr::write_stan_file(object@code@code), include_paths = stantva_path(), cpp_options = cpp_options)
    x <- m$sample(pdata, init = init, threads_per_chain = rstan_options("threads_per_chain"), save_warmup = 0L, ...)
    for(fp in x$output_files()) {
      fix_cmdstanr_output(fp)
    }
    f <- stancsv2stantvafit(x$output_files(), data, object@code)
  } else if(backend == "cmdstanr_mpi") {
    m <- cmdstanr::cmdstan_model(stan_file = cmdstanr::write_stan_file(object@code@code), include_paths = stantva_path(), cpp_options = cpp_options)
    x <- m$sample_mpi(pdata, init = init, save_warmup = 0L, ...)
    for(fp in x$output_files()) {
      fix_cmdstanr_output(fp)
    }
    f <- stancsv2stantvafit(x$output_files(), data, object@code)
  }
  .striplhspars("f")
  f
})

#'@describeIn optimizing method
#'@param object The StanTVA model object.
#'@param data The data to which the model should be fitted, usually a \code{data.frame}.
#'@param init How to initialize the individual chains, see \code{\link[rstan:optimizing]{rstan::optimizing()}}. Note that for \code{random}, any lower-level hierarchical (e.g., subject-level) parameters are initialized to zero.
#'@param ... Further arguments passed to \code{\link[rstan:optimizing]{rstan::optimizing()}}.
#'@return A list, representing the maximum-likelihood estimate, see \code{\link[rstan:optimizing]{rstan::optimizing()}}.
#'@export
setMethod("optimizing", c(object = "stantvamodel"), function(object, data, init, ...) {
  if(object@code@config$locations != ncol(data$S)) stop("Cannot fit a StanTVA model compiled for ",object@code@config$locations," location(s) to a data set with ",ncol(data$S)," location(s)!")
  formula_lhs <- attr(object@code@df, "formula_lhs")
  pdata <- prepare_data(data, object)
  pars <- NULL
  include <- TRUE
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


#' @title Log-likelihood
#' @description Returns the pointwise log-likelihood of a StanTVA fit.
#' @param object The StanTVA fit.
#' @return The pointwise log likelihood.
#' @examples
#' \donttest{
#' loglik <- logLik(model, data, params)
#' loglik
#' }
#'@export
setMethod("logLik", "stantvafit", function(object) {
  if(!isTRUE(object@stanmodel@code@config$save_log_lik)) stop("StanTVA model must be compiled with `save_log_lik` = TRUE in order to use logLik()!")
  extract(object, "log_lik")$log_lik
})


.strippars <- function(f, p) {
  eval(
    substitute(
      {
        f@sim$pars_oi <- f@model_pars <- setdiff(f@model_pars, p)
        we <- grepl(paste0("^(",paste0(p,collapse="|"),")\\["),f@sim$fnames_oi)
        f@sim$fnames_oi <- f@sim$fnames_oi[!we]
        f@par_dims <- f@par_dims[f@model_pars]
        for(i in seq_along(f@sim$samples)) {
          we <- grepl(paste0("^(",paste0(p,collapse="|"),")\\["),names(f@sim$samples[[i]]))
          f@sim$samples[[i]] <- f@sim$samples[[i]][names(f@sim$samples[[i]])[!we]]
        }
      },
      list(f = as.name(f), p = p)
    ),
    parent.frame(4)
  )
}

.striplhspars <- function(f) {
  eval(
    substitute(
      {
        p <- if(!is.null(attr(fx@stanmodel@code@df,"formula_lhs"))) attr(fx@stanmodel@code@df,"formula_lhs")[2,] else c()
        .strippars(fn, p)
      },
      list(fx = as.name(f), fn = as.character(f))
    ),
    parent.frame(1)
  )

}

#' @title Read StanTVA fit
#' @description Reads a StanTVA fit object from one or more files. If multiple files are given, the fits will be combined into a single fit object (e.g., combining separately fitted chains).
#' @param files The file names.
#' @return The StanTVA fit object.
#' @examples
#' \donttest{
#' fit <- read_stantva_fit(c("chain1.rds", "chain2.rds"))
#' fit
#' }
#' @export
read_stantva_fit <- function(files) {
  if(length(files) < 1) stop("`files` must contain at least one file name!")
  if(!is.character(files)) stop("`files` must be a character vector of `stantvafit` files!")
  if(any(!file.exists(files))) stop("At least one of the files does not exist!")
  message("Read ",files[[1]],"...")
  a <- readRDS(files[[1]])
  .striplhspars("a")
  if(!inherits(a, "stantvafit")) stop("Files must be `stantvafit` objects!")
  for(bf in files[-1]) {
    message("Read and append ",bf,"...")
    b <- readRDS(bf)
    .striplhspars("b")
    if(!inherits(b, "stantvafit")) stop("`fits` must be a list of `stantvafit` objects!")
    if(b@sim$iter != a@sim$iter) stop("All `stantvafit` objects must have the same number of iterations!")
    if(b@sim$warmup != a@sim$warmup) stop("All `stantvafit` objects must have the same number of warmup iterations!")
    a@sim$chains <- a@sim$chains + b@sim$chains
    a@sim$samples <- c(a@sim$samples, b@sim$samples)
    a@sim$n_save <- c(a@sim$n_save, b@sim$n_save)
    a@sim$warmup2 <- c(a@sim$warmup2, b@sim$warmup2)
    a@sim$permutation <- c(a@sim$permutation, b@sim$permutation)
    a@inits <- c(a@inits, b@inits)
    a@stan_args <- c(a@stan_args, b@stan_args)
    a@date <- max(a@date, b@date)
    b <- NULL
    gc(verbose = FALSE)
  }
  return(a)
}

#' @title Write StanTVA fit
#' @description Writes a StanTVA fit object to a file.
#' @param fit The StanTVA fit object.
#' @param file The file name.
#' @param ... Additional arguments passed to \code{\link[base:saveRDS]{saveRDS()}}.
#' @return No return value, called for side effects.
#' @examples
#' \donttest{write_stantva_fit(fit, "fit.rds")}
#' @export
write_stantva_fit <- function(fit, file, ...) if(inherits(fit, "stantvafit")) saveRDS(fit, file, ...) else stop("`fit` must be a `stantvafit` object!")


alias.stantvafit <- function(object) {
  attr(names(object), "alias")
}

#' @title Retrieve parameters aliases
#' @description Returns the StanTVA parameter aliases for the underlying RStan fit.
#' @param object The StanTVA fit object.
#' @return A character vector of parameter aliases.
#' @examples
#' \donttest{
#' al <- alias(fit)
#' al
#' }
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

#' @title Fixed effects
#' @description Returns the fixed effects for a StanTVA fit object.
#' @param object The StanTVA fit object.
#' @return The fixed effects.
#' @examples
#' \donttest{
#' fixed_effects <- fixef(fit)
#' fixed_effects
#' }
#'@export
setMethod("fixef", "stantvafit", fixef.stantvafit)


ranef.stantvafit <- function(object) {
  all_params <- bind_rows(lapply(object@stanmodel@code@config$formula, parse_formula)) %>% rename(name = .data$param) %>% mutate(dim = object@stanmodel@code@dim[.data$name], fdim = if_else(.data$name %in% names(object@stanmodel@code@df), object@stanmodel@code@df[.data$name], object@stanmodel@code@dim[.data$name]))
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
    gs <- unique(filter(g, .data$factor_txt == rf)$group)
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
          dimnames(bo)[[3]][m[,i]] <- sprintf("%s_%s[%d]", g$param[i], colnames(object@data[[sprintf("Z_%s", o)]])[m[,i]], g$index[i])
        }
      }
      b[,,y:(y+dim(bo)[[3]]-1L)] <- bo
      dimnames(b)[[3]][y:(y+dim(bo)[[3]]-1L)] <- dimnames(bo)[[3]]
      y <- y + dim(bo)[[3]]
    }
    b
  }, simplify = FALSE)
}

#' @title Random effects
#' @description Returns the random effects for a StanTVA fit object.
#' @param object The StanTVA fit object.
#' @return The fixed effects.
#' @examples
#' \donttest{
#' random_effects <- ranef(fit)
#' random_effects
#' }
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


#' @title Model coefficients
#' @description Returns the model coefficients (sum of fixed + random effects, grouped by random factor) for a StanTVA fit object.
#' @param object The StanTVA fit object.
#' @return The model coefficients, grouped by random factors.
#' @examples
#' \donttest{
#' fixef <- coef(fit)
#' fixef
#' }
#'@export
setMethod("coef", "stantvafit", coef.stanfit)

predict.stantvafit <- function(object, newdata, variables = names(object@stanmodel@code@df)) {
  newdata <- if(missing(newdata)) object@data else prepare_data(newdata, object@stanmodel, FALSE, attr(object@data, "contrasts"))
  fx <- bind_rows(lapply(object@stanmodel@code@config$formula, parse_formula))
  sapply(variables, function(parname) {
    which_formula <- match(parname, fx$param)
    if(is.na(which_formula)) {
      p <- extract(object, parname)[[1]]
      replicate(newdata$N, p)
    } else {
      #### !!!!!!!!
      bt <- fx$inverse_link[[which_formula]]
      rfs <- fx$random[[which_formula]]
      par_dim <- object@stanmodel@code@dim[parname]
      par_df <- object@stanmodel@code@df[parname]
      ps <- "b"
      for(i in seq_len(nrow(rfs))) {
        if(par_dim > 1) ps <- c(ps, paste0("w_",rfs$group[i],"_",seq_len(par_df)))
        else ps <- c(ps, paste0("w_",rfs$group[i]))
      }
      p <- extract(object, ps)
      r <- vapply(seq_len(par_dim), function(i) {
        m <- newdata[[if(par_dim > 1) paste0("map_",parname,"_",i) else paste0("map_",parname)]]
        if(i > par_df) return(matrix(1, (object@sim$iter-object@sim$warmup)*object@sim$chains, newdata$N))
        y <- tcrossprod(newdata$X[,m], p$b[,m])
        for(k in seq_len(nrow(rfs))) {
          mrf <- newdata[[if(par_dim > 1) paste0("map_",parname,"_",i,"_",rfs$group[k],"_",i) else paste0("map_",parname,"_",rfs$group[k])]]
          for(j in seq_len(newdata[[paste0("N_",rfs$factor_txt[k])]])) {
            J <- newdata[[rfs$factor_txt[k]]] == j
            #print(which(J))
            #print(if(par_dim > 1) paste0("w_",rfs$group[k],"_",i) else paste0("w_",rfs$group[k]))
            Z <- newdata[[if(par_dim > 1) paste0("Z_",rfs$group[k],"_",i) else paste0("Z_",rfs$group[k])]][J,mrf,drop=FALSE]
            #message(if(par_dim > 1) paste0("w_",rfs$group[k],"_",i) else paste0("w_",rfs$group[k]))
            w <- t(as.matrix(p[[if(par_dim > 1) paste0("w_",rfs$group[k],"_",i) else paste0("w_",rfs$group[k])]][,j,mrf]))
            W <- Z %*% w
            #print(if(par_dim > 1) paste0("Z_",rfs$group[k],"_",i) else paste0("Z_",rfs$group[k]))
            #print(names(newdata))
            #print(Z)
            y[J,] <- y[J,] + W
          }
        }
        t(fx$inverse_link[[which_formula]](y))
      }, matrix(NA_real_, (object@sim$iter-object@sim$warmup)*object@sim$chains, newdata$N))
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

#' @title Predict parameter values
#' @description Returns the predictions for latent model parameters.
#' @param object The StanTVA fit object.
#' @param newdata The new data (leave empty to use fitted data).
#' @param variables The names of the parameters to predict.
#' @return The predictions.
#' @examples
#' \donttest{
#' p <- predict(fit, variables = c("C","K"))
#' colMeans(p$C)
#' }
#'@export
setMethod("predict", "stantvafit", predict.stantvafit)


fitted.stantvafit <- function(object, variables = names(object@stanmodel@code@df)) {
  predict(object, variables)
}


#' @title Retrieve fitted parameter values
#' @description Returns the fitted values for latent model parameters. This is identical to calling \code{predict()} without new data.
#' @param object The StanTVA fit object.
#' @param variables The names of the parameters to retrieve.
#' @return The fitted values.
#' @examples
#' \donttest{
#' p <- fitted(fit, variables = c("C","K"))
#' colMeans(p$C)
#' }
#'@export
setMethod("fitted", "stantvafit", fitted.stantvafit)

#' @title Generate typical descriptive statistics for TVA reports
#' @description This function generates by-trial descriptive statistics, see `Value` below.
#' @param data The TVA report data as a \code{data.frame}.
#' @return The function returns a transmuted \code{data.frame}/\code{tibble} with columns \code{condition} (copied from \code{data}), \code{exposure} (copied from \code{data$T}), \code{n_items}, \code{n_targets}, \code{n_distractors}, and \code{score} (number of correctly reported items).
#' @examples
#' \donttest{tva_report(tva_recovery)}
#' @export
tva_report <- function(data) {
  data %>% transmute(
    condition = .data$condition,
    exposure = .data$T,
    score = as.integer(if(is.null(.data$D)) rowSums(.data$R & .data$S) else rowSums(.data$R & .data$S & !.data$D)),
    n_items = as.integer(rowSums(.data$S == 1L)),
    n_distractors = if(is.null(.data$D)) integer(n()) else as.integer(rowSums(.data$D))
  ) %>% mutate(n_targets = .data$n_items - .data$n_distractors)
}

#' @describeIn print-stantvafit-method Alias
#' @param object The StanTVA fit object.
#'@export
setMethod("show", "stantvafit", function(object) print(object))


#' Print StanTVA fit
#' @description Prints a StanTVA fit object.
#' @param x The StanTVA fit object.
#' @param digits_summary The number of significant digits to display in posterior summaries.
#' @param ... Currently not used.
#' @return Returns \code{x}. Usually called for its side effects (printing to the console).
#' @examples
#' \donttest{
#' print(fit)
#' }
#'@export
setMethod("print", "stantvafit", function(x, digits_summary = 2, ...) {

  par_names <- names(x)

  data <- x@data

  x <- as(x, "stanfit")

  sampler <- attr(x@sim$samples[[1]], "args")$sampler_t

  cat(col_cyan("StanTVA"), "model with", length(x@stanmodel@code@df), "free parameter(s), fitted with ")
  cat(x@sim$chains," ",sampler, " chains, each with iter=", x@sim$iter,
      "; warmup=", x@sim$warmup, "; thin=", x@sim$thin,"\n", sep = "")


  heading <- function(txt) cat("\n", style_underline(style_bold(txt)), "\n", sep = "")

  heading("Model configuration:")
  for(n in names(x@stanmodel@code@config)) {
    cat(n,"= ")
    if(n == "priors") cat(deparse1(deparse_prior(x@stanmodel@code@config$prior)))
    else cat(deparse1(x@stanmodel@code@config[[n]]))
    cat("\n")
  }

  fx <- bind_rows(lapply(x@stanmodel@code@config$formula, parse_formula))


  global_pars <- if(is.null(fx$param)) names(x@stanmodel@code@df) else setdiff(names(x@stanmodel@code@df), fx$param)

  not_converged <- character()

  if(length(global_pars) > 0) {

    heading("Global parameters:")

    g_summary <- rstan::summary(x, global_pars, use_cache = FALSE)$summary

    print(round(g_summary, digits_summary), max = prod(dim(g_summary)))

    not_converged <- c(not_converged, rownames(g_summary)[g_summary[,"Rhat"] >= 1.05])

  }


  if(nrow(fx) > 0L) {

    rfs <- fx$random %>% lapply(function(fxi) if(!is.null(fxi$param) && x@stanmodel@code@dim[fxi$param] > 1L) crossing(fxi, index = seq_len(x@stanmodel@code@df[fxi$param])) else bind_cols(fxi, index = NA_integer_)) %>% bind_rows(tibble(group = character(), custom = logical(), index = integer())) %>% mutate(group = if_else(.data$custom | is.na(.data$index), .data$group, paste0(.data$group,"_",.data$index)))

    random_factors_txt <- if(nrow(rfs) > 0L) unique(rfs$factor_txt) else character()

    b_summary <- rstan::summary(x, "b", use_cache = FALSE)$summary

    rownames(b_summary) <- attr(par_names, "alias")[match(rownames(b_summary), par_names)]

    heading("Fixed effects:")

    print(round(b_summary, digits_summary), max = prod(dim(b_summary)))

    not_converged <- c(not_converged, rownames(b_summary)[b_summary[,"Rhat"] >= 1.05])


    for(rf in random_factors_txt) {

      heading(paste0("Hyperparameters on random effects (", col_blue(deparse1(rfs$factor[[match(rf, rfs$factor_txt)]])), " level, N = ", data[[paste0("N_",rf)]] ,"):"))



      gs <- unique(rfs$group[rfs$factor_txt == rf])



      s_summary <- rstan::summary(x, paste0("s_", gs), use_cache = FALSE)$summary

      rownames(s_summary) <- attr(par_names, "alias")[match(rownames(s_summary), par_names)]

      for(g in gs) {

        M_rf <- x@par_dims[[sprintf("r_%s", g)]][1]
        if(M_rf >= 2L) {
          c_pars <- combn(M_rf, 2L)
          r_summary <- rstan::summary(x, paste0("r_", g, "[", c_pars[1,], ",", c_pars[2,],"]"), use_cache = FALSE)$summary
          rownames(r_summary) <- attr(par_names, "alias")[match(rownames(r_summary), par_names)]
          s_summary <- rbind(s_summary, r_summary)
        }
      }

      print(round(s_summary, digits_summary), max = prod(dim(s_summary)))



      not_converged <- c(not_converged, rownames(s_summary)[s_summary[,"Rhat"] >= 1.05])


      w_summary <- rstan::summary(x, sprintf("w_%s", gs), probs = double(0), use_cache = FALSE)$summary

      not_converged <- c(not_converged, attr(par_names, "alias")[match(rownames(w_summary)[w_summary[,"Rhat"] >= 1.05], par_names)])


    }


  }



  if(length(not_converged)) {
    warning("Model did not converge (Rhat >= 1.05) for ",length(not_converged)," parameter(s): ", paste(not_converged, collapse=", "))
  }

  invisible(x)

})

translate_names <- function(model, data, names) {

  ret <- names
  new_ret <- names

  fx <- bind_rows(lapply(model@code@config$formula, parse_formula))

  if(nrow(fx) > 0L) {
    for(i in seq_len(nrow(fx))) {
      param <- fx$param[i]
      if(model@code@dim[param] == 1L) {
        m <- data[[sprintf("map_%s", param)]]
        new_ret[match(sprintf("b[%d]", m), ret)] <- paste0(fx$param[i], "_", colnames(data$X)[m])
      } else for(j in seq_len(model@code@df[param])) {
        m <- data[[sprintf("map_%s_%d", param, j)]]
        new_ret[match(sprintf("b[%d]", m), ret)] <- paste0(fx$param[i], "_", colnames(data$X)[m], "[",j,"]")
      }
    }
    rfs <- if(is.null(fx$random) || nrow(bind_rows(fx$random)) == 0) tibble() else fx$random %>% lapply(function(fxi) if(!is.null(fxi$param) && model@code@dim[fxi$param] > 1L) crossing(fxi, index = seq_len(model@code@dim[fxi$param])) else bind_cols(fxi, index = NA_integer_)) %>% bind_rows(tibble(group = character(), custom = logical(), index = integer())) %>% mutate(group = if_else(.data$custom | is.na(.data$index), .data$group, paste0(.data$group,"_",.data$index)))
    for(i in seq_len(nrow(rfs))) {
      param <- rfs$param[i]
      if(model@code@dim[param] == 1L) {
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
        if(model@code@dim[param] == 1L) {
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

#'Retrieve model parameter names
#'@description Returns the names of the fitted model parameters.
#'@param x The StanTVA fit.
#'@return The list of parameter names and aliases.
#'@examples
#'\donttest{
#'f <- read_stantva_fit("fit.rds")
#'names(f)
#'}
#'@export
setMethod("names", "stantvafit", function(x) translate_names(x@stanmodel, x@data, x@sim$fnames_oi))


#'Extract samples from a fitted RStanTVA model
#'@description Returns posterior samples from a fitted RStanTVA model.
#'@param object The RStanTVA fit.
#'@param pars (Optional) A character vector of variable names to extract.
#'@param ... Additional arguments passed to \code{\link[rstan:extract]{rstan::extract()}}, e.g. \code{permuted} and \code{inc_warmup}.
#'@return See \code{\link[rstan:extract]{rstan::extract()}} for details.
#'@examples
#'\donttest{
#'f <- read_stantva_fit("fit.rds")
#'extract(f, "C_Intercept")
#'}
#'@export
setMethod("extract", c(object="stantvafit"), function(object, pars, ...) {
  f_names <- names(object)
  x <- if(missing(pars)) callNextMethod(object, ...) else callNextMethod(object, object@sim$fnames_oi[match(pars, alias(object))], ...)
  if(is.list(x)) names(x) <- attr(f_names,"alias")[match(names(x), f_names)]
  if(is.array(x)) rownames(x) <- attr(f_names,"alias")[match(dimnames(x)[1], f_names)]
  x
})

#'Summary method for RStanTVA fits
#'@description Summarize the distributions of estimated parameters and derived quantities using the posterior draws.
#'@param object The RStanTVA fit.
#'@param pars (Optional) A character vector of variable names to extract.
#'@param ... Additional arguments passed to \code{\link[rstan:summary,stanfit-method]{rstan::summary()}}, e.g. \code{probs} and \code{use_cache}.
#'@return See \code{\link[rstan:summary,stanfit-method]{rstan::summary()}} for details.
#'@examples
#'\donttest{
#'f <- read_stantva_fit("fit.rds")
#'summary(f, "C_Intercept", probs = c(.025, .975))
#'}
#'@export
setMethod("summary", c(object="stantvafit"), function(object, pars, ...) {
  f_names <- names(object)
  x <- if(missing(pars)) rstan::summary(as(object,"stanfit"), ...) else rstan::summary(as(object,"stanfit"), object@sim$fnames_oi[match(pars, alias(object))], ...)
  rownames(x$summary) <- attr(f_names,"alias")[match(rownames(x$summary), f_names)]
  rownames(x$c_summary) <- attr(f_names,"alias")[match(rownames(x$c_summary), f_names)]
  x
})
