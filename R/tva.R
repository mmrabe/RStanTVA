#'@importFrom rstan extract stan_model sampling optimizing gqs
#'@importFrom dplyr summarize mutate group_by %>% across if_else select bind_cols bind_rows
#'@importFrom tidyr pivot_longer pivot_wider
#'@importFrom readr read_table write_tsv
#'@importFrom methods formalArgs new as callNextMethod
#'@importFrom stats na.omit
#'@importFrom cli col_cyan col_magenta col_grey ansi_strwrap
#'@importFrom tibble tibble
#'@importFrom utils citation str
#'



#'@export
stantva_path <- function() {
  file.path(find.package("RStanTVA"), "StanTVA")
}

#'@export
read_tva_data <- function(file, ...) {
  if(inherits(file, "connection")) f <- file
  else f <- base::file(file, "rb")
  n <- as.integer(readLines(f, n = 1))
  dat <- read_table(f, col_names = c("condition", "exposure", "targets", "distractors", "report"), col_types = c("inccc"), n_max = n, ...)
  if(n != nrow(dat)) warning("Expected ",n," trials but only read ",nrow(dat)," lines!")
  close(f)
  dat %>% mutate(
    across(c(targets, distractors), ~do.call(rbind, strsplit(.x, "", TRUE) %>% lapply(function(x) if_else(x=="0",NA_character_,x)))),
    report = strsplit(if_else(report == "-", "", report), "", TRUE),
    S = (!is.na(targets) | !is.na(distractors))+0L,
    D = (!is.na(distractors) )+ 0L,
    R = t(vapply(seq_len(n()), function(i) targets[i,] %in% report[[i]] | distractors[i,] %in% report[[i]], logical(ncol(targets)))) + 0L,
    items = t(vapply(seq_len(n()), function(i) if_else(is.na(targets[i,]), distractors[i,], targets[i,]), character(ncol(targets))))
  ) %>% select(condition, S, D, items, T = exposure, R) %>% {
    N <- nrow(.)
    . <- as.list(.)
    .$N <- N
    .$locations <- ncol(.$S)
    .
  } %>% as("tvadata")
}


#'@export
write_tva_data <- function(data, file, ...) {
  if(inherits(file, "connection")) f <- file
  else f <- base::file(file, "w")
  stopifnot(data$N == nrow(data$S))
  stopifnot(data$N == nrow(data$R))
  stopifnot(is.null(data$D) || data$N == nrow(data$D))
  stopifnot(ncol(data$S) == ncol(data$R))
  stopifnot(is.null(data$D) || ncol(data$S) == ncol(data$D))
  if(is.null(data$items)) {
    data$items <- t(vapply(seq_len(data$N), function(i) sample(LETTERS, ncol(data$S)), character(ncol(data$S))))
  }
  if(is.null(data$condition)) {
    data$condition <- rep(1L, data$N)
  }
  if(is.null(data$D)) {
    data$D <- matrix(0L, nrow = data$N, ncol = ncol(data$S))
  }
  writeLines(as.character(data$N), f)
  bind_cols(
    condition = data$condition,
    exposure = data$T,
    targets = vapply(seq_len(data$N), function(i) paste0(if_else(data$S[i,] == 1L & data$D[i,] == 0L, data$items[i,], "0"), collapse = ""), character(1)),
    distractors = vapply(seq_len(data$N), function(i) paste0(if_else(data$S[i,] == 1L & data$D[i,] == 1L, data$items[i,], "0"), collapse = ""), character(1)),
    report = vapply(seq_len(data$N), function(i) if(sum(data$R[i,])>0) paste0(if_else(data$S[i,] == 1L & data$D[i,] == 0L & data$R[i,] == 1L, data$items[i,], ""), collapse = "") else "-", character(1))
  ) %>% write_tsv(f, col_names = FALSE, ...)
}

setGeneric("model_code", function(object, ...) {})

#'@export
setMethod("model_code", "stanmodel", function(object, type = c("stan","stan2","cpp")) {
  type <- match.arg(type)
  if(type == "stan") object@code
  else if(type == "stan2") object@model_code
  else if(type == "cpp") object@model_cpp
  else stop("Unknown type “", type,"”!")
})

#'@export
setMethod("model_code", "stanfit", function(object, type) {
  model_code(object@stanmodel, type)
})

#'@export
stantva_code <- function(locations, task = c("wr","pr"), regions = list(), C_mode = c("equal","locations","regions"), w_mode = c("equal","locations","regions"), t0_mode = c("constant", "gaussian"), K_mode = c("bernoulli", "free", "binomial", "hypergeometric"), parallel = FALSE, save_log_lik = FALSE, predict_scores = FALSE, priors = FALSE, sanity_checks = FALSE, simulate = FALSE, type = c("stan","stan2","cpp")) {

  task <- match.arg(task)
  C_mode <- match.arg(C_mode)
  t0_mode <- match.arg(t0_mode)
  K_mode <- match.arg(K_mode)
  w_mode <- match.arg(w_mode)
  type <- match.arg(type)

  if(isTRUE(parallel) && rstan_options("threads_per_chain") <= 1L) {
    warning("You requested a parallel model but you also need to configure RStan to run the model in parallel. Try `rstan_options(threads_per_chain = ...)` to set the appropriate number of parallel threads within each chain before compiling the model code! To use all available CPUs, try `rstan_options(threads_per_chain = parallel::detectCores())`.")
  }




  includeFile <- function(f) sprintf("#include %s", f)

  call_args <- character(0)
  call_args_list <- list()

  for(x in formalArgs(sys.function())) {
    val <- get(x)
    call_args[x] <- if(is.numeric(val) || is.character(val) || is.logical(val)) as.character(val) else deparse(val)
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
  data <- list()

  add_code <- function(name, ..., prepend = FALSE) code_blocks[[name]] <<- if(isTRUE(prepend))  c(..., code_blocks[[name]]) else c(code_blocks[[name]], ...)
  add_param <- function(name, ...) {
    parameters[[name]] <<- list(...)
    add_code(if(isTRUE(parameters[[name]]$transformed)) "transformed parameters" else "parameters", sprintf("%s %s;", parameters[[name]]$type, name))
  }
  add_data <- function(name, ...) {
    data[[name]] <<- list(...)
    add_code(if(isTRUE(data[[name]]$transformed)) "transformed data" else "data", sprintf("%s %s;", data[[name]]$type, name))
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

  datmap <- function(names = base::names(data)) {
    names_x_i <- Filter(function(name) !is.null(data[[name]]$class) && data[[name]]$class == "x_i", names)
    names_x_r <- Filter(function(name) !is.null(data[[name]]$class) && data[[name]]$class == "x_r", names)
    par_size_x_i <- vapply(data[names_x_i], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_i <- c(1L,1L+cumsum(par_size_x_i[-length(par_size_x_i)]))
    par_size_x_r <- vapply(data[names_x_r], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_r <- c(1L,1L+cumsum(par_size_x_r[-length(par_size_x_r)]))
    c(
      sprintf("array[N,%d] int x_i;", sum(par_size_x_i)),
      vapply(seq_along(names_x_i), function(i) if(par_size_x_i[i] == 1) sprintf("x_i[:,%d] = %s;", par_offset_x_i[i], names_x_i[i]) else sprintf("x_i[:,%d:%d] = %s;", par_offset_x_i[i], par_offset_x_i[i]+par_size_x_i[i]-1L, names_x_i[i]), character(1)),
      sprintf("array[N,%d] real x_r;", sum(par_size_x_r)),
      vapply(seq_along(names_x_r), function(i) if(par_size_x_r[i] == 1) sprintf("x_r[:,%d] = %s;", par_offset_x_r[i], names_x_r[i]) else sprintf("x_r[:,%d:%d] = %s;", par_offset_x_r[i], par_offset_x_r[i]+par_size_x_r[i]-1L, names_x_r[i]), character(1))
    )
  }

  datremap <- function(names = base::names(data), names_back = names) {
    names_x_i <- Filter(function(name) !is.null(data[[name]]$class) && "x_i" %in% data[[name]]$class, names)
    names_x_r <- Filter(function(name) !is.null(data[[name]]$class) && "x_r" %in% data[[name]]$class, names)
    which_names_back_x_i <- na.omit(match(names_back, names_x_i))
    which_names_back_x_r <- na.omit(match(names_back, names_x_r))
    par_size_x_i <- vapply(data[names_x_i], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_i <- c(1L,1L+cumsum(par_size_x_i[-length(par_size_x_i)]))
    par_size_x_r <- vapply(data[names_x_r], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
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
      vapply(seq_along(names_theta), function(i) if(par_size_theta[i] == 1) sprintf("theta[:,%d] = %s;", par_offset_theta[i], names_theta[i]) else sprintf("theta[:,%d:%d] = %s;", par_offset_theta[i], par_offset_theta[i]+par_size_theta[i]-1L, names_theta[i]), character(1))
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
      vapply(which_names_back_phi, function(i) if(isTRUE(types)) sprintf("%s %s", parameters[[names_phi[i]]]$rtype, names_phi[i]) else names_phi[i], character(1)),
      vapply(which_names_back_theta, function(i) if(isTRUE(types)) sprintf("%s %s", parameters[[names_theta[i]]]$rtype, names_theta[i]) else if(!is.null(index)) sprintf("%s[%s]", names_theta[i], index) else names_theta[i], character(1))
    )
  }

  datsig <- function(names = base::names(data), names_back = names, types = TRUE, index = NULL) {
    names_x_i <- Filter(function(name) !is.null(data[[name]]$class) && "x_i" %in% data[[name]]$class, names)
    names_x_r <- Filter(function(name) !is.null(data[[name]]$class) && "x_r" %in% data[[name]]$class, names)
    which_names_back_x_i <- na.omit(match(names_back, names_x_i))
    which_names_back_x_r <- na.omit(match(names_back, names_x_r))
    par_size_x_i <- vapply(data[names_x_i], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_i <- c(1L,1L+cumsum(par_size_x_i[-length(par_size_x_i)]))
    par_size_x_r <- vapply(data[names_x_r], function(item) if(is.null(item$dim)) 1L else as.integer(item$dim), integer(1))
    par_offset_x_r <- c(1L,1L+cumsum(par_size_x_r[-length(par_size_x_r)]))
    c(
      vapply(which_names_back_x_i, function(i) if(isTRUE(types)) sprintf("%s %s", data[[names_x_i[i]]]$rtype, names_x_i[i]) else if(!is.null(index)) sprintf("%s[%s]", names_x_i[i], index) else names_x_i[i], character(1)),
      vapply(which_names_back_x_r, function(i) if(isTRUE(types)) sprintf("%s %s", data[[names_x_r[i]]]$rtype, names_x_r[i]) else if(!is.null(index)) sprintf("%s[%s]", names_x_r[i], index) else names_x_r[i], character(1))
    )
  }


  add_code("functions", includeFile("tva.stan"))

  add_data(name = "N", type = "int<lower=1>")

  if(C_mode == "equal") {
    add_param(name = "C", type = "real<lower=machine_precision()>", ctype = "real", rtype = "real", prior = ~gamma(3.5,0.035))
    add_param(name = "s", type = sprintf("vector[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype = "vector", dim = locations, transformed = TRUE)
    add_code(
      "transformed parameters",
      sprintf("s = rep_vector(C, %1$d);", locations)
    )
  } else if(C_mode == "regions") {
    if(length(regions) == 0) stop("You must define regions if C_mode == `regions`!")
    add_param(name = "s", type = sprintf("vector[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype = "vector", dim = locations, transformed = TRUE)
    for(rname in names(regions)) {
      add_param(name = paste0("C_", rname), type = "real<lower=machine_precision()>", ctype = "real", rtype = "real", prior = ~gamma(3.5,0.035))
      add_code(
        "transformed parameters",
        sprintf("s[%d] = C_%s;", regions[[rname]], rname)
      )
      add_code(
        "generated quantities",
        paste0("real A_", rname," = (", paste0("s[", regions[[rname]], "]", collapse=" + "),")/sum(s);")
      )
    }
  } else if(C_mode == "locations") {
    add_param(name = "s", type = sprintf("vector<lower=machine_precision()>[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype="vector", dim = locations, prior = ~gamma(3.5,0.035))
    for(i in seq_along(regions)) {
      add_code(
        "generated quantities",
        paste0("real A_", names(regions)[i]," = (", paste0("s[", regions[[i]], "]", collapse=" + "),")/sum(s);")
      )
    }
  }

  if(w_mode == "equal") {
    add_param(name = "w", type = sprintf("vector[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype = "vector", dim = locations, transformed = TRUE)
    add_code(
      "transformed parameters",
      sprintf("w = rep_vector(1.0/%1$d.0, %1$d);", locations)
    )
  } else if(w_mode == "regions") {
    if(length(regions) == 0) stop("You must define regions if w_mode = “regions”!")
    add_param(name = "w", type = sprintf("vector[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype = "vector", dim = locations, transformed = TRUE)
    add_param(name = "b", type = sprintf("simplex[%d]", length(regions)), ctype=sprintf("vector[%d]", length(regions)), rtype = "vector", dim = length(regions), prior = substitute(~dirichlet(ones_vector(nR)), list(nR = length(regions))) )
    for(i in seq_along(regions)) {
      add_code(
        "transformed parameters",
        sprintf("w[%d] = b[%d]/%d.0;", regions[[i]], i, length(regions[[i]]))
      )
      add_code(
        "generated quantities",
        paste0("real w_", names(regions)[i]," = b[", i,"];")
      )
    }
  } else if(w_mode == "locations") {
    add_param(name = "w", type = sprintf("simplex[%d]", locations), ctype=sprintf("vector[%d]", locations), rtype ="vector", dim = locations, prior = substitute(~dirichlet(ones_vector(nS)), list(nS=locations)))
    for(i in seq_along(regions)) {
      add_code(
        "generated quantities",
        paste0("real w_", names(regions)[i]," = ", paste0("w[", regions[[i]], "]", collapse=" + "),";")
      )
    }
  }

  K_args <- "[]'"
  if(K_mode == "bernoulli") {
    add_code("functions", includeFile("bernoulliK.stan"))
    add_param(name = "K", class = c("phi", "K"), type = sprintf("real<lower=0,upper=%d>", locations), ctype="real", rtype="real", prior = substitute(~uniform(0,nS), list(nS = locations)))
    K_args <- "[K]'"
  } else if(K_mode == "free") {
    add_code("functions", includeFile("freeK.stan"))
    add_param(name = "pK", class = c("phi","K"), type = sprintf("simplex[%d]", locations+1L), ctype=sprintf("vector[%d]", locations+1L), rtype="vector", dim = locations+1L, prior = substitute(~dirichlet(ones_vector(nS)), list(nS = locations+1L)))
    add_code("generated quantities", paste0("real mK = ",paste(sprintf("%d * pK[%d]", seq_len(locations), seq_len(locations)+1L), collapse=" + "),";"));
    K_args <- "pK"
  } else if(K_mode == "binomial") {
    add_code("functions", includeFile("binomialK.stan"))
    # TODO add prior!
    add_param(name = "nK", class = c("phi", "K"), type = "real<lower=machine_precision()>", ctype="real", rtype="real")
    add_param(name = "pK", class = c("phi", "K"), type = "real<lower=machine_precision(),upper=1.0-machine_precision()>", ctype="real", rtype="real", prior = ~uniform(0,1))
    add_code("generated quantities", "real mK = nK * pK;")
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
    add_code("generated quantities", sprintf("real mK = %d * gK / (gK + bK);", locations))
    K_args <- "[gK, bK]'"
  }

  t0_args <- "[]'"
  if(t0_mode == "constant") {
    add_code("functions", includeFile("constantt0.stan"))
    add_data(name = "max_mu", type = "real", ctype = "real", rtype="real", transformed = TRUE)
    add_code("transformed data", "max_mu = max(T);", "for(i in 1:N) if(sum(R[i,]) && T[i] < max_mu) max_mu = T[i];")
    add_param(name = "t0", class = c("phi"), type = "real<upper=max_mu>", ctype="real", rtype="real", prior = ~normal(20, 30))
  } else if(t0_mode == "gaussian") {
    #add_code("transformed data", "real t0 = 0.0;")
    add_code("functions", includeFile("gaussiant0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real", ctype="real", rtype="real", prior = ~normal(20, 30))
    add_param(name = "sigma0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real", prior = ~gamma(2,0.08))
    t0_args <- "[mu0, sigma0]'"
  } else if(t0_mode == "exponential") {
    #add_code("transformed data", "real t0 = 0.0;")
    add_code("functions", includeFile("exponentialt0.stan"))
    add_param(name = "mu0", class = c("phi","t0"), type = "real<lower=machine_precision()>", ctype="real", rtype="real", prior = ~normal(20, 30))
    t0_args <- "[mu0]'"
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
    v_pars <- c("w","s")
    v_body <- c(
      "array[nS] int Ss = get_matches(S);",
      "vector[nS] v = s[Ss] .* w[Ss] / sum(w[Ss]);",
      "for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();",
      "return v/1000.0;"
    )
    l_data <- union(c("S","R","T","nS"), v_data)
    l_pars <- c(v_pars, if(!is.null(parameters$t0))"t0",Filter(function(p) any(c("t0","K") %in% parameters[[p]]$class), names(parameters)))
    l_body <- c(
      sprintf("vector[nS] v = calculate_v(%s);", paste(c(datsig(names_back = v_data, types = FALSE), parsig(v_pars, types = FALSE)), collapse=", ")),
      sprintf("return tva_wr_log(R, S, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
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
    v_pars <- c("w","s","alpha")
    v_body <- c(
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
      sprintf("return tva_pr_log(R, S, D, %s, %s, %s, v);", if(is.null(parameters$t0)) "T" else "T - t0", t0_args, K_args)
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



  add_code(
    "functions",
    paste0("real log_lik_single(",paste(c(datsig(names_back = l_data), parsig(l_pars)), collapse = ", "),") {"),
    paste0("\t", l_body),
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
        vapply(names(parameters), function(name) {
          item <- parameters[[name]]
          if(is.null(item$prior)) sprintf("// no prior for %s", name) else sprintf("%s ~ %s;", name, as.character(as.expression(item$prior[[length(item$prior)]])))
        }, character(1))
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


  add_code(
    "model",
    "// likelihood (only if prior != 0)",
    "if(target() != negative_infinity()) {",
    if(isTRUE(parallel)) {
      c(
        paste0("\t",parmap(l_pars)),
        paste0("\ttarget += map_rect(log_lik_rect, phi, theta, x_r, x_i);")
      )
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
        c(
          "{",
          paste0("\t",parmap(l_pars)),
          paste0("\tlog_lik = map_rect(log_lik_rect, phi, theta, x_r, x_i);"),
          "}"
        )
      } else {
        paste0("for(i in 1:N) log_lik[i] = log_lik_single(",paste(c(datsig(names_back = l_data, types = FALSE, index = "i"), parsig(l_pars, types = FALSE, index = "i")),collapse=", "),");")
      }
    )
  }

  if(isTRUE(predict_scores)) {
    add_code(
      "functions",
      paste0("vector predict_score_single(",paste(c(datsig(names_back = p_data), parsig(p_pars)), collapse = ", "),") {"),
      paste0("\t", p_body),
      "}"
    )
    if(isTRUE(parallel)) {
      add_code(
        "functions",
        "vector predict_score_rect(vector phi, vector theta, data array[] real x_r, data array[] int x_i) {",
        paste0("\treturn predict_score_single(",paste(c(datremap(names_back = p_data), parremap(p_pars)),collapse=", "),");"),
        "}"
      )
    }
    add_code(
      "generated quantities",
      "// scores",
      sprintf("matrix[N,%d] pred_scores;", locations+1),
      if(isTRUE(parallel)) {
        c(
          "{",
          paste0("\t",parmap(l_pars)),
          paste0("\tpred_scores = to_matrix(map_rect(predict_score_rect, phi, theta, x_r, x_i), N, ",locations+1,", 0);"),
          "}"
        )
      } else {
        paste0("for(i in 1:N) pred_scores[i,] = to_row_vector(predict_score_single(",paste(c(datsig(names_back = p_data, types = FALSE, index = "i"), parsig(p_pars, types = FALSE, index = "i")),collapse=", "),"));")
      }
    )
  }

  if(isTRUE(simulate)) {
    add_code(
      "functions",
      paste0("array[] int response_rng(",paste(c(datsig(names_back = s_data), parsig(s_pars)), collapse = ", "),") {"),
      paste0("\t", s_body),
      "}"
    )
    add_code(
      "generated quantities",
      "// simulation",
      sprintf("array[N,%d] int Rsim;", locations),
      "for(i in 1:N) {",
      paste0("\tRsim[i,] = response_rng(",paste(c(datsig(names_back = s_data, types = FALSE, index = "i"), parsig(s_pars, types = FALSE, index = "i")),collapse=", "),");"),
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
      if(!is.null(data$D)) "\t\tif(D[i,j] && !S[i,j]) reject(\"Inconsistency detected: According to the data, there should be a distractor at location #\",j,\", of trial \",i,\", but there was no item displayed at that location!\");",
      if(!is.null(data$D)) "\t\tif(R[i,j] && D[i,j]) reject(\"Inconsistency detected: According to the data, the item at location #\",j,\", of trial \",i,\", was reported as a target but it was a distractor!\");",
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

  if(type == "stan2") {
    ret <- stanc(model_code = ret, isystem = stantva_path())$model_code
  } else if(type == "cpp") {
    ret <- stanc(model_code = ret, isystem = stantva_path())$cppcode
  }

  new("stantvacode", code = ret, config = call_args_list, include_path = stantva_path())
}

#'@export
stantvacode <- setClass("stantvacode", slots = c("code" = "character", "config" = "list", "include_path" = "character"))


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
tvadata <- setClass("tvadata", contains = "list")

#'@export
setMethod("summary", "tvadata", function(object, ...) {
  tva_report(object, ...)
})

#'@export
setMethod("show", "tvadata", function(object) {
  if(is.null(object$D)) {
    cat(col_cyan("TVA"), "data containing",object$N,"whole-report trial(s)\n")
  } else {
    cat(col_cyan("TVA"), "data containing",object$N,"whole- and/or partial-report trial(s)\n")
  }
  str(object)
})





#'@importClassesFrom rstan stanmodel
#'@export
stantvamodel <- setClass("stantvamodel", contains = "stanmodel", slots = c("code" = "stantvacode"))

#'@importClassesFrom rstan stanfit
#'@export
stantvafit <- setClass("stantvafit", contains = "stanfit", slots = c("stanmodel" = "stantvamodel", "data" = "list"))



#'@export
setMethod("show", c(object="stantvamodel"), function(object) {
  cat(col_cyan("StanTVA"), "model with following configuration:\n")
  for(cname in names(object@code@config)) {
    cat("  -", col_magenta(cname), "=", deparse(object@code@config[[cname]]),"\n")
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

#'@export
setMethod("simulate", c(object = "stantvamodel"), function(object, nsim, data, params, seed = NULL) {
  if(!isTRUE(object@code@config$simulate)) stop("StanTVA model must be compiled with `simulate` = TRUE in order to simulate responses!")
  data$R <- generate(object, data, params[rep(seq_len(nrow(params)), nsim),,drop=FALSE], "Rsim", seed)$Rsim
  data
})


#'@export
setMethod("simulate", c(object = "stantvafit"), function(object, nsim, newdata, seed = NULL) {
  if(!isTRUE(object@stanmodel@code@config$simulate)) stop("StanTVA model must be compiled with `simulate` = TRUE in order to simulate responses!")
  if(missing(newdata) || is.null(newdata)) {
    simulate(object = object@stanmodel, nsim = nsim, data = object@data, params = as.matrix(object), seed = seed)
  } else {
    simulate(object = object@stanmodel, nsim = 1L, data = newdata, params = as.matrix(object), seed = seed)
  }
})


setGeneric("fit", function(object, ...) {})

#'@export
setMethod("sampling", c(object = "stantvamodel"), function(object, data, ...) {
  stopifnot(inherits(data, "tvadata"))
  if(object@code@config$locations != data$locations) stop("Cannot fit a StanTVA model compiled for ",object@code@config$locations," location(s) to a data set with ",data$locations," location(s)!")
  f <- callNextMethod()
  f@stanmodel <- object
  f <- as(f, "stantvafit")
  f@data <- data
  f
})

#'@export
setMethod("optimizing", c(object = "stantvamodel"), function(object, data, ...) {
  stopifnot(inherits(data, "tvadata"))
  if(object@code@config$locations != data$locations) stop("Cannot fit a StanTVA model compiled for ",object@code@config$locations," location(s) to a data set with ",data$locations," location(s)!")
  callNextMethod()
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

#'@export
setMethod("predict", "stantvamodel", function(object, data, params) {
  if(!isTRUE(object@code@config$predict_scores)) stop("StanTVA model must be compiled with `predict_scores` = TRUE in order to predict scores!")
  generate(object, data, params, "pred_scores")$pred_scores
})

#'@export
setMethod("predict", "stantvafit", function(object, newdata) {
  if(!isTRUE(object@stanmodel@code@config$predict_scores)) stop("StanTVA model must be compiled with `predict_scores` = TRUE in order to predict scores!")
  if(missing(newdata)) extract(object, "pred_scores")$pred_scores
  else {
    predict(object@stanmodel, newdata, extract(object))
  }
})


#'@export
tva_report <- function(data) {
  tibble(
    condition = data$condition,
    exposure = data$T,
    score = as.integer(if(is.null(data$D)) rowSums(data$R == 1L & data$S == 1L) else rowSums(data$R == 1L & data$S == 1L & data$D == 0L)),
    n_items = as.integer(rowSums(data$S == 1L)),
    n_distractors = if(is.null(data$D)) integer(data$N) else as.integer(rowSums(data$D))
  ) %>% mutate(n_targets = n_items - n_distractors)
}


