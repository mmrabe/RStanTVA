#'@importFrom ggplot2 aes facet_wrap geom_line ggplot labs
#'@importFrom dplyr rowwise
#'@importFrom rlang .data .env
#'
NULL

#' Compute TVA processing rates
#'
#' Computes processing rates \code{v} for given stimulus (\code{S}), distractor (\code{D}),
#' and exposure durations (\code{T}) using a fitted TVA model.
#'
#' @param S Integer matrix indicating stimulus presence.
#' @param D Integer matrix indicating distractor presence.
#' @param T Numeric vector of exposure durations.
#' @param tva_model TVA model object.
#' @param tva_fit Fitted TVA model parameters.
#'
#' @return A matrix of processing rates.
#' @export
tva_processing_rates <- function(S = matrix(1L, 1L, tva_model@code@config$locations), D = matrix(0L, nrow(S), ncol(S)), tva_model, tva_fit) {
  calc_v <- Vectorize(tva_model@initializers$calculate_v, c("nS","S","D"))
  tva_data <- tibble(S, D)
  calc_v_args <- list(nS = rowSums(tva_data$S))
  for(arg in setdiff(names(formals(tva_model@initializers$calculate_v)), c("pstream__","rng__","nS"))) {
    if(arg %in% names(tva_data)) {
      if(is.matrix(tva_data[[arg]])) {
        calc_v_args[[arg]] <- lapply(seq_len(nrow(tva_data[[arg]])), function(i) tva_data[[arg]][i,])
      } else {
        calc_v_args[[arg]] <- tva_data[[arg]]
      }
    } else {
      par_names <- names(tva_fit$par) == arg | startsWith(names(tva_fit$par), paste0(arg,"["))
      calc_v_args[[arg]] <- tva_fit$par[par_names]
    }
  }
  t(do.call(calc_v, calc_v_args))
}

#' Predict expected score
#'
#' Computes the expected score based on TVA score probabilities.
#'
#' @inheritParams tva_processing_rates
#'
#' @return Numeric vector of expected scores.
#' @export
tva_predict_score <- function(S = matrix(1L, tva_model@code@config$locations, length(T)), D = matrix(0L, tva_model@code@config$locations, length(T)), T, tva_model, tva_fit, scores = 1:tva_model@code@config$locations) {
  as.vector(tva_score_prob(S,D,T,tva_model,tva_fit,scores,FALSE) %*% scores)
}

#' Compute score probabilities
#'
#' Calculates the probability of observing specific scores.
#'
#' @inheritParams tva_processing_rates
#' @param score Integer vector of scores.
#' @param log_p Logical; return log-probabilities if TRUE.
#'
#' @return Matrix of probabilities (or log-probabilities).
#' @export
tva_score_prob <- function(S = matrix(1L, tva_model@code@config$locations, length(T)), D = matrix(0L, tva_model@code@config$locations, length(T)), T, tva_model, tva_fit, score = 0:ncol(S), log_p = FALSE) {
  ret <- do.call(cbind, lapply(score, function(score) {
    pred_score <- Vectorize(tva_model@initializers$tva_pr_score_log, c("S","D","t","v"))
    v <- tva_processing_rates(S, D, tva_model, tva_fit)
    pred_score(
      S = lapply(seq_len(nrow(S)), function(i) S[i,]),
      D = lapply(seq_len(nrow(D)), function(i) D[i,]),
      t = if("t0" %in% names(tva_fit$par)) T - tva_fit$par["t0"] else T,
      v = lapply(seq_len(nrow(v)), function(i) v[i,]),
      score = score,
      K_args = tva_fit$par[list("bernoulli" = "K", "free" = sprintf("pK[%d]", seq_len(tva_model@code@config$locations+1L)), "betabinomial" = c("aK","bK"), "binomial" = "pK", "hypergeometric" = c("gK","bK"), "probit" = c("mK","sK"))[[tva_model@code@config$K_mode]]],
      t0_args = tva_fit$par[list("constant" = character(0), "gaussian" = c("mu0","sigma0"), "exponential" = "mu0")[[tva_model@code@config$t0_mode]]]
    )
  }))
  colnames(ret) <- sprintf("score=%d",score)
  if(!isTRUE(log_p)) exp(ret) else ret
}

#' Simulate TVA response
#'
#' Simulates trial-wise TVA responses.
#'
#' @inheritParams tva_processing_rates
#' @param seed Random seed.
#'
#' @return Matrix of simulated responses.
#' @export
tva_simulate_response <- function(S = matrix(1L, tva_model@code@config$locations, length(T)), D = matrix(0L, tva_model@code@config$locations, length(T)), T, tva_model, tva_fit, seed = sample.int(.Machine$integer.max, 1)) {
  sim_trial <- Vectorize(tva_model@initializers$tva_pr_rng, c("S","D","t","v"))
  v <- tva_processing_rates(S, D, tva_model, tva_fit)
  t(sim_trial(
    S = lapply(seq_len(nrow(S)), function(i) S[i,]),
    D = lapply(seq_len(nrow(D)), function(i) D[i,]),
    t = if("t0" %in% names(tva_fit$par)) T - tva_fit$par["t0"] else T,
    v = lapply(seq_len(nrow(v)), function(i) v[i,]),
    base_rng__ = get_rng(seed),
    K_args = tva_fit$par[list("bernoulli" = "K", "free" = sprintf("pK[%d]", seq_len(tva_model@code@config$locations+1L)), "betabinomial" = c("aK","bK"), "binomial" = "pK", "hypergeometric" = c("gK","bK"), "probit" = c("mK","sK"))[[tva_model@code@config$K_mode]]],
    t0_args = tva_fit$par[list("constant" = character(0), "gaussian" = c("mu0","sigma0"), "exponential" = "mu0")[[tva_model@code@config$t0_mode]]]
  ))
}

#' Simulate TVA score
#'
#' Simulates scores by summing simulated responses.
#'
#' @inheritParams tva_simulate_response
#'
#' @return Numeric vector of simulated scores.
#' @export
tva_simulate_score <- function(S = matrix(1L, tva_model@code@config$locations, length(T)), D = matrix(0L, tva_model@code@config$locations, length(T)), T, tva_model, tva_fit, seed = sample.int(.Machine$integer.max, 1)) {
  rowSums(tva_simulate_response(S,D,T,tva_model,tva_fit,seed))
}

#' Compute response probability
#'
#' Computes the probability of a response pattern \code{R}.
#'
#' @param R Integer matrix of responses.
#' @inheritParams tva_processing_rates
#' @param log_p Logical; return log-probabilities if TRUE.
#'
#' @return Numeric vector of probabilities.
#' @export
tva_response_prob <- function(R, S = matrix(1L, tva_model@code@config$locations, length(T)), D = matrix(0L, tva_model@code@config$locations, length(T)), T, tva_model, tva_fit, log_p = FALSE) {
  prob_trial <- Vectorize(tva_model@initializers$tva_pr_log, c("S","D","R","t","v"))
  v <- tva_processing_rates(S, D, tva_model, tva_fit)
  ret <- prob_trial(
    R = lapply(seq_len(nrow(R)), function(i) R[i,]),
    S = lapply(seq_len(nrow(S)), function(i) S[i,]),
    D = lapply(seq_len(nrow(D)), function(i) D[i,]),
    t = if("t0" %in% names(tva_fit$par)) T - tva_fit$par["t0"] else T,
    v = lapply(seq_len(nrow(v)), function(i) v[i,]),
    K_args = tva_fit$par[list("bernoulli" = "K", "free" = sprintf("pK[%d]", seq_len(tva_model@code@config$locations+1L)), "betabinomial" = c("aK","bK"), "binomial" = "pK", "hypergeometric" = c("gK","bK"), "probit" = c("mK","sK"))[[tva_model@code@config$K_mode]]],
    t0_args = tva_fit$par[list("constant" = character(0), "gaussian" = c("mu0","sigma0"), "exponential" = "mu0")[[tva_model@code@config$t0_mode]]]
  )
  if(isTRUE(log_p)) ret else exp(ret)
}

tva_generate_matrix <- function(n, k) {
  t(combn(n, k, FUN = function(i) {
    r <- integer(n)
    r[i] <- 1L
    r
  }))
}

tva_generate_matrices <- function(locations, nS, nD) {
  stopifnot(all(nS <= locations), all(nD <= nS), length(nS) == length(nD))
  lapply(seq_along(nS), function(i) {
    S <- tva_generate_matrix(locations, nS[i])
    D <- tva_generate_matrix(nS[i], nD[i])
    crossing(S, D) %>% rowwise() %>% mutate(D = {
      L <- integer(locations)
      L[which(as.logical(as.vector(S)))[as.logical(D)]] <- 1L
      t(L)
    })
  }) %>% bind_rows()
}

#' Integrate TVA function over configurations
#'
#' Applies a function over all stimulus/distractor configurations and averages results.
#'
#' @param fun Function to evaluate.
#' @param locations Number of locations.
#' @param nS Number of stimuli.
#' @param nD Number of distractors.
#' @param T Exposure durations.
#' @param ... Additional arguments passed to \code{fun}.
#'
#' @return Numeric vector or matrix of averaged values.
#' @export
tva_integrate <- function(fun, locations = max(nS), nS, nD, T, ...) t(Vectorize(function(fun, locations, nS, nD, T, ...) {
  M <- tva_generate_matrices(locations, nS, nD)
  Y <- fun(S = M$S, D = M$D, T = T, ...)
  if(is.matrix(Y)) colMeans(Y) else mean(Y)
}, c("nS","nD","T"))(fun, locations, nS, nD, T, ...))


#' Compute item-level probabilities
#'
#' Computes probability of correctly reporting each item.
#'
#' @inheritParams tva_processing_rates
#' @param item Items to evaluate.
#' @param log_p Logical; return log-probabilities if TRUE.
#'
#' @return Matrix of probabilities.
#' @export
tva_item_prob <- function(S = matrix(1L, tva_model@code@config$locations, length(T)), D = matrix(0L, tva_model@code@config$locations, length(T)), T, tva_model, tva_fit, item = seq_len(ncol(S)), log_p = FALSE) {
  prob_trial <- Vectorize(tva_model@initializers$tva_pr_log, c("S","D","R","t","v"))
  v <- tva_processing_rates(S, D, tva_model, tva_fit)
  ret <- vapply(item, function(item) {
    prob_trial(
      R = lapply(seq_len(nrow(S)), function(i) as.integer(seq_len(ncol(S)) == item)),
      S = lapply(seq_len(nrow(S)), function(i) S[i,]),
      D = lapply(seq_len(nrow(S)), function(i) as.integer(D[i,] | seq_len(ncol(S)) != item)),
      t = if("t0" %in% names(tva_fit$par)) T - tva_fit$par["t0"] else T,
      v = lapply(seq_len(nrow(v)), function(i) v[i,]),
      K_args = tva_fit$par[list("bernoulli" = "K", "free" = sprintf("pK[%d]", seq_len(tva_model@code@config$locations+1L)), "betabinomial" = c("aK","bK"), "binomial" = "pK", "hypergeometric" = c("gK","bK"), "probit" = c("mK","sK"))[[tva_model@code@config$K_mode]]],
      t0_args = tva_fit$par[list("constant" = character(0), "gaussian" = c("mu0","sigma0"), "exponential" = "mu0")[[tva_model@code@config$t0_mode]]]
    )
  }, double(nrow(S)))
  colnames(ret) <- sprintf("Location %d", item)
  if(isTRUE(log_p)) ret else exp(ret)
}


#' Plot score densities
#'
#' Visualizes score probability densities across exposure durations.
#'
#' @param nS Number of stimuli.
#' @param nD Number of distractors.
#' @param tva_model TVA model.
#' @param tva_fit Fitted parameters.
#' @param locations Number of locations.
#' @param scores Scores to plot.
#' @param xlim Range of exposure durations.
#' @param n Number of grid points.
#'
#' @return A ggplot object.
#' @export
tva_plot_score_densities <- function(nS, nD, tva_model, tva_fit, scores = 0:tva_model@code@config$locations, xlim = c(0,200), n = 101) {
  T <- seq(xlim[1], xlim[2], length.out = n)
  panels <- tibble(nS, nD)
  P_scores <- bind_rows(lapply(seq_len(nrow(panels)), function(j) {
    P <- tva_integrate(fun=tva_score_prob,locations=tva_model@code@config$locations,nS=panels$nS[j],nD=panels$nD[j],T=T,tva_model=tva_model,tva_fit=tva_fit,score=scores)
    bind_cols(panels[j,], T = rep(T,length(scores)), Score = ordered(rep(scores, each=length(T))), y = as.vector(P))
  }))
  ggplot(P_scores) +
    geom_line(aes(x=.data$T,y=.data$y,color=.data$Score,group=.data$Score)) +
    facet_wrap(~sprintf("%dT%dD",.data$nS-.data$nD,.data$nD)) +
    labs(x = "Exposure duration", y = "Probability density")
}


#' Plot expected scores
#'
#' Visualizes expected score as a function of exposure duration.
#'
#' @inheritParams tva_plot_score_densities
#'
#' @return A ggplot object.
#' @export
tva_plot_predicted_score <- function(nS, nD, tva_model, tva_fit, scores = 1:tva_model@code@config$locations, xlim = c(0,200), n = 101) {
  T <- seq(xlim[1], xlim[2], length.out = n)
  panels <- tibble(nS, nD)
  P_scores <- bind_rows(lapply(seq_len(nrow(panels)), function(j) {
    P <- tva_integrate(fun=tva_predict_score,locations=tva_model@code@config$locations,scores=scores,nS=panels$nS[j],nD=panels$nD[j],T=T,tva_model=tva_model,tva_fit=tva_fit)
    bind_cols(panels[j,], T = T, y = as.vector(P))
  })) %>% mutate(Display = sprintf("%dT%dD",.data$nS-.data$nD,.data$nD))
  ggplot(P_scores) +
    geom_line(aes(x=.data$T,y=.data$y,color=.data$Display,group=.data$Display)) +
    labs(x = "Exposure duration", y = "Expected score")
}

