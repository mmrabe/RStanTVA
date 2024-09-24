## code to prepare `tva_recovery` dataset goes here


library(tidyverse)
library(TruncExpFam)
library(gtools)

set.seed(349593)
Nprior <- 100
N <- 400

max_nS <- 6L


conditions <- crossing(exposure = c(10,50,80,100,150,200,250,400,500,1000), nS = seq(2,max_nS,2), type = c("WR","PR"))


tva_recovery <- bind_rows(
  tibble(
    t0_mode = "constant",
    pK = rdirichlet(Nprior, rep(1, max_nS+1)),
    mu0 = runif(Nprior, -100, 100),
    sigma0 = NA_real_,
    C = as.double(rtruncgamma(Nprior, 2, 0.02, a=1, b=200)),
    w = rdirichlet(Nprior, rep(1, max_nS)),
    alpha = runif(Nprior, 0, 2)
  ),
  tibble(
    t0_mode = "gaussian",
    pK = rdirichlet(Nprior, rep(1, max_nS+1)),
    mu0 = runif(Nprior, -100, 100),
    sigma0 = runif(Nprior, 0, 50),
    C = as.double(rtruncgamma(Nprior, 2, 0.02, a=1, b=200)),
    w = rdirichlet(Nprior, rep(1, max_nS)),
    alpha = runif(Nprior, 0, 2)
  )
)

for(i in seq_len(nrow(tva_recovery))) {

  trials <- bind_rows(lapply(seq_len(N), \(k) {

    j <- sample.int(nrow(conditions), 1)

    # sample experimental exposure duration
    exposure <- conditions$exposure[j]

    # sample experimental stimuli
    nS <- conditions$nS[j]
    nD <- if(conditions$type[j]=="WR") 0L else nS %/% 2L
    S <- logical(max_nS)
    S[sample.int(max_nS, nS)] <- TRUE
    D <- logical(max_nS)
    D[which(S)[sample.int(nS, nD)]] <- TRUE

    stopifnot(sum(S) == nS, sum(D) == nD, !any(D & !S))

    # sample memory capacity
    K <- sample.int(max_nS+1, 1, prob = tva_recovery$pK[i,])-1L

    stopifnot(K <= max_nS)

    # calculate theoretical processing rates
    trial_w <- tva_recovery$w[i,S] * if_else(D[S], tva_recovery$alpha[i], 1)
    trial_w <- trial_w/sum(trial_w)
    v <- tva_recovery$C[i]/1000 * trial_w

    # generate trial t0
    t0 <- if(tva_recovery$t0_mode[i] == "constant") tva_recovery$mu0[i] else rnorm(1, tva_recovery$mu0[i], tva_recovery$sigma0[i])

    # generate processing times
    processing_times <- rexp(nS, v) + t0
    #processing_times <- vapply(v, function(v) mod$functions$tva_t_rng(v/1000, c(mu, sigma)), double(1))

    # items admitted to memory (R)
    R <- logical(max_nS)
    R[S] <- rank(processing_times) <= K & processing_times <= exposure & !D[S]

    stopifnot(!any(R & !S), !any(R & D), sum(R) <= K)

    tibble(
      condition = j,
      R = t(R) + 0L,
      S = t(S) + 0L,
      D = t(D) + 0L,
      exposure = exposure
    )
  }))

  f <- file(sprintf("inst/extdata/recovery/recovery_%d.dat", i), "w")
  writeLines(as.character(N), f)

  stimuli <- replicate(N, sample(LETTERS, max_nS))

  Tstr <- vapply(seq_len(N), \(i) paste0(if_else(trials$S[i,] & !trials$D[i,], stimuli[,i], "0"), collapse=""), "")
  Dstr <- vapply(seq_len(N), \(i) paste0(if_else(trials$S[i,] & trials$D[i,], stimuli[,i], "0"), collapse=""), "")
  Rstr <- vapply(seq_len(N), \(i) if(any(trials$R[i,])) paste0(stimuli[as.logical(trials$R[i,]),i], collapse="") else "-", "")

  writeLines(sprintf("%d\t%.0f\t%s\t%s\t%s", trials$condition, trials$exposure, Tstr, Dstr, Rstr), f)
  close(f)

  c(list(N=N,max_nS=max_nS), as.list(trials %>% select(-condition)))


}





usethis::use_data(tva_recovery, overwrite = TRUE, compress = "xz")
