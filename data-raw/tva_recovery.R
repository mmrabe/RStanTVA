
library(MASS)
library(tidyverse)
library(Matrix)
library(gtools)


set.seed(12345)

b <- c(
  C_Intercept = log(80),
  C_conditionhigh = log(100) - log(80),
  alpha_Intercept = log(0.7),
  alpha_conditionhigh = 0,
  `mK_Intercept` = 3.5,
  `mK_conditionhigh` = 0,
  `sK_Intercept` = -0.7,
  `sK_conditionhigh` = 0,
  sigma0_Intercept = 2.0,
  sigma0_conditionhigh = 0,
  mu0_Intercept = 20,
  mu0_conditionhigh = 0,
  `w_Intercept[1]` = 0,
  `w_conditionhigh[1]` = 0
)

sd_subj <- c(
  C_Intercept = 0.2,
  C_conditionhigh = 0.05,
  alpha_Intercept = 0.2,
  alpha_conditionhigh = 0.05,
  `mK_Intercept` = 0.75,
  `mK_conditionhigh` = 0,
  `sK_Intercept` = 0.1,
  `sK_conditionhigh` = 0,
  sigma0_Intercept = 0.2,
  sigma0_conditionhigh = 0,
  mu0_Intercept = 10,
  mu0_conditionhigh = 0,
  `w_Intercept[1]` = 0.1,
  `w_conditionhigh[1]` = 0
)

cor_subj <- diag(length(sd_subj))
dimnames(cor_subj) <- list(names(sd_subj), names(sd_subj))
cor_subj["C_Intercept","C_conditionhigh"] <- cor_subj["C_conditionhigh","C_Intercept"] <- -.2
cor_subj["C_Intercept","alpha_Intercept"] <- cor_subj["alpha_Intercept","C_Intercept"] <- -.3
cor_subj["mK_Intercept","sK_Intercept"] <- cor_subj["mK_Intercept","sK_Intercept"] <- 0.2

z_subj <- mvrnorm(50, rep(0, length(sd_subj)), cor_subj)

subject_ranefs <- z_subj * matrix(sd_subj, nrow = nrow(z_subj), ncol = length(sd_subj), byrow = TRUE)


tva_recovery <- crossing(
  bind_rows(
    crossing(
      subject = seq_len(nrow(subject_ranefs)),
      trial = 1:13,
      T = c(16.6,33.3,50,100,150,200),
      condition = factor(c("low","high")),
      nD = 0L
    ) %>% mutate(
      S = matrix(1L, nrow = n(), ncol = 6),
      D = matrix(0L, nrow = n(), ncol = 6)
    ),
    crossing(
      subject = seq_len(nrow(subject_ranefs)),
      trial = 1:13,
      T = 150,
      nD = 2:4,
      condition = factor(c("low","high"))
    ) %>% mutate(
      S = matrix(1L, nrow = n(), ncol = 6),
      D = t(vapply(nD, function(nD) sample(rep(0:1,c(6-nD,nD))), integer(6)))
    )
  )
) %>% mutate(
  R = matrix(NA_integer_, nrow = n(), ncol = 6)
) %>% group_by(subject) %>% sample_n(n()) %>% mutate(trial = seq_len(n()), true_values = tibble(t = matrix(NA_real_, n(), 6), v = matrix(NA_real_, n(), 6), t0 = NA_real_, K = NA_integer_)) %>% arrange(subject, trial)

contrasts(tva_recovery$condition) <- contr.treatment(levels(tva_recovery$condition), base = match("low", levels(tva_recovery$condition)))

subject_coefs <- crossing(
  j = 1:50,
  condition = unique(tva_recovery$condition)
) %>%
  mutate(
    x_condition = contrasts(tva_recovery$condition)[condition,1],
    C = exp((b["C_Intercept"] + subject_ranefs[j, "C_Intercept"]) + x_condition * (b["C_conditionhigh"] + subject_ranefs[j, "C_conditionhigh"])),
    alpha = exp((b["alpha_Intercept"] + subject_ranefs[j, "alpha_Intercept"]) + x_condition * (b["alpha_conditionhigh"] + subject_ranefs[j, "alpha_conditionhigh"])),
    w = vapply(1, function(i) exp((b[sprintf("w_Intercept[%d]",i)] + subject_ranefs[j, sprintf("w_Intercept[%d]",i)]) + x_condition * (b[sprintf("w_conditionhigh[%d]",i)] + subject_ranefs[j, sprintf("w_conditionhigh[%d]",i)])), double(n())),
    mu0 = (b["mu0_Intercept"] + subject_ranefs[j, "mu0_Intercept"]) + x_condition * (b["mu0_conditionhigh"] + subject_ranefs[j, "mu0_conditionhigh"]),
    sigma0 = exp(b["sigma0_Intercept"] + subject_ranefs[j, "sigma0_Intercept"]) + x_condition * (b["sigma0_conditionhigh"] + subject_ranefs[j, "sigma0_conditionhigh"]),
    mK = (b["mK_Intercept"] + subject_ranefs[j, "mK_Intercept"]) + x_condition * (b["mK_conditionhigh"] + subject_ranefs[j, "mK_conditionhigh"]),
    sK = exp((b["sK_Intercept"] + subject_ranefs[j, "sK_Intercept"]) + x_condition * (b["sK_conditionhigh"] + subject_ranefs[j, "sK_conditionhigh"])),
  ) %>%
  rename(subject = j) %>%
  dplyr::select(-x_condition) %>%
  mutate(across("w", ~cbind(.x,1)/(rowSums(.x)+1)))


tva_recovery$true_values <- bind_cols(tva_recovery$true_values, tva_recovery %>% dplyr::select(subject, condition) %>% left_join(subject_coefs, by = c("subject", "condition")) %>% ungroup() %>% dplyr::select(-subject,-condition))

for(i in seq_len(nrow(tva_recovery))) {

  C <- tva_recovery$true_values$C[i]
  alpha <- tva_recovery$true_values$alpha[i]
  w <- rep(tva_recovery$true_values$w[i,], c(3,3)) * if_else(!!tva_recovery$D[i,], alpha, 1)
  w <- w / sum(w)
  mu0 <- tva_recovery$true_values$mu0[i]
  sigma0 <- tva_recovery$true_values$sigma0[i]
  mK <- tva_recovery$true_values$mK[i]
  sK <- tva_recovery$true_values$sK[i]

  Ss <- which(tva_recovery$S[i,] == 1L)
  v <- C/1000 * w[Ss] / sum(w[Ss])
  K <- floor(rnorm(1, mK, sK))
  if(K < 0) K <- 0
  if(K > ncol(tva_recovery$S)) K <- ncol(tva_recovery$S)
  t0 <- rnorm(1, mu0, sigma0)
  processing_times <- rexp(length(v), v) + t0
  Rs <- Ss[rank(processing_times) <= K & processing_times <= tva_recovery$T[i] & !tva_recovery$D[i,]]
  tva_recovery$R[i,] <- as.integer(seq_along(w) %in% Rs)
  tva_recovery$true_values$t0[i] <- t0
  tva_recovery$true_values$K[i] <- K
  tva_recovery$true_values$v[i,] <- v[match(ncol(tva_recovery$S), Ss)]
  tva_recovery$true_values$t[i,] <- processing_times[match(ncol(tva_recovery$S), Ss)]
}

dimnames(cor_subj) <- list(names(sd_subj), names(sd_subj))
colnames(z_subj) <- names(sd_subj)

tva_recovery_true_params <- list(
  b = b,
  s_subject = sd_subj,
  r_subject = cor_subj,
  z_subject = z_subj,
  coef_subject = subject_coefs
)




usethis::use_data(tva_recovery, tva_recovery_true_params, overwrite = TRUE, compress = "xz")

