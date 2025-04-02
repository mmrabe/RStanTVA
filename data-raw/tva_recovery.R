
library(MASS)
library(tidyverse)
library(Matrix)
library(gtools)


set.seed(12345)

mu_log_C <- log(70)
mu_log_C_highvslow <- log(120) - mu_log_C
mu_log_alpha <- log(0.7)
mu_log_alpha_highvslow <- 0
mu_log_pK <- c(0,0,1,0.5)
mu_log_pK_highvslow <- c(0, 0, 0, 0)
mu_mu0 <- 30
mu_mu0_highvslow <- 0
mu_sigma0 <- 2.5
mu_sigma0_highvslow <- 0
mu_log_w <- rnorm(3, 0, 0.1)
mu_log_w_highvslow <- c(0, 0, 0)

sd_subj <- c(
  log_C_Intercept = 0.2,
  log_C_condition = 0.05,
  log_alpha_Intercept = 0.2,
  log_alpha_condition = 0.05,
  log_pK_Intercept1 = 0.2,
  log_pK_condition1 = 0,
  log_pK_Intercept2 = 0.2,
  log_pK_condition2 = 0,
  log_pK_Intercept3 = 0.2,
  log_pK_condition3 = 0,
  log_pK_Intercept4 = 0.2,
  log_pK_condition4 = 0,
  sigma0_Intercept = 0.2,
  sigma0_condition = 0,
  mu0_Intercept = 10,
  mu0_condition = 0,
  log_w_Intercept1 = 0.2,
  log_w_condition1 = 0,
  log_w_Intercept2 = 0.2,
  log_w_condition2 = 0,
  log_w_Intercept3 = 0.2,
  log_w_condition3 = 0
)

cor_subj <- diag(length(sd_subj))
dimnames(cor_subj) <- list(names(sd_subj), names(sd_subj))
cor_subj["log_C_Intercept","log_C_condition"] <- -.2
cor_subj["log_C_condition","log_C_Intercept"] <- -.2
cor_subj["log_C_Intercept","log_alpha_Intercept"] <- -.3
cor_subj["log_alpha_Intercept","log_C_Intercept"] <- -.3

z_subj <- mvrnorm(50, rep(0, length(sd_subj)), cor_subj)

subject_ranefs <- z_subj * matrix(sd_subj, nrow = nrow(z_subj), ncol = length(sd_subj), byrow = TRUE)


tva_recovery <- crossing(
  bind_rows(
    crossing(
      subject = seq_len(nrow(subject_ranefs)),
      trial = sprintf("WR1%02d", 1:24),
      T = c(10,50,100,200),
      condition = factor(c("low","high"))
    ) %>% mutate(
      S = matrix(1L, nrow = n(), ncol = 6),
      D = matrix(0L, nrow = n(), ncol = 6)
    ),
    crossing(
      subject = seq_len(nrow(subject_ranefs)),
      trial = sprintf("WR2%02d", 1:12),
      T = c(10,50,100,200),
      condition = factor(c("low","high"))
    ) %>% mutate(
      S = t(vapply(seq_len(n()), function(i) sample(c(0L,0L,0L,0L,1L,1L)), integer(6))),
      D = matrix(0L, nrow = n(), ncol = 6)
    ),
    crossing(
      subject = seq_len(nrow(subject_ranefs)),
      trial = sprintf("PR%03d", 1:32),
      T = c(50,100,200),
      condition = factor(c("low","high"))
    ) %>% mutate(
      S = matrix(1L, nrow = n(), ncol = 6),
      D = t(vapply(seq_len(n()), function(i) sample(c(0L,0L,1L,1L,1L,1L)), integer(6)))
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
    C = exp((mu_log_C + subject_ranefs[j, "log_C_Intercept"]) + x_condition * (mu_log_C_highvslow + subject_ranefs[j, "log_C_condition"])),
    alpha = exp((mu_log_alpha + subject_ranefs[j, "log_alpha_Intercept"]) + x_condition * (mu_log_alpha_highvslow + subject_ranefs[j, "log_alpha_condition"])),
    w = vapply(1:3, function(i) exp((mu_log_w[i] + subject_ranefs[j, paste0("log_w_Intercept",i)]) + x_condition * (mu_log_w_highvslow[i] + subject_ranefs[j, paste0("log_w_condition",i)])), double(n())),
    mu0 = (mu_mu0 + subject_ranefs[j, "mu0_Intercept"]) + x_condition * (mu_mu0_highvslow + subject_ranefs[j, "mu0_condition"]),
    sigma0 = exp(mu_sigma0 + subject_ranefs[j, "sigma0_Intercept"]) + x_condition * (mu_sigma0_highvslow + subject_ranefs[j, "sigma0_condition"]),
    pK = vapply(1:4, function(i) exp((mu_log_pK[i] + subject_ranefs[j, paste0("log_pK_Intercept",i)]) + x_condition * (mu_log_pK_highvslow[i] + subject_ranefs[j, paste0("log_pK_condition",i)])), double(n()))
  ) %>%
  rename(subject = j) %>%
  dplyr::select(-x_condition) %>%
  mutate(across(c("w","pK"), ~cbind(.x,1)/(rowSums(.x)+1)))


tva_recovery$true_values <- bind_cols(tva_recovery$true_values, tva_recovery %>% dplyr::select(subject, condition) %>% left_join(subject_coefs, by = c("subject", "condition")) %>% ungroup() %>% dplyr::select(-subject,-condition))

for(i in seq_len(nrow(tva_recovery))) {

  C <- tva_recovery$true_values$C[i]
  alpha <- tva_recovery$true_values$alpha[i]
  w <- tva_recovery$true_values$w[i,] * if_else(tva_recovery$D[i,] == 1L, alpha, 1)
  mu0 <- tva_recovery$true_values$mu0[i]
  sigma0 <- tva_recovery$true_values$sigma0[i]
  pK <- tva_recovery$true_values$pK[i,]

  Ss <- which(tva_recovery$S[i,] == 1L)
  v <- C/1000 * w[Ss] / sum(w[Ss])
  K <- sample.int(length(pK), 1, prob = pK) - 1L
  t0 <- rnorm(1, mu0, sigma0)
  processing_times <- rexp(length(v), v) + t0
  Rs <- Ss[rank(processing_times) <= K & processing_times <= tva_recovery$T[i] & tva_recovery$D[i,] == 0L]
  tva_recovery$R[i,] <- as.integer(seq_along(w) %in% Rs)
  tva_recovery$true_values$t0[i] <- t0
  tva_recovery$true_values$K[i] <- K
  tva_recovery$true_values$v[i,] <- v[match(ncol(tva_recovery$S), Ss)]
  tva_recovery$true_values$t[i,] <- processing_times[match(ncol(tva_recovery$S), Ss)]
}

names(sd_subj) <- c("C_Intercept", "C_conditionhigh", "alpha_Intercept", "alpha_conditionhigh", "pK_Intercept[1]", "pK_conditionhigh[1]", "pK_Intercept[2]", "pK_conditionhigh[2]", "pK_Intercept[3]", "pK_conditionhigh[3]", "pK_Intercept[4]", "pK_conditionhigh[4]", "sigma0_Intercept", "sigma0_conditionhigh", "mu0_Intercept", "mu0_conditionhigh", "w_Intercept[1]", "w_conditionhigh[1]", "w_Intercept[2]", "w_conditionhigh[2]", "w_Intercept[3]", "w_conditionhigh[3]")
dimnames(cor_subj) <- list(names(sd_subj), names(sd_subj))
colnames(z_subj) <- names(sd_subj)

tva_recovery_true_params <- list(
  b = c(
    C_Intercept = mu_log_C,
    C_conditionhigh = mu_log_C_highvslow,
    alpha_Intercept = mu_log_alpha,
    alpha_conditionhigh = mu_log_alpha_highvslow,
    `pK_Intercept[1]` = mu_log_pK[1],
    `pK_conditionhigh[1]` = mu_log_pK_highvslow[1],
    `pK_Intercept[2]` = mu_log_pK[2],
    `pK_conditionhigh[2]` = mu_log_pK_highvslow[2],
    `pK_Intercept[3]` = mu_log_pK[3],
    `pK_conditionhigh[3]` = mu_log_pK_highvslow[3],
    `pK_Intercept[4]` = mu_log_pK[4],
    `pK_conditionhigh[4]` = mu_log_pK_highvslow[4],
    sigma0_Intercept = mu_sigma0,
    sigma0_conditionhigh = mu_sigma0_highvslow,
    mu0_Intercept = mu_mu0,
    mu0_conditionhigh = mu_mu0_highvslow,
    `w_Intercept[1]` = mu_log_w[1],
    `w_conditionhigh[1]` = mu_log_w_highvslow[1],
    `w_Intercept[2]` = mu_log_w[2],
    `w_conditionhigh[2]` = mu_log_w_highvslow[2],
    `w_Intercept[3]` = mu_log_w[3],
    `w_conditionhigh[3]` = mu_log_w_highvslow[3]
  ),
  s_subject = sd_subj,
  r_subject = cor_subj,
  z_subject = z_subj,
  coef_subject = subject_coefs
)




usethis::use_data(tva_recovery, tva_recovery_true_params, overwrite = TRUE, compress = "xz")
