library(tidyverse)
library(lme4)
library(sandwich)
library(lmtest)
library(lmerTest)
library(lfe)
library(truncnorm)
library(purrr)
library(simhelpers) # citation("simhelpers")
library(psych)
library(future)
library(furrr)
rm(list = ls())

# Data Generating Model----------------------------------------------------
generate_dat <- function(assumption, ES, J, n_bar, ICC_k, ICC_jk){
  
  # model parameter fixed
  ICC_j = 0.05        
  sparse = .1         
  X_m = 0
  X_sd = sqrt(50)
  K = J * 3.5
  
  # set sigma and random effects based on ICC
  tau_J00 = ICC_j  # School
  tau_K00 = ICC_k  # Neighborhood
  tau_JK0 = ICC_jk # Interaction
  sigma = sqrt(1-tau_J00-tau_K00-tau_JK0) 
  
  # coefficients
  gamma_w = ES*(sigma/X_sd)      
  gamma_b_j = ES*(sigma/X_sd)      
  gamma_b_k = ES*(sigma/X_sd)      
  gamma_b_jk = ES*(sigma/X_sd)  
  
  # data assignment
  dat <- 
    rerun(.n = J,  
          sample(1:(sparse * K), 
                 round(rtruncnorm(1, a = 1, mean = n_bar, sd = 0.2*n_bar), 0), 
                 replace = TRUE)) %>% 
    map_df(~ tibble(neighid = .x), .id = "schid") %>% 
    as.data.frame() %>% 
    mutate(
      schid = as.numeric(schid),
      neighid = round(schid * (K / J)) + neighid,
      neighid = ifelse(neighid > K, neighid - K, neighid)
    ) 
  
  # neighborhood data
  # create between-neighborhood variance of X 
  X_bw_neigh <- 
    dat %>% 
    group_by(neighid) %>% 
    summarise() %>% 
    mutate(X_bw_neigh = rnorm(nrow(.), mean = X_m, sd = X_sd))
  
  if (assumption == "exogeneity") {
    
    # neighborhood random effect when the exogeneity assumption is violated
    r <- 0.4 # correlation between X_bw_neigh and b_0j0
    
    neighbordata <-
      X_bw_neigh %>% 
      mutate(
        v_00k = rnorm(nrow(.), mean = 0, sd = sqrt((1 - r^2) * tau_K00)),
        c_00k = r * sqrt(tau_K00 / (X_sd ^ 2)) * X_bw_neigh + v_00k
      ) %>% 
      dplyr::select(-c(v_00k))
    
  } else {
    
    # neighborhood random effect when all assumptions are met:
    neighbordata <-
      X_bw_neigh %>% 
      mutate(
        c_00k = rnorm(nrow(.), mean = 0, sd = sqrt(tau_K00)),  
      ) 
  }
  
  dat <- dat %>% left_join(neighbordata, by = "neighid") 
  
  # school data
  # create between-school variance of X
  schooldata <- 
    dat %>% 
    group_by(schid) %>%
    summarise() %>%
    mutate(
      b_0j0 = rnorm(nrow(.), mean = 0, sd = sqrt(tau_J00)),
      X_bw_school = rnorm(nrow(.), mean = X_m, sd = X_sd)
    )
  
  dat <- dat %>% left_join(schooldata, by = "schid")
  
  # interaction data
  intdata <- dat %>% 
    group_by(neighid, schid) %>% 
    summarise(., .groups = "drop") %>% 
    mutate(
      X_bw_int = rnorm(nrow(.), mean = X_m, sd = X_sd),
      d_0jk = rnorm(nrow(.), mean = 0, sd = sqrt(tau_JK0))
    )
  
  dat <- dat %>% left_join(intdata, by = c("neighid", "schid"))
  
  # student data
  dat <-
    dat %>%
    mutate(
      # student ID
      stuid = 1:nrow(.),
      
      # student-level X
      X_within = rnorm(nrow(.), mean = X_m, sd = X_sd),
      X = X_within + X_bw_neigh + X_bw_school + X_bw_int,
      
      # student-level residuals e
      e = rnorm(nrow(.), mean = 0, sd = sigma),
      
      # outcome variable
      y = gamma_w * X_within + gamma_b_j * X_bw_school + gamma_b_k * X_bw_neigh + gamma_b_jk * X_bw_int + b_0j0 + c_00k + d_0jk + e
    )
  
  # calc cluster means
  cell_id <- dat %>% group_by(schid, neighid) %>% count() %>% ungroup() %>% 
    mutate(cellid = row_number()) %>% ungroup() %>% 
    dplyr::select(-n)
  
  dat <- 
    dat %>%
    left_join(cell_id, by = c("schid", "neighid")) %>%
    group_by(cellid) %>%
    mutate(cell_mean = mean(X)) %>%
    group_by(schid) %>%
    mutate(cluster_school = mean(X)) %>%
    group_by(neighid) %>%
    mutate(cluster_neigh = mean(X)) %>% 
    ungroup() %>%
    mutate(
      grand = mean(X),
      param = gamma_w)
  
  return(dat)
}

dat <- generate_dat(   
  # design factor
  assumption = "exogeneity",   # exogeneity or met
  ES = 0.5,
  J = 20,              # school j = {1..J}
  n_bar = 10,          # average number of students per school
  ICC_k = 0.05,         # neighbor ICC
  ICC_jk = 0.05         # neighbor ICC
)

# Model-fitting/Estimation-------------------------------------------------
cluster_bw <- function(dat) {
  
  X_fit <- felm(X ~ 0 | schid + neighid, data = dat)
  dat$X_adapt <- as.numeric(residuals(X_fit))
  fit <- lmer(y ~ X_adapt + cluster_neigh + cluster_school + (1|neighid) + (1|schid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      cov = recode(cov, "X_adapt" = "X",
                   "cluster_neigh" = "BW_N",
                   "cluster_school" = "BW_S"), 
      method = "cluster_bw",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)"))
  
  return(fixed_est)
}

cell_bw <- function(dat) {
  
  X_fit <- felm(X ~ 0 | cellid, data = dat)
  dat$X_adapt_cell <- as.numeric(residuals(X_fit))
  cellmean_fit <- felm(cell_mean ~ 0 | schid + neighid, data = dat)
  dat$cellmean_adapt <- as.numeric(residuals(cellmean_fit))
  fit <- lmer(y ~ X_adapt_cell + cluster_neigh + cluster_school + cellmean_adapt + (1|neighid) + (1|schid) + (1|cellid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      cov = recode(cov, "X_adapt_cell" = "X",
                   "cluster_neigh" = "BW_N",
                   "cluster_school" = "BW_S",
                   "cellmean_adapt" = "BW_C"), 
      method = "cell_bw",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)"))
  
  return(fixed_est)
}


# bind_results
estimate <- function(dat) {
  
  CI_level = .95
  crit <- qnorm((1 + CI_level) / 2)
  
  results <- bind_rows(
    # cluster_bw(dat),
    cell_bw(dat),
    ) %>% 
    as_tibble() 
  return(results)
}

cluster_bw(dat)
cell_bw(dat)
estimate(dat = dat) 

results <-
  rerun(5, {
    dat <- generate_dat(   
      assumption = "exogeneity",   # exogeneity or met
      ES = 0.5,
      J = 20,              # school j = {1..J}
      n_bar = 30,          # average number of students per school
      ICC_k = 0.05,         # neighbor ICC
      ICC_jk = 0.01         # neighbor ICC
    ) %>%
    estimate()
  }) %>%
  bind_rows() %>% 
  filter(converged == TRUE) %>% 
  group_by(method, cov) %>% 
  summarise(mean = mean(est),
            se = mean(se),
            n = n())
results

# Simulation driver -------------------------------------------------------
run_sim <- function(iterations,
                    assumption, ES, J, n_bar, ICC_k, ICC_jk, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  options(dplyr.summarise.inform = FALSE)

  results <-
    rerun(iterations, {
      generate_dat(
        assumption = assumption,
        ES = ES,
        J = J,
        n_bar = n_bar,
        ICC_k = ICC_k,
        ICC_jk = ICC_jk
      ) %>%
        estimate()
    }) %>%
    bind_rows() 
}

results <- 
  run_sim(iterations = 1,
          assumption = "met",   # exogeneity or met
          ES = 0.5,
          J = 30,             # school j = {1..J}
          n_bar = 30,           # average number of students per school
          ICC_k = 0.15,         # neighbor ICC
          ICC_jk = 0.01,        # neighbor ICC
          seed = NULL)
results

# simulation parameters as vectors/lists
design_factors <- list(
  assumption = c("met", "exogeneity"),
  ES = c(0.01, 0.03, 0.05),
  J = c(200),  # J = school
  n_bar = 100,
  # n_bar = c(30, 100), 
  ICC_k = c(0.15, 0.25),
  # ICC_k = c(0.05, 0.15, 0.25), 
  ICC_jk = c(0)
  # ICC_jk = c(0, 0.01, 0.05) 
)

params <- 
  cross_df(design_factors) %>%
  mutate(
    iterations = 6000, 
    seed = 20230129 + 1:n()
  )

nrow(params)
head(params)

# pmap
options(error=recover)
plan("multisession") 

system.time(
  results <-
    params %>%
    mutate(res = future_pmap(., .f = run_sim, .options = furrr_options(seed=NULL))) %>% 
    unnest(cols = res)
)

results_final <- bind_rows(results_final, results)

# Check -------------------------------------------------------------------
# see if the relative parameter bias and relative bias of SE fit.
str(results_final)
head(results_final)
names(results_final)
nrow(results_final)

rel_crit <- results_final %>%
  filter(method == "cluster_bw" & cov == "X") %>% 
  mutate(param = ES*(sqrt(1-0.05-ICC_k-ICC_jk)/sqrt(50))) %>% 
  group_by(method, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>%
  group_modify(~ calc_relative(.x, estimates = est, true_param = param))

rel_crit %>%
  ggplot(aes(x = J, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.95, 1.05), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ES, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_crit <- results_final %>%
  filter(method == "cell_bw" & cov == "X") %>% 
  mutate(param = ES*(sqrt(1-0.05-ICC_k-ICC_jk)/sqrt(50))) %>% 
  group_by(method, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>%
  group_modify(~ calc_relative(.x, estimates = est, true_param = param))

rel_crit %>%
  ggplot(aes(x = J, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.95, 1.05), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ES, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 


# EScalc ------------------------------------------------------------------
cluster_BW_N <- results_final %>% 
  filter(method == "cluster_bw", cov == "BW_N" & converged == TRUE) %>%
  group_by(assumption, ES, n_bar, ICC_k, ICC_jk) %>% 
  summarise(n = n(),
            param_est = mean(est)) %>% ungroup() %>% 
  mutate(method = "cluster_bw", cov = "BW_N") %>% 
  select(method, cov, everything())

cluster_BW_S <- results_final %>% 
  filter(method == "cluster_bw", cov == "BW_S" & converged == TRUE) %>%
  group_by(assumption, ES, n_bar, ICC_k, ICC_jk) %>% 
  summarise(n = n(),
            param_est = mean(est)) %>% ungroup() %>% 
  mutate(method = "cluster_bw", cov = "BW_S") %>% 
  select(method, cov, everything())


cell_BW_N <- results_final %>% 
  filter(method == "cell_bw", cov == "BW_N" & converged == TRUE) %>%
  group_by(assumption, ES, n_bar, ICC_k, ICC_jk) %>% 
  summarise(n = n(),
            param_est = mean(est)) %>% ungroup() %>% 
  mutate(method = "cell_bw", cov = "BW_N") %>% 
  select(method, cov, everything())


cell_BW_S <- results_final %>% 
  filter(method == "cell_bw", cov == "BW_S" & converged == TRUE) %>%
  group_by(assumption, ES, n_bar, ICC_k, ICC_jk) %>% 
  summarise(n = n(),
            param_est = mean(est)) %>% ungroup() %>% 
  mutate(method = "cell_bw", cov = "BW_S") %>% 
  select(method, cov, everything())


cell_BW_C <- results_final %>% 
  filter(method == "cell_bw", cov == "BW_C" & converged == TRUE) %>%
  group_by(assumption, ES, n_bar, ICC_k, ICC_jk) %>% 
  summarise(n = n(),
            param_est = mean(est)) %>% ungroup() %>% 
  mutate(method = "cell_bw", cov = "BW_C") %>% 
  select(method, cov, everything()) 

param_list <- bind_rows(cluster_BW_N, cluster_BW_S,
                        cell_BW_N, cell_BW_S, cell_BW_C) %>% 
  select(-c(n))

save(results_final, file = "EScalc_raw.RData")
save(param_list, file = "EScalc_param.RData")

load("EScalc_raw.RData")
