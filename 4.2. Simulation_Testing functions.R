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
source("01. Simulation functions.R")

# Data Generating Model----------------------------------------------------
dat <- generate_dat(   
  # design factor
  assumption = "exogeneity",   # exogeneity or met
  ES = 0.5,
  J = 20,              # school j = {1..J}
  n_bar = 10,          # average number of students per school
  ICC_k = 0.05,         # neighbor ICC
  ICC_jk = 0.05         # neighbor ICC
) 

# Data structure-----------------------------------------------------------
dat %>% group_by(schid) %>% count() %>% ungroup() %>% 
  summarise(sample_size = mean(n)) # should match n_bar
dat %>% group_by(schid, neighid) %>% count() %>% 
  pivot_wider(names_from = "neighid", values_from = "n")

# Model-fitting/Estimation-------------------------------------------------
uncentered(dat)
grand(dat)
cluster_w(dat)
cluster_bw(dat)
cell_w(dat)
cell_bw(dat)
fe_crve(dat)
hybrid_sch_fe(dat)
hybrid_sch_re(dat)
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
  bind_rows() 

names(results)
results

# Simulation driver -------------------------------------------------------
results <- 
  run_sim(iterations = 5,
          assumption = "met",   # exogeneity or met
          ES = 0.5,
          J = 20,              # school j = {1..J}
          n_bar = 30,          # average number of students per school
          ICC_k = 0.15,         # neighbor ICC
          ICC_jk = 0.01,         # neighbor ICC
          seed = NULL)

# simulation parameters as vectors/lists
design_factors <- list(
  assumption = c("met"),
  ES = c(0.01, 0.03),
  J = c(70),  # J = school
  n_bar = c(30), 
  ICC_k = c(0.15), 
  ICC_jk = c(0.01) 
)

params <- 
  cross_df(design_factors) %>%
  mutate(
    iterations = 5, 
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

# Performance Criteria ----------------------------------------------------
load("EScalc_param.RData")

results <- results %>% 
  left_join(param_list, 
            by = c("method", "cov", "assumption", "ES", 
                   "n_bar", "ICC_k", "ICC_jk")) %>% 
  select(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk, param, param_est,
         everything()) %>% 
  mutate(param = ifelse(!is.na(param_est), param_est, param)) 

calc_performance(results)

# Check "met" condition --------------------------------------------------------
# model parameter fixed
ICC_j = 0.05        
ICC_jk = 0.01        
sparse = .1         
X_m = 0
X_sd = sqrt(50)

dat <- generate_dat(   
  # design factor
  assumption = "met",   # exogeneity or met
  ES = 0.5,
  J = 100,              # school j = {1..J}
  n_bar = 30,          # average number of students per school
  ICC_k = 0.15,         # neighbor ICC
  ICC_jk = 0.01         # neighbor ICC
) 

model <- lmer(y ~ 1 + X_within + cluster_neigh + cluster_school + cluster_cell + (1 | schid) + (1 | neighid)  + (1 | cellid), data = dat)
summary(model) # Variance components and fixed effects should closely match parameters

# neighborhood effects
dat_neigh <- 
  dat %>% 
  select(neighid, c_00k, X_bw_neigh) %>%
  distinct()
dat_neigh %>%
  summarise(
    across(-neighid, list(M = ~ mean(.x), SD = ~ sd(.x)))
  )
sqrt(ICC_k)  # should match sd of c_00k
X_sd         # should match sd of X_bw_neigh

# student effects
sd(dat$e)            
sqrt(1 - ICC_j - ICC_k - ICC_jk) # should match sd of u
lm_fit <- lm(e ~ X, data = dat)
plot(lm_fit)

# Check "endogeneity" condition -------------------------------------------
dat_endo <- generate_dat(
  assumption = "exogeneity",   # exogeneity or met
  ES = 0.5,
  J = 100,              # school j = {1..J}
  n_bar = 30,          # average number of students per school
  ICC_k = 0.15,        # neighbor ICC
  ICC_jk = 0.01        # neighbor ICC
)  

# neighborhood effects
dat_neigh <- 
  dat_endo %>% 
  select(neighid, c_00k, X_bw_neigh) %>%
  distinct()
dat_neigh %>%
  summarise(
    across(-neighid, list(M = ~ mean(.x), SD = ~ sd(.x)))
  )
sqrt(ICC_k)          # should match sd of c_00k
X_sd                 # should match sd of X_bw_neigh

cor(dat_neigh) # should be correlation of 0.4 between X_bw_neigh and c_00k

# school effects
dat_sch <- 
  dat_endo %>% 
  select(schid, b_0j0, X_bw_school) %>%
  distinct()
dat_sch %>%
  summarise(
    across(-schid, list(M = ~ mean(.x), SD = ~ sd(.x)))
  )
sqrt(ICC_j)      # should match sd of c_00k
X_sd             # should match sd of X_bw_schol

cor(dat_sch) # should all be near zero

# student effects
sd(dat_endo$e)
sqrt(1 - ICC_j - ICC_k - ICC_jk) # should match sd of u
lm_fit <- lm(e ~ X, data = dat_endo)
plot(lm_fit)
