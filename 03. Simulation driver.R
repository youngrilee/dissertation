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

# simulation parameters as vectors/lists
design_factors <- list(
  assumption = c("met", "exogeneity"),
  ES = c(0.01, 0.03, 0.05),
  J = c(20, 70, 150),  # J = school
  n_bar = c(30, 100), 
  ICC_k = c(0.05, 0.15, 0.25), 
  ICC_jk = c(0, 0.05, 0.15) 
)

params <- 
  cross_df(design_factors) %>%
  mutate(
    iterations = 1000, 
    seed = 20230129 + 1:n()
  )
  
nrow(params)
head(params)

#--------------------------------------------------------
# run simulations in parallel - future + furrr workflow
#--------------------------------------------------------
# pmap
options(error=recover)
plan("multisession") 

system.time(
  results <-
    params %>%
    mutate(res = future_pmap(., .f = run_sim, .options = furrr_options(seed=NULL))) %>% 
    unnest(cols = res)
) 


# Check for convergence
results %>% filter(converged == TRUE) %>% 
  group_by(method, cov) %>% count() # 216 * N_iteration

load("EScalc_param.RData")

results_raw <- results
results <- calc_performance(results)

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------
session_info <- sessionInfo() 
run_date <- date() 
save(results, params, session_info, run_date, file = "sim_results.Rdata")
