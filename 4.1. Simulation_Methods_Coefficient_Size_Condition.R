# Empirical Correlations between the Covariates and the Corresponding Errors
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
  X_sd = sqrt(25)
  K = J * 3.5
  
  # set sigma and random effects based on ICC
  tau_J00 = ICC_j  # School
  tau_K00 = ICC_k  # Neighborhood
  tau_JK0 = ICC_jk # Interaction
  sigma = sqrt(1-tau_J00-tau_K00-tau_JK0) 
  
  # coefficients
  gamma_w = gamma_b_j = gamma_b_k = gamma_b_jk = ES*(1/10)
  
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
  # X_bw_neigh <- 
  #   dat %>% 
  #   group_by(neighid) %>% 
  #   summarise() %>% 
  #   mutate(X_bw_neigh = rnorm(nrow(.), mean = X_m, sd = X_sd))
  
  if (assumption == "exogeneity") {
    
    # neighborhood random effect when the exogeneity assumption is violated
    r <- 0.4 # correlation between X_bw_neigh and b_0j0
    
    neighbordata <-
      # X_bw_neigh %>% 
      dat %>%
      group_by(neighid) %>%
      summarise() %>%
      mutate(
        X_bw_neigh = rnorm(nrow(.), mean = X_m, sd = X_sd),
        v_00k = rnorm(nrow(.), mean = 0, sd = sqrt((1 - r^2) * tau_K00)),
        c_00k = r * sqrt(tau_K00 / (X_sd ^ 2)) * X_bw_neigh + v_00k
      ) %>% 
      dplyr::select(-c(v_00k))
    
  } else {
    
    # neighborhood random effect when all assumptions are met:
    neighbordata <-
      # X_bw_neigh %>%
      dat %>%
      group_by(neighid) %>%
      summarise() %>%
      mutate(
        c_00k = rnorm(nrow(.), mean = 0, sd = sqrt(tau_K00)),  
        X_bw_neigh = rnorm(nrow(.), mean = X_m, sd = X_sd)
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
      gamma = ES*(1/10),
      ES_j = gamma_b_j*X_sd/sqrt(tau_J00),
      ES_k = gamma_b_k*X_sd/sqrt(tau_K00),
      ES_jk = gamma_b_jk*X_sd/sqrt(tau_JK0)
    ) 
  
  
  
  return(dat)
}

# test
dat <- generate_dat(   
  # design factor
  assumption = "met",   # exogeneity or met
  ES = 0.5,
  J = 20,              # school j = {1..J}
  n_bar = 10,          # average number of students per school
  ICC_k = 0.05,         # neighbor ICC
  ICC_jk = 0.01         # neighbor ICC
)


# calc ES -----------------------------------------------------------------
cor_school <- function(dat) {
  
  # coefficients
  cor_school <-
    dat %>%
    group_by(schid) %>%
    summarise(X_bw_school = mean(X_bw_school),
              y_bw_school = mean(y))
  cor_school <- cor(cor_school$X_bw_school, cor_school$y_bw_school)
  return(cor_school)
}


cor_neigh <- function(dat) {
  
  cor_neigh <-
    dat %>%
    group_by(neighid) %>%
    summarise(X_bw_neigh = mean(X_bw_neigh),
              y_bw_neigh = mean(y))
  cor_neigh <- cor(cor_neigh$X_bw_neigh, cor_neigh$y_bw_neigh)
  return(cor_neigh)
}

cor_int <- function(dat) {
  
  cor_int <-
    dat %>%
    group_by(cellid) %>%
    summarise(X_bw_int = mean(X_bw_int),
              y_bw_int = mean(y))
  cor_int <- cor(cor_int$X_bw_int, cor_int$y_bw_int)
  return(cor_int)
}

# bind_results
estimate <- function(dat) {
  
  results <- dat %>% 
    mutate(ES_j_est = gamma * sd(X_bw_school) / sd(b_0j0),
           ES_k_est = gamma * sd(X_bw_neigh) / sd(c_00k),
           ES_jk_est = gamma * sd(X_bw_int) / sd(d_0jk),
           cor_j = cor_school(dat),
           cor_k = cor_neigh(dat),
           cor_jk = cor_int(dat)) %>% 
    summarise(ES_j = mean(ES_j),
              ES_k = mean(ES_k),
              ES_jk = mean(ES_jk),
              ES_j_est = mean(ES_j_est),
              ES_k_est = mean(ES_k_est),
              ES_jk_est = mean(ES_jk_est),
              cor_j = mean(cor_j),
              cor_k = mean(cor_k),
              cor_jk = mean(cor_jk)) %>% 
  as_tibble() 
  
  return(results)
}

# test
estimate(dat = dat) 

results <-
  rerun(5, {
    dat <- generate_dat(   
      assumption = "met",   # exogeneity or met
      ES = 0.5,
      J = 200,              # school j = {1..J}
      n_bar = 30,          # average number of students per school
      ICC_k = 0.05,         # neighbor ICC
      ICC_jk = 0.01         # neighbor ICC
    ) %>%
    estimate()
  }) %>%
  bind_rows() 

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
          J = 200,             # school j = {1..J}
          n_bar = 30,           # average number of students per school
          ICC_k = 0.15,         # neighbor ICC
          ICC_jk = 0.01,        # neighbor ICC
          seed = NULL)
results

# simulation parameters as vectors/lists
design_factors <- list(
  assumption = c("met", "exogeneity"),
  ES = c(0.1, 0.2, 0.4),
  J = c(200),  # J = school
  n_bar = c(30, 100), 
  ICC_k = c(0.05, 0.15, 0.25), 
  ICC_jk = c(0, 0.05, 0.15) 
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
  results_final <-
    params %>%
    mutate(res = future_pmap(., .f = run_sim, .options = furrr_options(seed=NULL))) %>% 
    unnest(cols = res)
)

save(results_final, file = "calc_cor_raw.RData")

# calc_cor ------------------------------------------------------------------
calc_cor <- results_final %>% 
  group_by(assumption, ES, n_bar, ICC_k, ICC_jk) %>% 
  summarise(ES_j = mean(ES_j),
            ES_k = mean(ES_k),
            ES_jk = mean(ES_jk),
            ES_j_est = mean(ES_j_est),
            ES_k_est = mean(ES_k_est),
            ES_jk_est = mean(ES_jk_est),
            cor_j = mean(cor_j),
            cor_k = mean(cor_k),
            cor_jk = mean(cor_jk)) %>% ungroup() %>% 
  mutate(gamma = ES * 0.1) %>% 
  select(assumption, ES, gamma, everything()) %>% 
  arrange(desc(assumption))
  
write.csv(calc_cor, "calc_cor.csv", row.names = FALSE)

# graph -------------------------------------------------------------------
calc_cor <- read.csv("calc_cor.csv")

es <- calc_cor %>% 
  select(assumption, ES, gamma, n_bar, ICC_k, ICC_jk, ES_j, ES_k, ES_jk) %>% 
  pivot_longer(ES_j:ES_jk, names_to = "dimension", values_to = "es")

es_est <- calc_cor %>% 
  select(assumption, ES, gamma, n_bar, ICC_k, ICC_jk, 
         ES_j_est, ES_k_est, ES_jk_est) %>% 
  rename(ES_j = ES_j_est,
         ES_k = ES_k_est,
         ES_jk = ES_jk_est) %>% 
  pivot_longer(ES_j:ES_jk, names_to = "dimension", values_to = "es_est")

dat_fig <- es %>% 
  left_join(es_est, 
            by = c("assumption", "ES", "gamma", "n_bar", "ICC_k", "ICC_jk",
                   "dimension")) %>% 
  pivot_longer(es:es_est, names_to = "group", values_to = "est") %>% 
  mutate(gamma = recode(gamma,
                        "0.01" = "0.01",
                        "0.02" = "0.02",
                        "0.04" = "0.04"),
         gamma = factor(gamma, levels = c("0.01", "0.02", "0.04")),
         dimension = recode(dimension,
                            "ES_j" = "School",
                            "ES_k"= "Neighborhood",
                            "ES_jk" = "Cell-Interaction"),
         dimension = factor(dimension, 
                            levels = c("School", 
                                       "Neighborhood", 
                                       "Cell-Interaction")),
         group = recode(group, 
                        "es" = "Parameter",
                        "es_est" = "Estimate"),
         group = factor(group, levels = c("Parameter", "Estimate")),
         ICC_k = as.factor(ICC_k), 
         ICC_k = recode(ICC_k,
                        "0.05" = "ICC(N) = 0.05",
                        "0.15" = "ICC(N) = 0.15",
                        "0.25" = "ICC(N) = 0.25"),
         ICC_jk = as.factor(ICC_jk), 
         ICC_jk = recode(ICC_jk,
                         "0" = "ICC(C) = 0",
                         "0.05" = "ICC(C) = 0.05",
                         "0.15" = "ICC(C) = 0.15")) %>% 
  filter_all(all_vars(!is.infinite(.)))

# anova
library(sjstats)
fit1 <- aov(est ~ as.factor(assumption)+as.factor(gamma)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk), data = dat_fig)
anova_stats(fit1) %>% as_tibble %>% 
  select(term, partial.etasq) %>% 
  mutate(size = ifelse(partial.etasq >= 0.14, "(large)", NA),
         size = ifelse(partial.etasq >= 0.06 & partial.etasq < 0.14, "(medium)", size),
         size = ifelse(partial.etasq >= 0.01 & partial.etasq < 0.06, "(small)", size)) %>% 
  mutate_if(is.numeric, round, 3)

# using bar
bar <- dat_fig %>% 
  ggplot(aes(x = gamma, y = est, fill = group)) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  facet_grid(dimension ~ ICC_k + ICC_jk, scales = "free_y") +
  labs(x = paste0("Coefficient"), y = "Correlation") +
  scale_fill_manual(values=c('darkgray','lightgray'))+
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))

png("Figures/Correlation.png", units="in", width=9, height=5.25, res=300)
bar 
dev.off()

# using line
# dat_fig %>% 
#   group_by(gamma, ICC_k, ICC_jk, dimension, group) %>%
#   summarise(est = mean(est)) %>% ungroup() %>%
#   ggplot(aes(x = gamma, y = est, group = group)) + 
#   geom_line(aes(color = group), size = 1.0) +
#   geom_point(aes(color = group), size = 1.5) +
#   facet_grid(dimension ~ ICC_k + ICC_jk, scales = "free_y") +
#   labs(x = paste0("Coefficient Condition"), y = "Correlation") +
#   theme_bw() +
#   theme(text = element_text(size = 13),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         plot.caption=element_text(hjust = 0),
#         axis.text.x = element_text(angle=90, hjust=1))