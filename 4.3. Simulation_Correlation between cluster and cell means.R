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
    ungroup() 
  
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


data.frame(cor_s_c = cor(dat$cluster_school, dat$cell_mean),
           cor_n_c = cor(dat$cluster_neigh, dat$cell_mean))

# calc cor -----------------------------------------------------------------
# bind_results
estimate <- function(dat) {
  
  results <- data.frame(cor_s_c = cor(dat$cluster_school, dat$cell_mean),
                    cor_n_c = cor(dat$cluster_neigh, dat$cell_mean))
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

# pmap
options(error=recover)
plan("multisession") 

system.time(
  results_final <-
    params %>%
    mutate(res = future_pmap(., .f = run_sim, .options = furrr_options(seed=NULL))) %>% 
    unnest(cols = res)
)

save(results_final, file = "cor_multicollinearity.RData")






# descriptive stat --------------------------------------------------------
rm(list = ls())
load("cor_multicollinearity.RData")

summary(results_final$cor_s_c)
summary(results_final$cor_n_c)


# anova -------------------------------------------------------------------
library(sjstats)
fit1 <- aov(cor_s_c ~ as.factor(assumption)+as.factor(ES)+as.factor(J) + as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk), 
            data = results_final)
anova_stats(fit1) %>% as_tibble %>% 
  select(term, partial.etasq) %>% 
  mutate(size = ifelse(partial.etasq >= 0.14, "(large)", NA),
         size = ifelse(partial.etasq >= 0.06 & partial.etasq < 0.14, "(medium)", size),
         size = ifelse(partial.etasq >= 0.01 & partial.etasq < 0.06, "(small)", size)) %>% 
  mutate_if(is.numeric, round, 3) 

fit1 <- aov(cor_n_c ~ as.factor(assumption)+as.factor(ES)+as.factor(J) + as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk), 
            data = results_final)
anova_stats(fit1) %>% as_tibble %>% 
  select(term, partial.etasq) %>% 
  mutate(size = ifelse(partial.etasq >= 0.14, "(large)", NA),
         size = ifelse(partial.etasq >= 0.06 & partial.etasq < 0.14, "(medium)", size),
         size = ifelse(partial.etasq >= 0.01 & partial.etasq < 0.06, "(small)", size)) %>% 
  mutate_if(is.numeric, round, 3) 



# graph -------------------------------------------------------------------
cor <- results_final %>% 
  select(-c(iterations, seed)) %>% 
  pivot_longer(7:8, names_to = "cor_type", values_to = "cor") %>% 
  mutate(J = as.factor(J),
         cor_type = ifelse(cor_type == "cor_s_c", "School", "Neighborhood"),
         cor_type = factor(cor_type, levels = c("School", "Neighborhood"))) %>% 
  ggplot(aes(x = J, y = cor, fill = cor_type, color = cor_type)) + 
  geom_boxplot(alpha = .7, lwd = .1) +
  labs(x = "The Number of Schools", y = "Correlation") + 
  theme_bw() +
  # scale_y_continuous(limits = c(0, 1)) +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0)) 

png("Figures/multicollinearity.png", units="in", width=9, height=4, res=300)
cor
dev.off()