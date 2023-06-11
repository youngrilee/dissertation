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
library(patchwork)
rm(list = ls())

# Basic condition 
assumption = "met"   # exogeneity or met
ES = 0.5
J = 20               # school j = {1..J}
n_bar = 10           # average number of students per school
ICC_k = 0.05         # neighbor ICC
ICC_jk = 0.01    

# Data Generating Model----------------------------------------------------
generate_dat <- function(assumption, ES, J, n_bar, ICC_k, ICC_jk){
  
  # model parameter fixed
  ICC_j = 0.05        
  sparse = .1         
  K = J * 3.5

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
  
  # calc cluster means
  dat <- dat %>% group_by(schid, neighid) %>% 
    summarise(n = n()) %>% ungroup() %>% 
    summarise(min = min(n),
              q1 = quantile(n, .25),
              q2 = quantile(n, .5),
              q3 = quantile(n, .75),
              max = max(n),
              median = median(n),
              mean = mean(n),
              sd = sd(n))
    
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

results <-
  rerun(5, {
    dat <- generate_dat(   
      assumption = "met",   # exogeneity or met
      ES = 0.5,
      J = 200,              # school j = {1..J}
      n_bar = 30,          # average number of students per school
      ICC_k = 0.05,         # neighbor ICC
      ICC_jk = 0.01         # neighbor ICC
    ) 
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
      ) 
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

save(results_final, file = "cell_n_summary.RData")

# graph -------------------------------------------------------------------
load("cell_n_summary.Rdata")
temp <- results_final %>% 
  group_by(J, n_bar) %>% 
  summarise(min = min(min),
            q1 = mean(q1),
            q2 = mean(q2),
            q3 = mean(q3),
            max = max(max),
            median = mean(median),
            mean = mean(mean),
            sd = mean(sd)) %>% 
  ungroup() %>% 
  mutate(J = as.factor(J))

temp
  
n_bar_30 <- temp %>% filter(n_bar == 30) %>% 
  ggplot(aes(x = J, group = J,
             ymin = min, 
             lower = q1, 
             middle = q2,
             upper = q3,
             ymax = max)) +
  geom_boxplot(stat = "identity", fill = "dark gray") +
  labs(title = "Number of Students per School = 30",
       x = "Schools") +
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13))

n_bar_100 <- temp %>% filter(n_bar == 100) %>% 
  ggplot(aes(x = J, group = J,
             ymin = min, 
             lower = q1, 
             middle = q2,
             upper = q3,
             ymax = max)) +
  geom_boxplot(stat = "identity", fill = "dark gray") +
  labs(title = "Number of Students per School = 100",
       x = "Schools") +
  theme_bw() + 
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13))


n_bar_30 + n_bar_100

png("Figures/cell_n_distribution.png", units="in", width=8, height=4, res=300)
n_bar_30 + n_bar_100
dev.off()