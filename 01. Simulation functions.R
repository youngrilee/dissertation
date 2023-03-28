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
    mutate(cluster_cell = mean(X)) %>%
    group_by(schid) %>%
    mutate(cluster_school = mean(X)) %>%
    group_by(neighid) %>%
    mutate(cluster_neigh = mean(X)) %>% 
    ungroup() %>%
    mutate(
      grand = mean(X),
      param = gamma_w,
      J = J)
  
  return(dat)
}

# Model-fitting/Estimation-------------------------------------------------
uncentered <- function(dat) {
  
  fit <- lmer(y ~ X + (1|neighid) + (1|schid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      method = "uncentered",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)"))
  
  return(fixed_est)
}

grand <- function(dat) {
  
  fit <- lmer(y ~ I(X - grand) + (1|neighid) + (1|schid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      cov = ifelse(cov == "I(X - grand)", "X", cov), 
      method = "grand",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)")) 
  
  return(fixed_est)
}

cluster_w <- function(dat) {
  
  X_fit <- felm(X ~ 0 | schid + neighid, data = dat)
  dat$X_adapt <- as.numeric(residuals(X_fit))
  fit <- lmer(y ~ X_adapt + (1 | schid) + (1 | neighid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      cov = ifelse(cov == "X_adapt", "X", cov), 
      method = "cluster_w",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)"))
  
  return(fixed_est)
}

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

cell_w <- function(dat) {
  
  X_fit <- felm(X ~ 0 | schid + neighid + cellid, data = dat)
  dat$X_adapt_cell <- as.numeric(residuals(X_fit))
  fit <- lmer(y ~ X_adapt_cell + (1|neighid) + (1|schid) + (1|cellid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      cov = recode(cov, "X_adapt_cell" = "X"), 
      method = "cell_w",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)"))
  
  return(fixed_est)
}

cell_bw <- function(dat) {
  
  X_fit <- felm(X ~ 0 | schid + neighid + cellid, data = dat)
  dat$X_adapt_cell <- as.numeric(residuals(X_fit))
  fit <- lmer(y ~ X_adapt_cell + cluster_neigh + cluster_school + cluster_cell + (1|neighid) + (1|schid) + (1|cellid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      cov = recode(cov, "X_adapt_cell" = "X",
                   "cluster_neigh" = "BW_N",
                   "cluster_school" = "BW_S",
                   "cluster_int" = "BW_C"), 
      method = "cell_bw",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(!cov %in% c("(Intercept)"))
  
  return(fixed_est)
}

fe_crve <- function(dat) {

    fit <- felm(y ~ X | schid + neighid, data = dat)
    
    fixed_est <- coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
      dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
      mutate(method = "fe_crve",
        converged = TRUE
      )
    
    return(fixed_est)
}

    
# hybrid model: sch-FE, nei-RE
hybrid_sch_fe <- function(dat) {
  
  fit <- lmer(y ~ X + factor(schid) + (1 | neighid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      method = "hybrid_sch_fe",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(cov %in% c("X"))
  
  return(fixed_est)
}

# hybrid model: sch-RE, nei-FE
hybrid_sch_re <- function(dat) {
  
  fit <- lmer(y ~ X + factor(neighid) + (1 | schid), data = dat)
  fixed_est <- 
    coef(summary(fit)) %>% as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      method = "hybrid_sch_re",
      converged = (is.na(is.na(fit@optinfo$conv$lme4)[1]))
    ) %>% 
    filter(cov %in% c("X"))
  
  return(fixed_est)
}

# bind_results
estimate <- function(dat) {
  
  CI_level = .95
  crit <- qnorm((1 + CI_level) / 2)
  neigh_n <- length(unique(dat$neighid))
  K <- dat$J[1]*3.5
  
  results <- bind_rows(
    uncentered(dat),
    grand(dat),
    cluster_w(dat),
    cluster_bw(dat),
    cell_w(dat),
    cell_bw(dat),
    fe_crve(dat),
    hybrid_sch_fe(dat),
    hybrid_sch_re(dat)
    ) %>% 
    mutate(
      var = se^2,
      lower_bound = est - crit * se,
      upper_bound = est + crit * se,
      param = dat$param[1],
      neigh_n = length(unique(dat$neighid)),
      neigh_pct = neigh_n/(dat$J[1]*3.5)
    ) %>% 
    as_tibble() 
  return(results)
}


# Performance calculations ------------------------------------------------
calc_performance <- function(results) {
  
  load("EScalc_param.RData")
  
  results <- results %>% 
    left_join(param_list, 
              by = c("method", "cov", "assumption", "ES", 
                     "n_bar", "ICC_k", "ICC_jk")) %>% 
    select(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk, param, param_est,
           everything()) %>% 
    mutate(param = ifelse(!is.na(param_est), param_est, param)) 
  
  abs_crit <- results %>% 
    group_by(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>%
    group_modify(~ calc_absolute(.x, estimates = est, true_param = param))
  
  rel_crit <- results %>%
    group_by(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>%
    group_modify(~ calc_relative(.x, estimates = est, true_param = param)) 
  
  rel_crit_val <- results %>% 
    group_by(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>%
    group_modify(~ calc_relative_var(.x, estimates = est, var_estimates = var))
  
  convergence <- results %>% 
    group_by(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>% 
    summarise(convergence_rate = sum(converged)/n(), .groups = "drop")
  
  
  performance_measures <- 
    abs_crit %>% 
    left_join(rel_crit, 
              by = c("method", "cov", "assumption", "ES", "J", "n_bar",
              "ICC_k", "ICC_jk", "K")) %>% 
    left_join(rel_crit_val, 
              by = c("method", "cov", "assumption", "ES", "J", "n_bar",
                                   "ICC_k", "ICC_jk", "K")) %>% 
    left_join(convergence, by = c("method", "cov", "assumption", "ES", "J", "n_bar",
                                  "ICC_k", "ICC_jk")) %>% 
    rename(N = K)
  
  return(performance_measures)
}

# Simulation driver -------------------------------------------------------
run_sim <- function(iterations,
                    assumption, ES, J, n_bar, ICC_k, ICC_jk, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

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

  # calc_performance(results)
}



