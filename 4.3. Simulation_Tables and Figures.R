# Figures
library(tidyverse)
library(sjstats)
library(knitr)
library(kableExtra)
library(simhelpers)
library(patchwork)

rm(list = ls())
load("sim_results_raw.Rdata")

# Cleaning ----------------------------------------------------------------
results <- results %>% 
  mutate(method = recode(method,
                         "uncentered" = "Uncentered",
                         "grand" = "Grand-mean",
                         "cluster_w" = "Within-Cluster",
                         "cluster_bw" = "BW-Cluster",
                         "cell_w" = "Within-Cell",
                         "cell_bw" = "BW-Cell",
                         "fe_crve" = "FE-CRVE",
                         "hybrid_sch_fe" = "Hybrid (School FE/ Neighbor RE)",
                         "hybrid_sch_re" = "Hybrid (School RE/ Neighbor FE)"),
         method = factor(method, 
                         levels = c("Uncentered", "Grand-mean", 
                                    "Within-Cluster", "BW-Cluster", 
                                    "Within-Cell", "BW-Cell",
                                    "FE-CRVE",
                                    "Hybrid (School FE/ Neighbor RE)",
                                    "Hybrid (School RE/ Neighbor FE)")),
         coefficient = recode(ES, "0.1" = "\u03B3 = 0.01",
                              "0.2" = "\u03B3 = 0.02",
                              "0.4" = "\u03B3 = 0.04"),
         coefficient = factor(coefficient, levels = c("\u03B3 = 0.01", 
                                                      "\u03B3 = 0.02", 
                                                      "\u03B3 = 0.04")),
         J = as.factor(J),
         J_f = recode(J, 
                      "20" = "Schools = 20",
                      "70" = "Schools = 70",
                      "150" = "Schools = 150"),
         J_f = factor(J_f, 
                      levels = c("Schools = 20",
                                 "Schools = 70",
                                 "Schools = 150")),
         n_bar = as.factor(n_bar),
         n_bar_f = recode(n_bar, 
                          "30" = "Students/school = 30",
                          "100" = "Students/school = 100"),
         n_bar_f = factor(n_bar_f, 
                          levels = c("Students/school = 30",
                                     "Students/school = 100")),
         ICC_jk = as.factor(ICC_jk)) 

# ANOVA functions ---------------------------------------------------------------
fml <- " ~ as.factor(assumption)+as.factor(method)+as.factor(coefficient)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk)+
        as.factor(assumption)*(as.factor(method)+as.factor(coefficient)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+ 
        as.factor(method)*(as.factor(coefficient)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+
        as.factor(coefficient)*(as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+
        as.factor(J)*(as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+
        as.factor(n_bar)*(as.factor(ICC_k)+as.factor(ICC_jk))+
        as.factor(ICC_k)*as.factor(ICC_jk)"


tidy_func <- function(obj){
  anova_stats(obj) %>% as_tibble %>% 
    select(term, partial.etasq) %>% 
    mutate(size = ifelse(partial.etasq >= 0.14, "(large)", NA),
           size = ifelse(partial.etasq >= 0.06 & partial.etasq < 0.14, "(medium)", size),
           size = ifelse(partial.etasq >= 0.01 & partial.etasq < 0.06, "(small)", size)) %>% 
    mutate_if(is.numeric, round, 3)
}

table_anova <- function(outcome, cov){
  fit1 <- aov(as.formula(paste0(outcome, fml)), data = cov)
  fit1 <- tidy_func(fit1)
  options(knitr.kable.NA = '')
  fit1 <- fit1[1:7,]
  fit1 %>%
    kable(digits = 3) %>%
    kable_styling(full_width = F) 
}

dat_anova <- function(outcome, cov){
  fit1 <- aov(as.formula(paste0(outcome, fml)), data = cov)
  fit1 <- tidy_func(fit1)
  options(knitr.kable.NA = '')
  fit1 <- fit1[1:7,]
}

table_anova_bw <- function(outcome, cov1, cov2){
  fit1 <- aov(as.formula(paste0(outcome, fml)), data = cov1)
  fit1 <- tidy_func(fit1)
  fit1 <- fit1[1:7,]
  
  fit2 <- aov(as.formula(paste0(outcome, fml)), data = cov2)
  fit2 <- tidy_func(fit2)
  fit2 <- fit2[1:7,]
  
  options(knitr.kable.NA = '')
  fit1 %>% left_join(fit2, by = c("term")) %>% 
    kable(digits = 3, 
          col.names = c(" ", "partial.etasq", "size",
                        "partial.etasq", "size")) %>%
    kable_styling(full_width = F) %>% 
    add_header_above(c(" " = 1, "school" = 2, "neighborhood" = 2))
}
# Rate of Convergence -----------------------------------------------------
# descriptive stat
convergence <- results %>% filter(converged == TRUE) %>% 
  group_by(method, cov, assumption, coefficient, J, n_bar, ICC_k, ICC_jk) %>%
  count() %>% 
  summarise(convergence_rate = n/1000, .groups = "drop")

convergence %>% 
  group_by(method) %>% 
  summarise(min = min(convergence_rate)*100,
            Q1 = quantile(convergence_rate, .25)*100,
            mean = mean(convergence_rate)*100, 
            median = median(convergence_rate)*100, 
            Q3 = quantile(convergence_rate, .75)*100,
            max = max(convergence_rate)*100) %>% 
  kable(digits = 2) %>%
  kable_styling(full_width = F)

table_anova("convergence_rate", convergence) # large coefficient: method, ICC_jk

convergence <- results %>% filter(converged == TRUE) %>%
  group_by(method, cov, ICC_jk) %>% 
  count() %>% 
  mutate(n = n/(2 * 3 * 3 * 2* 3), # conditions: assumption * coefficient * J * n_bar * ICC_k
         pct = n/1000) 

fig_converge <- convergence %>%
  ggplot(aes(x = ICC_jk, y = pct, group =  method)) + 
  geom_hline(yintercept = c(1)) +
  geom_line(aes(color = method), size = 1) +
  geom_point(aes(color = method), size = 1.5) +
  # facet_grid(assumption ~ coefficient, scales = "free_y") +
  labs(x = "Cell IUCC", y = "Rate of Convergence") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0)) 
fig_converge

png("Figures/converge.png", units="in", width=8, height=4, res=300)
fig_converge
dev.off()

# Neighborhood Generation Success Rates  ----------------------------------
neigh_n <- results %>% 
  group_by(method, cov, assumption, coefficient, J, n_bar, ICC_k, ICC_jk) %>% 
  summarise(neigh_n = mean(neigh_n),
            neigh_pct = mean(neigh_pct))

table_anova("neigh_pct", neigh_n) # large coefficient: J, n_bar; small coefficient: coefficient

neigh_n <- results %>% 
  group_by(method, cov, assumption, coefficient, J, n_bar_f, ICC_k, ICC_jk) %>% 
  summarise(neigh_n = mean(neigh_n),
            neigh_pct = mean(neigh_pct))

fig_neigh_rate <- neigh_n %>%
  ggplot(aes(x = J, y = neigh_pct)) + 
  geom_hline(yintercept = c(1), linetype=2) +
  geom_boxplot(alpha = .7, lwd = .1, color = "black", fill = "darkgray") + 
  # geom_line(aes(color = method), size = 1) +
  # geom_point(aes(color = method), size = 1.5) +
  facet_wrap(~ n_bar_f) +
  labs(x = "School", y = "Neighborhood Generation Rate") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0)) 
fig_neigh_rate

png("Figures/neigh_rate.png", units="in", width=8, height=4, res=300)
fig_neigh_rate
dev.off()

# Cleaning 2 --------------------------------------------------------------
rm(list = ls())
load("sim_results_1000_performance.Rdata")

results <- results %>% 
  mutate(method_origin = method,
         method = recode(method,
                         "uncentered" = "Uncentered/Grand-mean",
                         "grand" = "Uncentered/Grand-mean",
                         "cluster_w" = "Cluster-mean",
                         "cluster_bw" = "Cluster-mean",
                         "cell_w" = "Cell-mean",
                         "cell_bw" = "Cell-mean",
                         "fe_crve" = "FE-CRVE",
                         "hybrid_sch_fe" = "Hybrid (School FE/ Neighbor RE)",
                         "hybrid_sch_re" = "Hybrid (School RE/ Neighbor FE)"),
         method = factor(method, 
                         levels = c("Uncentered/Grand-mean", 
                                    "Within-cluster mean", 
                                    "Cluster-mean",
                                    "Cell-mean",
                                    "FE-CRVE",
                                    "Hybrid (School FE/ Neighbor RE)",
                                    "Hybrid (School RE/ Neighbor FE)")),
         assumption = recode(assumption,
                             "met" = "Assumptions met",
                             "exogeneity" = "Endogenous N."),
         assumption = factor(assumption, 
                             levels = c("Assumptions met", "Endogenous N.")),
         coefficient = as.factor(ES),
         coefficient_f = recode(ES, "0.1" = "\u03B3 = 0.01",
                                "0.2" = "\u03B3 = 0.02",
                                "0.4" = "\u03B3 = 0.04"),
         coefficient_f = factor(coefficient_f, levels = c("\u03B3 = 0.01", 
                                                          "\u03B3 = 0.02", 
                                                          "\u03B3 = 0.04")),
         J = as.factor(J),
         J_f = recode(J, 
                      "20" = "Schools = 20",
                      "70" = "Schools = 70",
                      "150" = "Schools = 150"),
         J_f = factor(J_f, 
                      levels = c("Schools = 20",
                                 "Schools = 70",
                                 "Schools = 150")),
         n_bar = as.factor(n_bar),
         n_bar_f = recode(n_bar, 
                          "30" = "Students/school = 30",
                          "100" = "Students/school = 100"),
         n_bar_f = factor(n_bar_f, 
                          levels = c("Students/school = 30",
                                     "Students/school = 100")),
         ICC_k = as.factor(ICC_k), 
         ICC_k = recode(ICC_k,
                        "0.05" = "IUCC=0.05",
                        "0.15" = "IUCC=0.15",
                        "0.25" = "IUCC=0.25"),
         ICC_jk = as.factor(ICC_jk), 
         ICC_jk = recode(ICC_jk,
                         "0" = "IUCC (cell) = 0",
                         "0.05" = "IUCC (cell) = 0.05",
                         "0.15" = "IUCC (cell) = 0.15")) 

X <- results %>% filter(cov == "X")       # level-1 cov within
BW_S <- results %>% filter(cov == "BW_S") # level-2 cov between-school
BW_N <- results %>% filter(cov == "BW_N") # level-2 cov between-neighborhood
# Within-Cluster ----------------------------------------------------------

## rel_bias & bias --------------------------------------------------------
aov_rel_bias <- dat_anova("rel_bias", X) 
# assumption, method, coefficient, n_bar (medium)
aov_bias <- dat_anova("bias", X) 
# assumption, method, n_bar, ICC_jk (small - not presented)
options(knitr.kable.NA = '')
aov_rel_bias %>% left_join(aov_bias, by = c("term")) %>% 
  kable(digits = 3, 
        col.names = c(" ", "partial.etasq", "size",
                      "partial.etasq", "size")) %>%
  kable_styling(full_width = F) %>% 
  add_header_above(c(" " = 1, 
                     "Relative Parameter Bias" = 2, 
                     "Absolute Parameter Bias" = 2))

summary(X$rel_bias_mcse) # MCSE
summary(X$bias_mcse)

rel_bias <- X %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.95, 1.05), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias <- X %>%
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "B. Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

fig_bias <- rel_bias / bias + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# tiff("rel_bias.tiff", units="in", width=8, height=6, res=300)
png("Figures/sim1_bias.png", units="in", width=8, height=10, res=300)
fig_bias
dev.off()

## rmse --------------------------------------------------------
# Within-cluster effect
table_anova("rmse", X) # assumption, method, J, n_bar, ICC_k (rather small)

summary(X$rmse_mcse)

rmse <- X %>%
  ggplot(aes(x = n_bar, y = rmse, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ J_f, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 
  # scale_x_discrete(limits=c("20","70","150"))
rmse

png("Figures/sim1_rmse.png", units="in", width=8, height=6, res=300)
rmse
dev.off()


## rel_bias_var --------------------------------------------------------
# Within-cluster effect
table_anova("rel_bias_var", X) # method, J, ICC_jk, n_bar(medium)

summary(X$rel_bias_var_mcse)

rel_bias_var <- X %>%
  ggplot(aes(x = n_bar, y = rel_bias_var, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .3) + 
  facet_grid(ICC_jk ~ J_f, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Relative Bias of SE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))
rel_bias_var

png("Figures/sim1_rel_bias_var.png", units="in", width=8, height=6, res=300)
rel_bias_var
dev.off()


# Between-Cluster ----------------------------------------------------------
BW_S <- BW_S %>% 
  mutate(method = recode(method, 
                         "Cluster-mean" = "CRE",
                         "Cell-mean" = "Correlated-cell RE"))
BW_N <- BW_N %>% 
  mutate(method = recode(method, 
                         "Cluster-mean" = "CRE",
                         "Cell-mean" = "Correlated-cell RE"))
## rel_bias & bias ---------------------------------------------------------

table_anova_bw("rel_bias", BW_S, BW_N) # assumption, J, ICC_k
table_anova_bw("bias", BW_S, BW_N) # assumption, coefficient, J

summary(BW_S$rel_bias_mcse)
summary(BW_N$rel_bias_mcse)

summary(BW_S$bias_mcse)
summary(BW_N$bias_mcse)

### School (Col: J)
rel_bias_bw_s <- BW_S %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) +
  facet_grid(assumption ~ J_f, scales = "free_y") + 
  labs(title = "A. School: Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))

bias_bw_s <- BW_S %>% 
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ J_f, scales = "free_y") +  
  labs(title = "B. School: Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_s / bias_bw_s + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_bias_J_sch.png", units="in", width=8, height=9.5, res=300)
rel_bias_bw_s / bias_bw_s + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

### Neighborhood (Col: J)
rel_bias_bw_n <- BW_N %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ J_f, scales = "free_y") +  
  labs(title = "A. Neighborhood: Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_n <- BW_N %>% 
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ J_f, scales = "free_y") +  
  labs(title = "B. Neighborhood: Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 


rel_bias_bw_n / bias_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_bias_J_ngh.png", units="in", width=8, height=9.5, res=300)
rel_bias_bw_n / bias_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()


### School (Col: IUCC)
rel_bias_bw_s <- BW_S %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "A. School: Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_s <- BW_S %>% 
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "B. School: Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_s / bias_bw_s + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_bias_iucc_sch.png", units="in", width=8, height=9.5, res=300)
rel_bias_bw_s / bias_bw_s + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

### Neighborhood (Col: IUCC)
rel_bias_bw_n <- BW_N %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "A. Neighborhood: Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_n <- BW_N %>% 
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "B. Neighborhood: Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_n / bias_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_bias_iucc_ngh.png", units="in", width=8, height=9.5, res=300)
rel_bias_bw_n / bias_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()


### School (Col: coefficient)
rel_bias_bw_s <- BW_S %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. School: Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_s <- BW_S %>% 
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "B. School: Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_s / bias_bw_s + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_bias_coef_sch.png", units="in", width=8, height=9.5, res=300)
rel_bias_bw_s / bias_bw_s + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

### Neighborhood (Col: coefficient)
rel_bias_bw_n <- BW_N %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. Neighborhood: Relative Parameter Bias",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_n <- BW_N %>% 
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "B. Neighborhood: Absolute Parameter Bias",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_n / bias_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_bias_coef_ngh.png", units="in", width=8, height=9.5, res=300)
rel_bias_bw_n / bias_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

## rmse --------------------------------------------------------
# Within-cluster effect
table_anova_bw("rmse", BW_S, BW_N)  #assumption, ICC_k, n_bar

summary(BW_S$rmse_mcse)
summary(BW_N$rmse_mcse)

### School (Col: coefficient)
rmse_bw_s <- BW_S %>%
  ggplot(aes(x = n_bar, y = rmse, fill = method, color = method)) + 
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. School",
       x = "Number of Students per School", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))


### Neighborhood (Col: ICC_k)
rmse_bw_n <- BW_N %>%
  ggplot(aes(x = n_bar, y = rmse, fill = method, color = method)) + 
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "B. Neighborhood",
       x = "Number of Students per School", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))

rmse_bw_s / rmse_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

png("Figures/sim2_rmse.png", units="in", width=8, height=10, res=300)
rmse_bw_s / rmse_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

### School (Col: J)
rmse_bw_s <- BW_S %>%
  ggplot(aes(x = n_bar, y = rmse, fill = method, color = method)) + 
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ J_f, scales = "free_y") +  
  labs(title = "A. School",
       x = "Number of Students per School", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))

### Neighborhood (Col: J)
rmse_bw_n <- BW_N %>%
  ggplot(aes(x = n_bar, y = rmse, fill = method, color = method)) + 
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ J_f, scales = "free_y") +  
  labs(title = "B. Neighborhood",
       x = "Number of Students per School", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))

rmse_bw_s / rmse_bw_n

png("Figures/sim2_rmse_J.png", units="in", width=8, height=10, res=300)
rmse_bw_s / rmse_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

## rel_bias_var --------------------------------------------------------
table_anova_bw("rel_bias_var", BW_S, BW_N) # coefficient

summary(BW_S$rel_bias_var_mcse)
summary(BW_N$rel_bias_var_mcse)

rel_bias_var_bw_s <- BW_S %>%
  ggplot(aes(x = J, y = rel_bias_var, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .3) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. School",
       x = "Number of Schools", y = "Relative Bias of SE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_x_discrete(limits=c("20","70","150"))

rel_bias_var_bw_n <- BW_N %>%
  ggplot(aes(x = J, y = rel_bias_var, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .3) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "B. Neighborhood",
       x = "Number of Schools", y = "Relative Bias of SE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_x_discrete(limits=c("20","70","150")) 

rel_bias_var_bw_s / rel_bias_var_bw_n

png("Figures/sim2_rel_bias_var.png", units="in", width=8, height=10, res=300)
rel_bias_var_bw_s / rel_bias_var_bw_n + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()

