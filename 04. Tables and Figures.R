library(tidyverse)
library(sjstats)
library(knitr)
library(kableExtra)
library(simhelpers)
library(patchwork)
library(wesanderson)
rm(list = ls())

# Cleaning: first 1000 rep ------------------------------------------------
load("sim_results_raw.Rdata")

results <- results %>% 
  mutate(method = recode(method,
                         "uncentered" = "Uncentered",
                         "grand" = "Grand-mean",
                         "cluster_w" = "Within-Cluster",
                         "cluster_bw" = "CRE-Cluster",
                         "cell_w" = "Within-Cell",
                         "cell_bw" = "CRE-Cell",
                         "fe_crve" = "FE-CRVE",
                         "hybrid_sch_fe" = "Hybrid (School FE)",
                         "hybrid_sch_re" = "Hybrid (School RE)"),
         method = factor(method, 
                         levels = c("Uncentered", "Grand-mean", 
                                    "Within-Cluster", "CRE-Cluster", 
                                    "Within-Cell", "CRE-Cell",
                                    "FE-CRVE",
                                    "Hybrid (School FE)",
                                    "Hybrid (School RE)")),
         J = as.factor(J),
         n_bar_f = recode(n_bar, 
                          "30" = "Students per school = 30",
                          "100" = "Students per school = 100"),
         n_bar_f = factor(n_bar_f, levels = c("Students per school = 30",
                                              "Students per school = 100")),
         n_bar_f = as.factor(n_bar_f),
         ICC_jk = as.factor(ICC_jk)) 

# ANOVA functions ---------------------------------------------------------------
fml <- " ~ as.factor(method)+as.factor(assumption)+as.factor(coefficient)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk)+as.factor(method)*(as.factor(assumption)+as.factor(coefficient)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+as.factor(assumption)*(as.factor(coefficient)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+as.factor(coefficient)*(as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+as.factor(J)*(as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+as.factor(n_bar)*(as.factor(ICC_k)+as.factor(ICC_jk))+as.factor(ICC_k)*as.factor(ICC_jk)"


tidy_func <- function(obj){
  anova_stats(obj) %>% as_tibble %>% 
    select(term, partial.etasq) %>% 
    mutate(size = ifelse(partial.etasq >= 0.14, "(large)", NA),
           size = ifelse(partial.etasq >= 0.06 & partial.etasq < 0.14, "(medium)", size),
           size = ifelse(partial.etasq >= 0.01 & partial.etasq < 0.06, "(small)", size)) %>% 
    mutate_if(is.numeric, round, 3)
}

table_anova <- function(outcome, cov){
  cov <- cov %>% rename(coefficient = ES)
  fit1 <- aov(as.formula(paste0(outcome, fml)), data = cov)
  fit1 <- tidy_func(fit1)
  options(knitr.kable.NA = '')
  fit1 <- fit1[1:7,]
  fit1 %>%
    kable(digits = 3) %>%
    kable_styling(full_width = F) 
}

table_anova_bw <- function(outcome, cov1, cov2){
  cov1 <- cov1 %>% rename(coefficient = ES)
  cov2 <- cov2 %>% rename(coefficient = ES)
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

## Rate of Convergence ----------------------------------------------------
convergence <- results %>% filter(converged == TRUE) %>% 
  group_by(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>% count() %>% 
  summarise(convergence_rate = n/1000, .groups = "drop")

table_anova("convergence_rate", convergence) # large ES: method, ICC_jk

convergence <- results %>% filter(converged == TRUE) %>% 
  group_by(method, cov, ICC_jk) %>% 
  count() %>% ungroup() %>%  
  mutate(n = n/(2 * 3 * 3 * 2 * 3), # assumption * coefficient * J * n_bar *ICC_k
         pct = n/1000) 
convergence

convergence %>% filter(method %in% c("Within-Cell", "CRE-Cell")) %>% 
  summarise(max = max(pct),
            min = min(pct))

# boxplot
fig_converge1 <-
  convergence %>%
  ggplot(aes(x = method, y = pct)) +
  geom_boxplot(fill = "darkgray", color = "black") +
  geom_hline(yintercept = c(1)) +
  labs(x = "", y = "Rate of Convergence") + 
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.title.y = element_text(vjust = +3), # make it left
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text()) +
  scale_y_continuous(limits = c(0.4, 1)) +
  guides(colour = guide_legend(nrow = 3))
  
fig_converge1
png("Figures/converge1.png", units="in", width=5, height=3, res=300)
fig_converge1
dev.off()

# line plot
fig_converge2 <-
  convergence %>%
  ggplot(aes(x = ICC_jk, y = pct, group = method)) +
  geom_line(aes(color = method), size = 1) +
  geom_point(aes(color = method), size = 1.5) +
  geom_hline(yintercept = c(1)) +
  labs(x = "Cell IUCC", y = "Rate of Convergence") + 
  theme_minimal() +
  theme(text = element_text(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.title.y = element_text(vjust = +3), # make it left
        axis.title.x = element_text(vjust = -0.75)) + # make it lower
  scale_y_continuous(limits = c(0.4, 1)) +
  guides(colour = guide_legend(nrow = 3))
fig_converge2

png("Figures/converge2.png", units="in", width=6, height=4, res=300)
fig_converge2
dev.off()


## Neighborhood Generation Success Rates  ----------------------------------
neigh_n <- results %>% 
  group_by(method, cov, assumption, ES, J, n_bar, ICC_k, ICC_jk) %>% 
  summarise(neigh_n = mean(neigh_n),
            neigh_pct = mean(neigh_pct))

table_anova("neigh_pct", neigh_n) # large ES: J, n_bar; small ES: ES

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

# Cleaning 2: perfect 1000 rep-------------------------------------------
load("sim_results_1000_performance.Rdata")

results <- results %>% 
  mutate(method = recode(method,
                         "uncentered" = "Uncentered/Grand-mean",
                         "grand" = "Uncentered/Grand-mean",
                         "cluster_w" = "Cluster-mean",
                         "cluster_bw" = "CRE-Cluster",
                         "cell_w" = "Within-Cell",
                         "cell_bw" = "CRE-Cell",
                         "fe_crve" = "FE-CRVE",
                         "hybrid_sch_fe" = "Hybrid (School FE)",
                         "hybrid_sch_re" = "Hybrid (School RE)"),
         method = factor(method, 
                         levels = c("Uncentered/Grand-mean", 
                                    "Within-cluster mean", 
                                    "Cluster-mean",
                                    "Cell-mean",
                                    "FE-CRVE",
                                    "Hybrid (School FE)",
                                    "Hybrid (School RE)")),
         assumption = recode(assumption,
                             "met" = "Assumptions met",
                             "exogeneity" = "Exogeneity"),
         assumption = factor(assumption, 
                             levels = c("Assumptions met", "Exogeneity")),
         coefficient_f = recode(ES, "0.1" = "\u03B3 = 0.01",
                     "0.2" = "\u03B3 = 0.02",
                     "0.4" = "\u03B3 = 0.04"),
         coefficient_f = factor(coefficient_f, levels = c("\u03B3 = 0.01", 
                                    "\u03B3 = 0.02", 
                                    "\u03B3 = 0.04")),
         J = as.factor(J),
         n_bar_f = recode(n_bar, 
                          "30" = "Students/school = 30",
                          "100" = "Students/school = 100"),
         n_bar_f = factor(n_bar_f, 
                          levels = c("Students/school = 30",
                                     "Students/school = 100")),
         n_bar = as.factor(n_bar),
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

## rel_bias ---------------------------------------------------------------
table_anova("rel_bias", X) # assumption, method, coefficients, n_bar (medium)

summary(X$rel_bias_mcse) # MCSE

rel_bias <- X %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.95, 1.05), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 
rel_bias

# tiff("rel_bias.tiff", units="in", width=8, height=6, res=300)
png("Figures/sim1_1_rel_bias.png", units="in", width=8, height=6, res=300)
rel_bias
dev.off()

## bias --------------------------------------------------------
# Within-cluster effect
table_anova("bias", X) # assumption, method, n_bar, ICC_jk (small - not presented)

summary(X$bias_mcse)

bias <- X %>%
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_wrap(~ assumption) +  
  labs(x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 
bias

png("Figures/sim1_2_bias.png", units="in", width=8, height=5, res=300)
bias
dev.off()

## rel_bias_var --------------------------------------------------------
# Within-cluster effect
table_anova("rel_bias_var", X) # ICC_jk, J, method (medium), n_bar(medium)

summary(X$rel_bias_var_mcse)

rel_bias_var <- X %>%
  ggplot(aes(x = J, y = rel_bias_var, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .3) + 
  facet_grid(n_bar_f ~ ICC_jk, scales = "free_y") +  
  labs(x = "Number of Schools", y = "Relative Bias of SE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_x_discrete(limits=c("20","70","150")) 
rel_bias_var

png("Figures/sim1_3_rel_bias_var.png", units="in", width=8, height=6, res=300)
rel_bias_var
dev.off()

## rmse --------------------------------------------------------
# Within-cluster effect
table_anova("rmse", X) # assumption, method, J, n_bar, ICC_k

summary(X$rmse_mcse)

rmse <- X %>%
  ggplot(aes(x = J, y = rmse, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ n_bar_f, scales = "free_y") +  
  labs(x = "Number of Schools", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_x_discrete(limits=c("20","70","150"))
rmse

png("Figures/sim1_4_rmse.png", units="in", width=8, height=6, res=300)
rmse
dev.off()


# Between-Cluster ----------------------------------------------------------


## rel_bias ----------------------------------------------------------------
table_anova_bw("rel_bias", BW_S, BW_N) # assumption, J, ICC_k

summary(BW_S$rel_bias_mcse)
summary(BW_N$rel_bias_mcse)

rel_bias_bw_s <- BW_S %>%
  ggplot(aes(x = J, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "A. School",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_n <- BW_N %>%
  ggplot(aes(x = J, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(title = "B. Neighborhood",
       x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

rel_bias_bw_s + rel_bias_bw_n

png("Figures/sim2_1_rel_bias_bw.png", units="in", width=8, height=6, res=300)
rel_bias_bw_s + rel_bias_bw_n
dev.off()

## bias --------------------------------------------------------------------
table_anova_bw("bias", BW_S, BW_N) # assumption, coefficients, J

summary(BW_S$bias_mcse)
summary(BW_N$bias_mcse)

bias_bw_s <- BW_S %>% 
  ggplot(aes(x = J, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. School",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_n <- BW_N %>% 
  ggplot(aes(x = J, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "B. Neighborhood",
       x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) 

bias_bw_s + bias_bw_n

png("Figures/sim2_2_bias_bw.png", units="in", width=8, height=6, res=300)
bias_bw_s + bias_bw_n
dev.off()

## rel_bias_var --------------------------------------------------------
table_anova_bw("rel_bias_var", BW_S, BW_N) # coefficients

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

rel_bias_var_bw_s + rel_bias_var_bw_n

png("Figures/sim2_3_rel_bias_var_bw.png", units="in", width=8, height=6, res=300)
rel_bias_var_bw_s + rel_bias_var_bw_n
dev.off()

## rmse --------------------------------------------------------
# Within-cluster effect
table_anova_bw("rmse", BW_S, BW_N)  #assumption, ICC_k, n_bar

summary(BW_S$rmse_mcse)
summary(BW_N$rmse_mcse)

rmse_bw_s <- BW_S %>%
  ggplot(aes(x = J, y = rmse, fill = method, color = method)) + 
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ coefficient_f, scales = "free_y") +  
  labs(title = "A. School",
       x = "Number of Schools", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        plot.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1))

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

rmse_bw_s + rmse_bw_n

png("Figures/sim2_4_rmse_bw.png", units="in", width=8, height=6, res=300)
rmse_bw_s + rmse_bw_n
dev.off()
