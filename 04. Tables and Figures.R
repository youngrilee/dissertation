# Figures
library(tidyverse)
library(sjstats)
library(knitr)
library(kableExtra)
library(simhelpers)

rm(list = ls())
load("sim_results_first1000.RData")
load("sim_results_first1000_performance.RData")
names(results)
head(results)
table(results$method)
results %>% filter(method %in% c("cell_bw", "cluster_bw")) %>% 
  filter(cov == "BW_N") %>% View()
table(results$assumption)
results <- 
  results %>% 
  mutate(method = recode(method,
                         "uncentered" = "Uncentered",
                         "grand" = "Grand-mean",
                         "cluster_w" = "Within-cluster mean",
                         "cluster_bw" = "Between-cluster mean",
                         "cell_w" = "Within-cell mean",
                         "cell_bw" = "Between-cell mean"),
         method = factor(method, 
                         levels = c("Uncentered", "Grand-mean", 
                                    "Within-cluster mean", 
                                    "Between-cluster mean",
                                    "Within-cell mean", 
                                    "Between-cell mean")),
         ES = recode(ES, "0.01" = "r = 0.1",
                     "0.03" = "r = 0.3",
                     "0.05" = "r = 0.5"),
         ES = factor(ES, levels = c("r = 0.1", "r = 0.3", "r = 0.5")),
         J = as.factor(J),
         n_bar_f = recode(n_bar, 
                          "30" = "Students per school = 30",
                          "100" = "Students per school = 100"),
         n_bar_f = factor(n_bar_f, levels = c("Students per school = 30",
                                              "Students per school = 100")),
         n_bar = as.factor(n_bar),
         ICC_k = as.factor(ICC_k), 
         ICC_k = recode(ICC_k,
                          "0.05" = "IUCC (neighbor) = 0.05",
                          "0.15" = "IUCC (neighbor) = 0.15",
                          "0.25" = "IUCC (neighbor) = 0.25"),
         ICC_jk = as.factor(ICC_jk), 
         ICC_jk = recode(ICC_jk,
                        "0.01" = "IUCC (cell) = 0.01",
                        "0.05" = "IUCC (cell) = 0.05"),
         assumption = recode(assumption,
                             "met" = "Assumptions met",
                             "exogeneity" = "Exogeneity"),
         assumption = factor(assumption, 
                             levels = c("Assumptions met", "Exogeneity"))) 
table(results$cov)
X <- results %>% filter(cov == "X")       # level-1 cov within
BW_N <- results %>% filter(cov == "BW_N") # level-2 cov between-neighborhood
BW_S <- results %>% filter(cov == "BW_S") # level-2 cov between-school

# ANOVA functions ---------------------------------------------------------------
fml <- " ~ as.factor(assumption)+as.factor(method)+as.factor(ES)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk)+
        as.factor(assumption)*(as.factor(method)+as.factor(ES)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+ 
        as.factor(method)*(as.factor(ES)+as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+
        as.factor(ES)*(as.factor(J)+as.factor(n_bar)+as.factor(ICC_k)+as.factor(ICC_jk))+
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
  fit1 %>%
    kable(digits = 3) %>%
    kable_styling(full_width = F) 
}

head(results)
table(results$method, results$cov)

# WITHIN ----------------------------------------------------------

# rel_bias --------------------------------------------------------
# Within-cluster effect
table_anova("rel_bias", X) # assumption, method, ES, n_bar (medium)

summary(X$rel_bias_mcse)

rel_bias <- X %>%
  ggplot(aes(x = n_bar, y = rel_bias, fill = method, color = method)) + 
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
        axis.text.x = element_text(angle=90, hjust=1)) +
  # scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))
    
rel_bias

tiff("rel_bias.tiff", units="in", width=8, height=6, res=300)
png("rel_bias.png", units="in", width=8, height=6, res=300)
rel_bias
dev.off()

# bias --------------------------------------------------------
# Within-cluster effect
table_anova("bias", X) # assumption, method, n_bar, ICC_k (small)

summary(X$bias_mcse)

bias <- X %>%
  ggplot(aes(x = n_bar, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ICC_k, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  # scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))

bias

png("bias.png", units="in", width=8, height=6, res=300)
bias
dev.off()

# rel_bias_var --------------------------------------------------------
# Within-cluster effect
table_anova("rel_bias_var", X) # method, J, ICC_jk, n_bar(medium)

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
  scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))

rel_bias_var

png("rel_bias_var.png", units="in", width=8, height=6, res=300)
rel_bias_var
dev.off()

# rmse --------------------------------------------------------
# Within-cluster effect
table_anova("rmse", X) # assumption, method, J, n_bar, ICC_k (rather small)

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
  scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))

rmse

png("rmse.png", units="in", width=8, height=6, res=300)
rmse
dev.off()

converge <- X %>%
  ggplot(aes(x = method, y = convergence_rate, fill = method, color = method)) + 
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) +
  labs(x = " ", y = "Convergence Rate") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x=element_blank()) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))

converge

png("converge.png", units="in", width=8, height=4, res=300)
converge
dev.off()

# BETWEEN -----------------------------------------------------------------

# rel_bias ----------------------------------------------------------------
table_anova("rel_bias", BW_N) # assumption, method, ES, n_bar (medium)

rel_bias <- BW_N %>%
  ggplot(aes(x = J, y = rel_bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ES, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Relative Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  # scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#00BFC4", "#F564E3")) +
  scale_color_manual(values = c("#00BFC4", "#F564E3"))

rel_bias

# tiff("rel_bias.tiff", units="in", width=8, height=6, res=300)
png("rel_bias_bw.png", units="in", width=8, height=6, res=300)
rel_bias
dev.off()


# bias --------------------------------------------------------------------
table_anova("bias", BW_N) # assumption, method, n_bar, ICC_k (small)

bias <- BW_N %>%
  ggplot(aes(x = J, y = bias, fill = method, color = method)) + 
  geom_hline(yintercept = c(0)) +
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ES, scales = "free_y") +  
  labs(x = "Number of Students per School", y = "Parameter Bias") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  # scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#00BFC4", "#F564E3")) +
  scale_color_manual(values = c("#00BFC4", "#F564E3"))

bias

png("bias_bw.png", units="in", width=8, height=6, res=300)
bias
dev.off()

# rel_bias_var --------------------------------------------------------
table_anova("rel_bias_var", BW_N) # method, J, ICC_jk, n_bar(medium)

rel_bias_var <- BW_N %>%
  ggplot(aes(x = J, y = rel_bias_var, fill = method, color = method)) + 
  geom_hline(yintercept = c(.9, 1.1), linetype = "dashed") +
  geom_hline(yintercept = c(1)) +
  geom_boxplot(alpha = .6, lwd = .3) + 
  facet_grid(n_bar_f ~ ES, scales = "free_y") +  
  labs(x = "Number of Schools", y = "Relative Bias of SE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#00BFC4", "#F564E3")) +
  scale_color_manual(values = c("#00BFC4", "#F564E3"))

rel_bias_var

png("rel_bias_var_bw.png", units="in", width=8, height=6, res=300)
rel_bias_var
dev.off()

# rmse --------------------------------------------------------
# Within-cluster effect
table_anova("rmse", BW_N) # assumption, method, J, n_bar, ICC_k (rather small)

rmse <- BW_N %>%
  ggplot(aes(x = J, y = rmse, fill = method, color = method)) + 
  geom_boxplot(alpha = .6, lwd = .1) + 
  facet_grid(assumption ~ ES, scales = "free_y") +  
  labs(x = "Number of Schools", y = "RMSE") + 
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.caption=element_text(hjust = 0),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_x_discrete(limits=c("20","70","150")) +
  # guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))

rmse

png("rmse_bw.png", units="in", width=8, height=6, res=300)
rmse
dev.off()
