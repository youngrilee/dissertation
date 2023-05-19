library(tidyverse)
library(lfe)
library(lme4)
library(plm)
library(insight)
rm(list = ls())

dat <- read_csv("Raudenbush & Bryk (2002) school-neighborhood.csv")
dat <- dat %>% mutate(attain = attain*10) # outcome \times 10

# Hausman Test
fixed <- plm(attain ~ p7read, data = dat, index = c("neighid", "schid"),
             model = "within")
random <- plm(attain ~ p7read, data = dat, index=c("neighid", "schid"),
              model="random")
phtest(fixed, random)

# Descriptive Statistics
per_sch <- dat %>% group_by(schid) %>% count()
per_ngh <- dat %>% group_by(neighid) %>% count()
ngh_per_sch <- dat %>% group_by(schid, neighid) %>% count() %>% 
  group_by(schid) %>% count()
desc <- data.frame(
  Variable = c("per_schools", "per_neighborhoods",
               "ngh_per_sch",
               "Education Attainment",
               "Reading"),
  mean = c(mean(per_sch$n), mean(per_ngh$n), mean(ngh_per_sch$n),
           mean(dat$attain), mean(dat$p7read)),
  sd = c(sd(per_sch$n), sd(per_ngh$n), sd(ngh_per_sch$n),
         sd(dat$attain), sd(dat$p7read)),
  min = c(min(per_sch$n), min(per_ngh$n), min(ngh_per_sch$n),
          min(dat$attain), min(dat$p7read)),
  Q1 = c(quantile(per_sch$n, c(.25)),
         quantile(per_ngh$n, c(.25)),
         quantile(ngh_per_sch$n, c(.25)),
         quantile(dat$attain, c(.25)),
         quantile(dat$p7read, c(.25))),
  Q2 = c(quantile(per_sch$n, c(.5)),
         quantile(per_ngh$n, c(.5)),
         quantile(ngh_per_sch$n, c(.5)),
         quantile(dat$attain, c(.5)),
         quantile(dat$p7read, c(.5))),
  Q3 = c(quantile(per_sch$n, c(.75)),
         quantile(per_ngh$n, c(.75)),
         quantile(ngh_per_sch$n, c(.75)),
         quantile(dat$attain, c(.75)),
         quantile(dat$p7read, c(.75))),
  max = c(max(per_sch$n), max(per_ngh$n), max(ngh_per_sch$n),
          max(dat$attain), max(dat$p7read)))

desc

# Approximate centering of p7read
cell_id <- dat %>% group_by(schid, neighid) %>% count() %>% ungroup() %>% 
  mutate(cellid = row_number()) %>% ungroup() %>% 
  dplyr::select(-n)

# calc grand/cluster/cell means
dat_cent <- 
  dat %>%
  left_join(cell_id, by = c("schid", "neighid")) %>%
  group_by(cellid) %>%
  mutate(cell_mean = mean(p7read)) %>%
  group_by(schid) %>%
  mutate(sch_mean = mean(p7read)) %>%
  group_by(neighid) %>%
  mutate(neighbor_mean = mean(p7read)) %>%
  select(neighid, schid, cellid, everything()) %>% 
  ungroup() %>%
  mutate(
    grand_mean = mean(p7read)
    # read_approx = p7read - sch_mean - neighbor_mean + grand_mean,
    # read_approx_cell = p7read - cell_mean
  )

# calc adaptively centered-cluster mean
fit <- felm(p7read ~ 0 | schid + neighid, data = dat_cent)
dat_cent$adapt_cluster <- residuals(fit)

# calc adaptively centered-cell mean
fit <- felm(p7read ~ 0 | cellid, data = dat_cent)
dat_cent$adapt_cell <- residuals(fit)

# calc adaptively calculated cell mean
fit <- felm(cell_mean ~ 0 | schid + neighid, data = dat_cent)
dat_cent$cell_mean_adapt <- as.numeric(residuals(fit))

# sparsity calculation
JK <- max(cell_id$cellid)
J <- dat_cent %>% group_by(schid) %>% count() %>% nrow()
K <- dat_cent %>% group_by(neighid) %>% count() %>% nrow()
JK/(J*K)

# No Centering
CCREM_un <- lmer(attain ~ p7read + (1 | schid) + (1 | neighid), data = dat_cent)
summary(CCREM_un)
cor(dat_cent$neighid, dat_cent$neighbor_mean)

# Grand-mean centering
CCREM_grand <- lmer(attain ~ I(p7read - grand_mean) + (1 | schid) + (1 | neighid), data = dat_cent)
summary(CCREM_grand)

# Within-RE 
CCREM_cluster <- lmer(attain ~ adapt_cluster + (1 | schid) + (1 | neighid), data = dat_cent)
summary(CCREM_cluster)

# CRE
CCREM_cluster_CRE <- lmer(attain ~ adapt_cluster + sch_mean + neighbor_mean + (1 | schid) + (1 | neighid), data = dat_cent)
summary(CCREM_cluster_CRE)

# Within-RE Cell
CCREM_cell <- lmer(attain ~ adapt_cell + (1 | schid) + (1 | neighid) + (1 | cellid), data = dat_cent)
summary(CCREM_cell) 

# CRE Cell
CCREM_cell_CRE <- lmer(attain ~ adapt_cell + sch_mean + neighbor_mean + cell_mean_adapt + (1 | schid) + (1 | neighid) + (1 | cellid), data = dat_cent)
summary(CCREM_cell_CRE)

# FE
felm(attain ~ p7read | schid + neighid, data = dat) %>%
  summary()

# hybrid model: sch-FE, nei-RE
CCREM_hybrid <- lmer(attain ~ p7read + factor(schid) + (1 | neighid), data = dat_cent)
summary(CCREM_hybrid) 

# hybrid model: sch-RE, nei-FE
CCREM_hybrid <- lmer(attain ~ p7read + factor(neighid) + (1 | schid), data = dat_cent)
summary(CCREM_hybrid) 

# variance decomposition of p7read
model <- lmer(attain ~ adapt_cell + (1 | schid) + (1 | neighid) + (1 | cellid), data = dat_cent)
summary(model)

## f1
gamma <- fixed.effects(model)["adapt_cell"] %>% as.numeric()
cov <- dat_cent$adapt_cell %>% as.matrix()
phi1 <- var(cov)
f1 <- t(gamma) %*% phi1 %*% gamma # (gamma^2) * phi1
## sch_var
sch_var <- as.numeric(unlist(VarCorr(model))["schid"])
## neigh_var
neigh_var <- as.numeric(unlist(VarCorr(model))["neighid"])
## cell_var
cell_var <- as.numeric(unlist(VarCorr(model))["cellid"])
## sigma
sigma <- sigma(model)^2
## double-check
get_variance(model)

data.frame(
  component = c("fixed", "school_random", "neighborhood_random", "cell_random", "residual"),
  var = c(f1, sch_var, neigh_var, cell_var, sigma)) %>% 
  mutate(total = sum(var),
         var_pct = var / total,
         var_pct = round(var_pct, 4))
