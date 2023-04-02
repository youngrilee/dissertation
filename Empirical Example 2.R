library(tidyverse)
library(lfe)
library(lme4)
library(plm)
rm(list = ls())

dat <- read.csv("Paterson (1991) primary-secondary.csv")
dat <- dat %>% mutate(ATTAIN = ATTAIN*10) # outcome \times 10

# Hausman Test
fixed <- plm(ATTAIN ~ VRQ, data = dat, index = c("PID", "SID"), model = "within")
random <- plm(ATTAIN ~ VRQ, data = dat, index = c("PID", "SID"), model="random")
phtest(fixed, random)

# Descriptive Statistics
per_PID <- dat %>% group_by(PID) %>% count()
per_SID <- dat %>% group_by(SID) %>% count()
PID_per_SID <- 
  dat %>% group_by(PID, SID) %>% count() %>% 
  group_by(SID) %>% count()
desc <- data.frame(
  Variable = c("per_PID", "per_SID",
               "PID_per_SID",
               "ATTAIN",
               "VRQ"),
  mean = c(mean(per_PID$n), mean(per_SID$n), mean(PID_per_SID$n),
           mean(dat$ATTAIN), mean(dat$VRQ)),
  sd = c(sd(per_PID$n), sd(per_SID$n), sd(PID_per_SID$n),
         sd(dat$ATTAIN), sd(dat$VRQ)),
  min = c(min(per_PID$n), min(per_SID$n), min(PID_per_SID$n),
          min(dat$ATTAIN), min(dat$VRQ)),
  Q1 = c(quantile(per_PID$n, c(.25)),
         quantile(per_SID$n, c(.25)),
         quantile(PID_per_SID$n, c(.25)),
         quantile(dat$ATTAIN, c(.25)),
         quantile(dat$VRQ, c(.25))),
  Q2 = c(quantile(per_PID$n, c(.5)),
         quantile(per_SID$n, c(.5)),
         quantile(PID_per_SID$n, c(.5)),
         quantile(dat$ATTAIN, c(.5)),
         quantile(dat$VRQ, c(.5))),
  Q3 = c(quantile(per_PID$n, c(.75)),
         quantile(per_SID$n, c(.75)),
         quantile(PID_per_SID$n, c(.75)),
         quantile(dat$ATTAIN, c(.75)),
         quantile(dat$VRQ, c(.75))),
  max = c(max(per_PID$n), max(per_SID$n), max(PID_per_SID$n),
          max(dat$ATTAIN), max(dat$VRQ)))

desc

# Approximate centering of VRQ
cell_id <- dat %>% group_by(PID, SID) %>% count() %>% ungroup() %>% 
  mutate(cellid = row_number()) %>% ungroup() %>% 
  dplyr::select(-n)

# calc grand/cluster/cell means
dat_cent <- 
  dat %>%
  left_join(cell_id, by = c("PID", "SID")) %>%
  group_by(cellid) %>% 
  mutate(cell_mean = mean(VRQ)) %>% 
  group_by(PID) %>%
  mutate(PID_mean = mean(VRQ)) %>%
  group_by(SID) %>%
  mutate(SID_mean = mean(VRQ)) %>%
  select(PID, SID, cellid, everything()) %>% 
  ungroup() %>%
  mutate(
    grand_mean = mean(VRQ),
    # VRQ_approx = VRQ - PID_mean - SID_mean + grand_mean,
    # VRQ_approx_cell = VRQ - cell_mean
  )

# calc adaptively centered-cluster mean
fit <- felm(VRQ ~ 0 | PID + SID, data = dat_cent)
dat_cent$adapt_cluster <- residuals(fit)

# calc adaptively centered-cell mean
fit <- felm(VRQ ~ 0 | cellid, data = dat_cent)
dat_cent$adapt_cell <- residuals(fit)

# calc adaptively calculated cell mean
fit <- felm(cell_mean ~ 0 | PID + SID, data = dat_cent)
dat_cent$cell_mean_adapt <- as.numeric(residuals(fit))


# sparsity calculation
JK <- max(cell_id$cellid)
J <- dat_cent %>% group_by(PID) %>% count() %>% nrow()
K <- dat_cent %>% group_by(SID) %>% count() %>% nrow()
JK/(J*K)

# No Centering
CCREM_un <- lmer(ATTAIN ~ VRQ + (1 | PID) + (1| SID), data = dat_cent)
summary(CCREM_un)

# Grand-mean centering
CCREM_grand <- lmer(ATTAIN ~ I(VRQ - grand_mean) + (1 | PID) + (1 | SID), data = dat_cent)
summary(CCREM_grand)

# Within-RE 
CCREM_adapt <- lmer(ATTAIN ~ adapt_cluster + (1 | PID) + (1 | SID), data = dat_cent)
summary(CCREM_adapt)

# CRE
CCREM_adapt <- lmer(ATTAIN ~ adapt_cluster + PID_mean + SID_mean + (1 | PID) + (1 | SID), data = dat_cent)
summary(CCREM_adapt)

# Within-RE Cell
CCREM_adapt <- lmer(ATTAIN ~ adapt_cell + (1 | PID) + (1 | SID) + (1 | cellid), data = dat_cent)
summary(CCREM_adapt) 

# CRE Cell
CCREM_WB <- lmer(ATTAIN ~ adapt_cell + PID_mean + SID_mean + cell_mean_adapt + (1 | PID) + (1 | SID) , data = dat_cent)
summary(CCREM_WB)

# FE
felm(ATTAIN ~ VRQ | PID + SID, data = dat_cent) %>% summary()

# hybrid model: PID-FE, SID-RE
CCREM_hybrid <- lmer(ATTAIN ~ VRQ + factor(PID) + (1 | SID), data = dat_cent)
summary(CCREM_hybrid) 

# hybrid model: PID-RE, SID-FE
CCREM_hybrid <- lmer(ATTAIN ~ VRQ + factor(SID) + (1 | PID), data = dat_cent)
summary(CCREM_hybrid) 

# Adaptive centering of VRQ
# VRQ_fit <- felm(VRQ ~ 0 | SID, data = dat)
# dat_cent$VRQ_adapt <- residuals(VRQ_fit)
# CCREM_hybrid <- lmer(ATTAIN ~ VRQ_adapt + factor(PID) + (1 | SID), data = dat_cent)
# summary(CCREM_hybrid) 

# Adaptive centering of VRQ
# VRQ_fit <- felm(VRQ ~ 0 | PID + SID, data = dat)
# dat_cent$VRQ_adapt <- residuals(VRQ_fit)
# CCREM_hybrid <- lmer(ATTAIN ~ VRQ_adapt + factor(PID) + (1 | SID), data = dat_cent)
# summary(CCREM_hybrid) 

