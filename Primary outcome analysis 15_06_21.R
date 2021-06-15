library(haven)
library(tidyverse)
library(psych)
library(nlme)
library(sjlabelled)
library(sjPlot)
library(ggeffects)
library(extrafont)
library(brms)
loadfonts(device = "win")

nac <- read_sav("C:\\Users\\lachy\\OneDrive - The University of Melbourne\\Unimelb\\NAC placebo\\datafiles\\COGNAC data file (2) MERGED 28-Apr-2021.sav")



#### Primary outcome ####

## check YBOCS variables for errors/outliers
## recalculate total scores 

bl_ybocs <- names(nac %>% 
        select(starts_with("YBOCS_BL2_"), -YBOCS_BL2_total, -YBOCS_BL2_11_insight))


nac$YBOCS_BL2_total2 <- rowSums(nac[,bl_ybocs])



## Week 4

w4_ybocs <- names(nac %>% 
        select(starts_with("YBOCS_W4_"), -YBOCS_W4_total, -YBOCS_W4_11_insight))

rowSums(is.na(nac[,w4_ybocs])) # check for missing items

nac$YBOCS_W4_total2 <- rowSums(nac[,w4_ybocs])



## week 8

w8_ybocs <- names(nac %>% 
                    select(starts_with("YBOCS_W8_"), -YBOCS_W8_total, -YBOCS_W8_11_insight))

rowSums(is.na(nac[,w8_ybocs])) # check for missing items


nac$YBOCS_W8_total2 <- rowSums(nac[,w8_ybocs])



## week 12


w12_ybocs <- names(nac %>% 
                    select(starts_with("YBOCS_W12_"), 
                           -YBOCS_W12_total, -YBOCS_W12_11_insight, -YBOCS_W12_severity))

rowSums(is.na(nac[,w12_ybocs])) # check for missing items

nac$YBOCS_W12_total2 <- rowSums(nac[,w12_ybocs])


## week 16

w16_ybocs <- names(nac %>% 
                     select(starts_with("YBOCS_W16_"), 
                            -YBOCS_W16_total, -YBOCS_W16_severity, -YBOCS_W16_insight))

rowSums(is.na(nac[,w16_ybocs])) # check for missing items

nac$YBOCS_W16_total2 <- rowSums(nac[,w16_ybocs])



## week 20

w20_ybocs <- names(nac %>% 
                     select(starts_with("YBOCS_W20_"), 
                            -YBOCS_W20_total, -YBOCS_W20_severity, -YBOCS_W20_insight))

rowSums(is.na(nac[,w20_ybocs])) # check for missing items

nac$YBOCS_W20_total2 <- rowSums(nac[,w20_ybocs])



## Treatment variable 

nac$Group_allocation <- as_label(nac$Group_allocation)

nac$Group_allocation <- relevel(nac$Group_allocation, "B") # placebo reference group


### Change data to long format ###


ybocs_vars <- c("YBOCS_W4_total2", "YBOCS_W8_total2", 
                "YBOCS_W12_total2","YBOCS_W16_total2", "YBOCS_W20_total2")

nac_long <- nac %>% 
  pivot_longer(cols = all_of(ybocs_vars), 
               names_to = "time", 
               values_to = "YBOCS")

nac_long <- nac_long %>% 
  mutate(time = if_else(time == "YBOCS_W4_total2", 0,
                              if_else(time == "YBOCS_W8_total2", 1,
                                      if_else(time == "YBOCS_W12_total2", 2, 
                                              if_else(time == "YBOCS_W16_total2", 3, 4)))))




#### linear mixed models 

## test linearity (continuous time vs log time vs categorical time) ##

nac_long$cat_time <- as.factor(nac_long$time) # categorical time

nac_long$log_time <- log(nac_long$time + 1)


# LR test to compare fit 

m1 <- lme(YBOCS ~ YBOCS_BL2_total2 + time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long, method = "ML") # best model 

m2 <- lme(YBOCS ~ YBOCS_BL2_total2 + log_time,
          random = ~ log_time | Participant_ID,
          na.action = na.omit, data = nac_long, method = "ML")

m3 <- lme(YBOCS ~ YBOCS_BL2_total2 + cat_time,
          random = ~ 1 | Participant_ID,
          na.action = na.omit, data = nac_long, method = "ML")

anova(m1, m2, m3)


## testing covariance structure 

m1 <- lme(YBOCS ~ YBOCS_BL2_total2 + time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long) # best model 

m2 <- lme(YBOCS ~ YBOCS_BL2_total2 + time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long,
          correlation = corAR1()) # autoregressive

m3 <- lme(YBOCS ~ YBOCS_BL2_total2 + time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long,
          correlation = corCompSymm()) # compound symmetry

anova(m1, m2, m3)


# zero time variable at final visit

nac_long$ctime <- nac_long$time - 4



# primary outcome model (conditional growth model) 


model1 <- lme(YBOCS ~ YBOCS_BL2_total2 + time*Group_allocation,
              random = ~ time | Participant_ID,
              na.action = na.omit, data = nac_long)

summary(model1)

em <- emmeans(model1,
              pairwise ~ Group_allocation | time, 
              at = list(time = 4))

em

#### Figure 1 ####


preds <- ggpredict(model1, 
                   terms = c("time","Group_allocation","YBOCS_BL2_total2 [22.5]"))


# add baseline data to predictions

bl <- data.frame(x = c(-1,-1), 
                 predicted = c(22.5, 22.5), 
                 std.error = c(NA,NA),
                 conf.low = c(NA,NA),
                 conf.high = c(NA,NA),
                 group = c("B","A"),
                 facet = c(22.5,22.5))

preds <- rbind(bl, preds)

preds$Treatment <- as.factor(preds$group)

levels(preds$Treatment) <- c("N-acetylcysteine", "Placebo")



# set plot theme 

theme_set(theme_linedraw() +
            theme(text = element_text(family = "Times New Roman"),
                  panel.grid = element_blank()))

# plot 

ggplot(preds, aes(x = x, y = predicted, colour = Treatment)) +
  geom_line(aes(linetype = Treatment)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0.1),
                  fatten = 2) +
  scale_x_continuous(labels = c("BL","W4","W8","W12","W16","W20")) +
  xlab("Visit") + ylab("YBOCS ± 95% CI") + 
  ylim(c(12.5,25)) +
  scale_color_manual(values=c("#999999", "#E69F00"))

# save plot to png
ggsave("NAC ybocs.png", device = "png", path = "C:/Users/lachy/OneDrive - The University of Melbourne/Unimelb/NAC placebo/Plots")




#### sensitivity pattern mixture model ####

# create drop_out indicators
# drop_early for dropout before W12, drop_late for dropout W12 or later

nac <- nac %>% 
  mutate(drop_early = if_else(is.na(YBOCS_W4_total2) | is.na(YBOCS_W8_total2) &
                                is.na(YBOCS_W12_total2) & is.na(YBOCS_W16_total2) &
                                is.na(YBOCS_W20_total2), 1, 0)) %>% 
  mutate(drop_late = if_else(drop_early == 0 & is.na(YBOCS_W20_total2), 1, 0))
                                                                 

# long format 

nac_long <- nac %>% 
  pivot_longer(cols = all_of(ybocs_vars), 
               names_to = "time", 
               values_to = "YBOCS")

nac_long <- nac_long %>% 
  mutate(time = if_else(time == "YBOCS_W4_total2", 0,
                        if_else(time == "YBOCS_W8_total2", 1,
                                if_else(time == "YBOCS_W12_total2", 2, 
                                        if_else(time == "YBOCS_W16_total2", 3, 4)))))


## Bayesian pattern mixture model

# priors 

priors <- prior(normal(0,2), class = b) # weak shrinkage priors so model is estimable

pm1 <- brm(YBOCS ~ YBOCS_BL2_total2 + time*Group_allocation*drop_early + 
             time*Group_allocation*drop_late + (1 + time | Participant_ID),
           data = nac_long,
           prior = priors,
           seed = 1,
           cores = 2, chains = 4, iter = 3500, warmup = 1000,
           control = list(adapt_delta = .95))

# average over missing data patterns 

em1 <- emmeans(pm1,
               pairwise ~ Group_allocation | time, 
               CIs = F,
               at = list(time = 4),
               weights = "proportional")

em1


## simple pattern mixture
# Just one dropout indicator (dropout or no)

nac_long$drop_any <- ifelse(nac_long$drop_early == 1 | nac_long$drop_late == 1, 1, 0)

# model
pm2 <- lme(YBOCS ~ YBOCS_BL2_total2 + time*Group_allocation*drop_any,
              random = ~ time | Participant_ID,
              na.action = na.omit, data = nac_long)

summary(pm2)

# average over dropout pattern
em <- emmeans(pm2,
              pairwise ~ Group_allocation | time, 
              at = list(time = 4),
              weights = "proportional")

em
