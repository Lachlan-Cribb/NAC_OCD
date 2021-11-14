# Load required packages
library(haven)
library(tidyverse)
library(lme4)
library(sjPlot)
library(psych)
library(ggeffects)
library(extrafont)
library(emmeans)
library(mice)
library(qwraps2)
library(brms)
loadfonts(device = "win")

# read data

nac <- read_sav("PATH HERE")

# remove -9 and -999 SPSS missing data indicators 

nac <- nac %>% 
  mutate(across(where(is.numeric), ~ na_if(., -9))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -99)))

# drop unneeded variable with problematic formatting

nac <- nac[,-890]

## Attended any follow-up? ##

nac <- nac %>% 
  mutate(any_fu = if_else(is.na(YBOCS_W4_total) & is.na(YBOCS_W8_total) &
                            is.na(YBOCS_W12_total) & is.na(YBOCS_W16_total) & 
                            is.na(YBOCS_W20_total), 0, 1)) %>% 
  mutate(any_fu = factor(any_fu, labels = c("No", "Yes")))


## Retain only participants who attended at least one follow-up 

nac <- nac %>% 
  filter(any_fu == "Yes")

## Create drop-out variable 

nac <- nac %>% 
  mutate(dropout = if_else(is.na(YBOCS_W20_total), 1, 0))

## Treatment variable 

nac$Group_allocation <- sjlabelled::as_label(nac$Group_allocation)

nac$Group_allocation <- relevel(nac$Group_allocation, "B") # placebo reference group

levels(nac$Group_allocation) <- c("Placebo","NAC")


#### Compliance ####

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(median = median(Compliance_percentage_overall, na.rm=T), 
            IQR = IQR(Compliance_percentage_overall, na.rm=T))

nac %>% 
  filter(!is.na(Compliance_percentage_overall)) %>% 
  mutate(good_compliance = if_else(Compliance_percentage_overall > 75, 1, 0)) %>% 
  group_by(good_compliance) %>% 
  tally() %>% 
  mutate(percent = n /sum(n) * 100)


## Titration

nac %>% 
  select(contains("Titration"))

nac <- nac %>% 
  rowwise() %>% 
  mutate(total_titrations = sum(c(Titration_W8,Titration_W12,Titration_W16), na.rm=T)) %>% 
  mutate(any_titration = if_else(total_titrations > 0, "Yes", "No")) %>% 
  mutate(any_titration = factor(any_titration)) %>% 
  ungroup()


nac %>% 
  group_by(Group_allocation, any_titration) %>%
  tally()

# Does titration differ between the two groups?

chisq.test(x = nac$Group_allocation, y = nac$any_titration)


## Check for outliers and missing data in covariates 

# alcohol

sum(is.na(nac$BL_standard_drinks_week))

hist(nac$BL_standard_drinks_week, breaks = 20)

# gender 

nac$Gender <- as.factor(nac$Gender)

levels(nac$Gender) <- c("Male","Female")

nac$Gender <- fct_rev(nac$Gender) # Female as reference category 

# age

sum(is.na(nac$Age))

hist(nac$Age)

# Weight

sum(is.na(nac$Weight_kg_BL)) # 2 missing obs

# smoking

sum(is.na(nac$BL_cigarettes_amount))

hist(nac$BL_cigarettes_amount, breaks = 20)

nac %>% 
  mutate(smoking = if_else(BL_cigarettes_amount > 0, 1, 0)) %>% 
  group_by(smoking) %>% 
  tally()

# site

nac$Site <- as.factor(nac$Site)

levels(nac$Site) <- c("TMC","RBH","NICM")

# standardize continuous covariates to Z scores

nac$zAge <- scale(nac$Age)

nac$zBL_standard_drinks_week <- scale(nac$BL_standard_drinks_week)

### Comorbid psychiatric illness 

nac %>%
  select(-OCD) %>% 
  rowwise() %>% 
  mutate(
    total_comorbid = sum(c_across(MDE_current:Gambling), na.rm=T)) %>%
  mutate(any_comorbid = if_else(total_comorbid > 0, 1, 0)) %>% 
  group_by(any_comorbid) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

nac %>% 
  group_by(GAD) %>% 
  tally() %>% 
  mutate(p = n/sum(n) * 100)

median(nac$Chronicity_OCD_dx_years)


#### Primary outcome ####

## recalculate total scores and check for missing data

bl_ybocs <- nac %>% 
  select(starts_with("YBOCS_BL2_"), -contains("total"), -YBOCS_BL2_11_insight) %>% 
  names()

nac$YBOCS_BL2_total2 <- rowSums(nac[,bl_ybocs])

w4_ybocs <- nac %>% # week 4
  select(starts_with("YBOCS_W4_"), -contains("total"), -YBOCS_W4_11_insight) %>% 
  names()

rowSums(is.na(nac[,w4_ybocs])) # check for missing items

nac$YBOCS_W4_total2 <- rowSums(nac[,w4_ybocs])

w8_ybocs <- nac %>% # week 8
  select(starts_with("YBOCS_W8_"), -contains("total"), -YBOCS_W8_11_insight) %>% 
  names()

rowSums(is.na(nac[,w8_ybocs])) # check for missing items

nac$YBOCS_W8_total2 <- rowSums(nac[,w8_ybocs])

w12_ybocs <- nac %>% # week 12
  select(starts_with("YBOCS_W12_"), -contains("total"), -YBOCS_W12_11_insight, -YBOCS_W12_severity) %>% 
  names()

rowSums(is.na(nac[,w12_ybocs])) # check for missing items

nac$YBOCS_W12_total2 <- rowSums(nac[,w12_ybocs])

w16_ybocs <- nac %>% # week 16
  select(starts_with("YBOCS_W16_"), -contains("total"), -YBOCS_W16_severity, -YBOCS_W16_insight) %>% 
  names()

rowSums(is.na(nac[,w16_ybocs])) # check for missing items

nac$YBOCS_W16_total2 <- rowSums(nac[,w16_ybocs])

w20_ybocs <- nac %>% # week 20
  select(starts_with("YBOCS_W20_"), -contains("total"), -YBOCS_W20_severity, -YBOCS_W20_insight) %>% 
  names()

rowSums(is.na(nac[,w20_ybocs])) # check for missing items

nac$YBOCS_W20_total2 <- rowSums(nac[,w20_ybocs])


# raw means #

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c("YBOCS_BL2_total2", "YBOCS_W20_total2"),
                   list(mean = ~ mean(., na.rm=T), 
                        sd = ~ sd(., na.rm=T))))

## Histogram of YBOCS variables

nac %>% 
  pivot_longer(cols = ends_with("total2"), 
               values_to = "YBOCS", names_to = "Visit") %>% 
  mutate(Visit = factor(Visit, levels = unique(Visit))) %>% 
  ggplot(aes(x = YBOCS, fill = Group_allocation)) +
  geom_histogram(aes(x = YBOCS, colour = Group_allocation), 
               alpha = 0.5, bins = 15) +
  facet_wrap(~ Visit) +
  labs(x = "Total YBOCS") +
  ggthemes::theme_stata() +
  ggthemes::scale_fill_colorblind() +
  ggthemes::scale_color_colorblind()
  
## data to long format ##

ybocs_vars <- c("YBOCS_W4_total2", "YBOCS_W8_total2", 
                "YBOCS_W12_total2","YBOCS_W16_total2", "YBOCS_W20_total2")

nac_long <- nac %>% 
  pivot_longer(cols = all_of(ybocs_vars), 
               names_to = "time", 
               values_to = "YBOCS")

nac_long$time <- factor(nac_long$time, levels = unique(nac_long$time))

nac_long$time <- as.numeric(nac_long$time) - 1

### linear mixed models 

## test linearity (continuous time vs log time vs categorical time) ##

nac_long$cat_time <- as.factor(nac_long$time) # categorical time

nac_long$log_time <- log(nac_long$time + 1)

# LR test to compare fit 

m1 <- lmer(YBOCS ~ YBOCS_BL2_total2*time + Group_allocation*time + # linear
          (1 + time | Participant_ID),
          data = nac_long, REML = F)

m2 <- lmer(YBOCS ~ YBOCS_BL2_total2*log_time + Group_allocation*log_time + # log
          (1 + log_time | Participant_ID), 
          data = nac_long, REML = F)

m3 <- lmer(YBOCS ~ YBOCS_BL2_total2*cat_time + Group_allocation*cat_time + 
          (1 | Participant_ID), # categorical
          data = nac_long, REML = F)

anova(m1, m2, m3)

## testing covariance structure 

m1 <- nlme::lme(YBOCS ~ YBOCS_BL2_total2*time + Group_allocation*time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long) # none

m2 <- nlme::lme(YBOCS ~ YBOCS_BL2_total2*time + Group_allocation*time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long,
          correlation = corAR1()) # autoregressive

m3 <- nlme::lme(YBOCS ~ YBOCS_BL2_total2*time + Group_allocation*time,
          random = ~ time | Participant_ID,
          na.action = na.omit, data = nac_long,
          correlation = corCompSymm()) # compound symmetry

anova(m1, m2, m3)

## primary outcome model (conditional growth model) 

nac_long$zYBOCS_BL2_total2 <- scale(nac_long$YBOCS_BL2_total2)

model1 <- lmer(YBOCS ~ zYBOCS_BL2_total2 + time*Group_allocation + Site + zAge +
               Gender + zBL_standard_drinks_week + (1 + time | Participant_ID),
               data = nac_long)

summary(model1)

emmeans(model1,
        "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

## Cohen's d

c(0.525, -2.18, 3.23) / sd(nac$YBOCS_BL2_total2)


#### Figure 1 ####

preds <- ggemmeans(model1, 
                   terms = c("time","Group_allocation"))

preds

# add baseline data to predictions

bl <- data.frame(x = c(-1,-1), 
                 predicted = c(22.6, 22.6), # baseline mean
                 std.error = c(NA,NA),
                 conf.low = c(NA,NA),
                 conf.high = c(NA,NA),
                 group = c("Placebo","NAC"))

preds <- rbind(bl, preds)

# plot 

ggplot(preds, aes(x = x, y = predicted, colour = group)) +
  geom_line() + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high,
                      shape = group), 
                  position = position_dodge(width = 0.2),
                  fatten = 4) +
  scale_x_continuous(labels = c("BL","W4","W8","W12","W16","W20")) +
  xlab("Visit") + ylab("Total YBOCS, 95% CI") + 
  ylim(c(12.5,25)) +
  ggthemes::scale_colour_colorblind() +
  ggthemes::theme_stata()

# save plot to png

ggsave("NAC ybocs v2.png", device = "png", 
       width = 5, height = 5.5,
       path = "C:/Users/lachy/OneDrive - The University of Melbourne/Unimelb/NAC placebo/Plots")

#### responder analysis ####

## response 

response_dat <- nac %>% 
  filter(!is.na(YBOCS_W20_total2)) %>% 
  mutate(change_YBOCS_percent = (YBOCS_W20_total2 - YBOCS_BL2_total2) / YBOCS_BL2_total2 * 100)  %>% 
  mutate(response = if_else(change_YBOCS_percent <= -35 & CGI_I_W20 %in% c(1,2), 1, 0))
  
response_dat %>%
  group_by(Group_allocation, response) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# risk ratio for response

response_mod <- glm(response ~ Group_allocation, data = response_dat, family = binomial(link = "log"))

tidy(response_mod, conf.int =T) %>% 
  mutate(RR = exp(estimate), 
         RR_low = exp(conf.low), 
         RR_high = exp(conf.high)) %>% 
  select(term, RR, RR_low, RR_high, p.value)


## Jacobson and Truax

# extract estimates from linear mixed model

preds <- predict(model1, newdata = nac_long) %>% 
  as_tibble()

nac_long2 <- cbind(nac_long, preds) # merge predictions into data

# Calculate model predicted week 20 score - baseline

nac_long2 <- nac_long2 %>% 
  rowwise() %>% 
  mutate(estimated_change = value - YBOCS_BL2_total2) %>% 
  ungroup()

# response

nac_long2 <- nac_long2 %>% 
  filter(time == 4) %>% 
  mutate(JT_response = if_else(estimated_change <= -10 & value <= 14, 1, 0))

nac_long2 %>% 
  filter(!is.na(JT_response)) %>% 
  group_by(Group_allocation, JT_response) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)  

# Risk ratio

JT_mod <- glm(JT_response ~ Group_allocation, data = nac_long2, 
              family = binomial(link = "log"))

tidy(JT_mod, conf.int =T) %>% 
  mutate(RR = exp(estimate), 
         RR_low = exp(conf.low), 
         RR_high = exp(conf.high)) %>% 
  select(term, RR, RR_low, RR_high, p.value)


#### sensitivity analyses ####

### Pattern mixture model ###

pm1 <- lmer(YBOCS ~ zYBOCS_BL2_total2 + time*Group_allocation*dropout + Site +
             zAge + Gender + zBL_standard_drinks_week + (time | Participant_ID),
            data = nac_long)

summary(pm1)

# average over dropout pattern

emmeans(pm1,
        "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

### Complier only analysis ###

comp_mod <- lmer(YBOCS ~ zYBOCS_BL2_total2 + time*Group_allocation + zAge + Site +
                Gender + zBL_standard_drinks_week + (1 + time | Participant_ID),
                data = filter(nac_long, Compliance_percentage_overall > 75))

summary(comp_mod)

emmeans(comp_mod,
        "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

#### YBOCS obsessions ####

# create total scores

bl_obs <- nac %>% select(YBOCS_BL2_01:YBOCS_BL2_05) %>% names()

nac$obs_BL2_total <- rowSums(nac[,bl_obs])

W4_obs <- nac %>% select(YBOCS_W4_01:YBOCS_W4_05) %>% names()

nac$obs_W4_total <- rowSums(nac[,W4_obs])

W8_obs <- nac %>% select(YBOCS_W8_01:YBOCS_W8_05) %>% names()

nac$obs_W8_total <- rowSums(nac[,W8_obs])

W12_obs <- nac %>% select(YBOCS_W12_01:YBOCS_W12_05) %>% names()

nac$obs_W12_total <- rowSums(nac[,W12_obs])

W16_obs <- nac %>% select(YBOCS_W16_01:YBOCS_W16_05) %>% names()

nac$obs_W16_total <- rowSums(nac[,W16_obs])

W20_obs <- nac %>% select(YBOCS_W20_01:YBOCS_W20_05) %>% names()

nac$obs_W20_total <- rowSums(nac[,W20_obs])


## Check histograms

nac %>% 
  select(starts_with("obs_")) %>% 
  multi.hist()

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(obs_BL2_total, obs_W20_total),
            list(mean = ~ mean(.,na.rm=T),
                 sd = ~ sd(., na.rm=T))))


# long format data

obs_vars <- c("obs_W4_total", "obs_W8_total", "obs_W12_total","obs_W16_total","obs_W20_total")

obs_long <- nac %>% 
  pivot_longer(cols = all_of(obs_vars), 
               names_to = "time", 
               values_to = "obs")

obs_long$time <- as.numeric(factor(obs_long$time, levels = unique(obs_long$time))) - 1

# Model

obs_long$zobs_BL2_total <- scale(obs_long$obs_BL2_total)

obs_mod <- lmer(obs ~ zobs_BL2_total + time*Group_allocation + zAge + Site +
                Gender + zBL_standard_drinks_week + (time | Participant_ID),
                data = obs_long)

summary(obs_mod)

emmeans(obs_mod, 
        "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

# cohens d

c(0.0566, -1.33, 1.44) / sd(nac$obs_BL2_total, na.rm=T)


#### YBOCS compulsions ####

# create total scores

bl_comp <- nac %>% select(YBOCS_BL2_06:YBOCS_BL2_10) %>% names()

nac$comp_BL2_total <- rowSums(nac[,bl_comp])

W4_comp <- nac %>% select(YBOCS_W4_06:YBOCS_W4_10) %>% names()

nac$comp_W4_total <- rowSums(nac[,W4_comp])

W8_comp <- nac %>% select(YBOCS_W8_06:YBOCS_W8_10) %>% names()

nac$comp_W8_total <- rowSums(nac[,W8_comp])

W12_comp <- nac %>% select(YBOCS_W12_06:YBOCS_W12_10) %>% names()

nac$comp_W12_total <- rowSums(nac[,W12_comp])

W16_comp <- nac %>% select(YBOCS_W16_06:YBOCS_W16_10) %>% names()

nac$comp_W16_total <- rowSums(nac[,W16_comp])

W20_comp <- nac %>% select(YBOCS_W20_06:YBOCS_W20_10) %>% names()

nac$comp_W20_total <- rowSums(nac[,W20_comp])

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(comp_BL2_total, comp_W20_total),
                   list(mean = ~ mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

# long format data

comp_vars <- c("comp_W4_total", "comp_W8_total", "comp_W12_total","comp_W16_total","comp_W20_total")

comp_long <- nac %>% 
  pivot_longer(cols = all_of(comp_vars), 
               names_to = "time", 
               values_to = "comp")

comp_long$time <- as.numeric(factor(comp_long$time, levels = unique(comp_long$time))) - 1

# Model

comp_long$zcomp_BL2_total <- scale(comp_long$comp_BL2_total)

comp_mod <- lmer(comp ~ zcomp_BL2_total + time*Group_allocation + zAge + Site +
                 Gender + zBL_standard_drinks_week + (time | Participant_ID),
                 data = comp_long)

summary(comp_mod)

emmeans(comp_mod, "Group_allocation", 
       at = list(time = 4)) %>% 
  pairs(reverse = T, infer = T)

# Cohen's d

c(0.598, -0.908, 2.1) / sd(nac$comp_BL2_total, na.rm=T)


#### DOCS ####

# imputation function

pmm_impute <- function(vars, max_missing, data){
  
  skip <- which(rowSums(is.na(data[,vars]))>max_missing) # maximum missing for imputation
  
  imp <- mice(data = data[-skip,vars], m = 1, method = "pmm", seed = 100) # impute
  
  imp <- complete(imp) # save as dataframe 
  
  data[-skip,vars] <- imp # replace data with imputed data 
  
  data
}


# remove -3 (SPSS missing data indicator)

nac <- nac %>% mutate(across(contains("DOCS"), ~ na_if(., -3)))

# check distributions
nac %>% 
  select(starts_with("DOCS_BL")) %>% 
  multi.hist(global = F)

nac %>% 
  select(starts_with("DOCS_W4")) %>% 
  multi.hist(global = F)

nac %>% 
  select(starts_with("DOCS_W8")) %>% 
  multi.hist(global = F)

nac %>% 
  select(starts_with("DOCS_W12")) %>% 
  multi.hist(global = F)

nac %>% 
  select(starts_with("DOCS_W16")) %>% 
  multi.hist(global = F)

nac %>% 
  select(starts_with("DOCS_W20")) %>% 
  multi.hist(global = F)
  
# remove outlier 

nac <- nac %>% mutate(across(contains("DOCS"), ~na_if(., 101)))

## Impute missing DOCS items ##

# remove labels for MICE imputation

nac_imp <- nac %>% 
  mutate(across(contains("DOCS"), haven::zap_labels))


# BL imputation - skipping participants with more than 50% missing #

bl_docs <- nac %>% 
  select(DOCS_BL_contamination_01:DOCS_BL_symmetry_05) %>% 
  select(-DOCS_BL_contamination_total, -DOCS_BL_harm_total, -DOCS_BL_taboo_total) %>% 
  names()

nac_imp <- pmm_impute(vars = bl_docs, max_missing = 12, data = nac_imp)

# week 4

W4_docs <- nac %>% 
  select(DOCS_W4_contamination_01:DOCS_W4_symmetry_05) %>% 
  select(-DOCS_W4_contamination_total, -DOCS_W4_harm_total, -DOCS_W4_taboo_total) %>% 
  names()

nac_imp <- pmm_impute(vars = W4_docs, max_missing = 12, data = nac_imp)

# week 8 

W8_docs <- nac %>% 
  select(DOCS_W8_contamination_01:DOCS_W8_symmetry_05) %>% 
  select(-DOCS_W8_contamination_total, -DOCS_W8_harm_total, -DOCS_W8_taboo_total) %>% 
  names()

nac_imp <- pmm_impute(vars = W8_docs, max_missing = 12, data = nac_imp)

# week 12

W12_docs <- nac %>% 
  select(DOCS_W12_contamination_01:DOCS_W12_symmetry_05) %>% 
  select(-DOCS_W12_contamination_total, -DOCS_W12_harm_total, -DOCS_W12_taboo_total) %>% 
  names()

nac_imp <- pmm_impute(vars = W12_docs, max_missing = 12, data = nac_imp)

# Week 16 

W16_docs <- nac %>% 
  select(DOCS_W16_contamination_01:DOCS_W16_symmetry_05) %>% 
  select(-DOCS_W16_contamination_total, -DOCS_W16_harm_total, -DOCS_W16_taboo_total) %>% 
  names()

nac_imp <- pmm_impute(vars = W16_docs, max_missing = 12, data = nac_imp)

# week 20

W20_docs <- nac %>% 
  select(DOCS_W20_contamination_01:DOCS_W20_symmetry_05) %>% 
  select(-DOCS_W20_contamination_total, -DOCS_W20_harm_total, -DOCS_W20_taboo_total) %>% 
  names()

nac_imp <- pmm_impute(vars = W20_docs, max_missing = 12, data = nac_imp)

## recalculate total scores

nac_imp$DOCS_BL_total2 <- rowSums(nac_imp[,bl_docs])

nac_imp$DOCS_W4_total2 <- rowSums(nac_imp[,W4_docs])

nac_imp$DOCS_W8_total2 <- rowSums(nac_imp[,W8_docs])

nac_imp$DOCS_W12_total2 <- rowSums(nac_imp[,W12_docs])

nac_imp$DOCS_W16_total2 <- rowSums(nac_imp[,W16_docs])

nac_imp$DOCS_W20_total2 <- rowSums(nac_imp[,W20_docs])

# raw means

nac_imp %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(DOCS_BL_total2, DOCS_W20_total2),
                   list(mean = ~ mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

## DOCS mixed model 

# data to long format

docs_vars <- c("DOCS_W4_total2", "DOCS_W8_total2", "DOCS_W12_total2", "DOCS_W16_total2", "DOCS_W20_total2")

nac_long_imp <- nac_imp %>% 
  pivot_longer(cols = all_of(docs_vars), 
               names_to = "time", 
               values_to = "DOCS")

# check DOCS histogram

hist(nac_long_imp$DOCS, breaks = 30)

# time variable

nac_long_imp$time <- factor(nac_long_imp$time, levels = unique(nac_long_imp$time))

nac_long_imp$time <- as.numeric(nac_long_imp$time) - 1

# model 

docs_m1 <- lmer(DOCS ~ DOCS_BL_total2 + time*Group_allocation + zAge + Site +
               Gender + zBL_standard_drinks_week + (time | Participant_ID),
               data = nac_long_imp)

summary(docs_m1)

emmeans(docs_m1, "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

# Cohen's d

c(-1.08, -5.55, 3.39) / sd(nac_imp$DOCS_BL_total2, na.rm=T)


#### SIGHD ####

# remove -3 (missing)

nac <- nac %>% mutate(across(contains("SIGHD"), ~na_if(., -3)))

# check histograms for outliers

nac %>% 
  select(contains("SIGHD")) %>% 
  select(contains("BL"), contains("W4")) %>% 
  multi.hist(global = F)

nac %>% 
  select(contains("SIGHD")) %>% 
  select(contains("W8"), contains("W12")) %>% 
  multi.hist(global = F)

nac %>% 
  select(contains("SIGHD")) %>% 
  select(contains("W16"), contains("W20")) %>% 
  multi.hist(global = F)

# missed items at BL, W4, W8, W20
# impute missing items

nac_imp2 <- nac %>% 
  mutate(across(contains("SIGHD"), haven::zap_labels))

# BL

bl_sighd <- nac_imp2 %>% select(SIGHD_BL_01:SIGHD_BL_17) %>% names()

nac_imp2 <- pmm_impute(vars = bl_sighd, max_missing = 8, data = nac_imp2)

# W4

W4_sighd <- nac_imp2 %>% select(SIGHD_W4_01:SIGHD_W4_17) %>% names()

nac_imp2 <- pmm_impute(vars = W4_sighd, max_missing = 8, data = nac_imp2)

# W8

W8_sighd <- nac_imp2 %>% select(SIGHD_W8_01:SIGHD_W8_17) %>% names()

nac_imp2 <- pmm_impute(vars = W8_sighd, max_missing = 8, data = nac_imp2)

# W12

W12_sighd <- nac_imp2 %>% select(SIGHD_W12_01:SIGHD_W12_17) %>% names()

nac_imp2 <- pmm_impute(vars = W12_sighd, max_missing = 8, data = nac_imp2)

# W16

W16_sighd <- nac_imp2 %>% select(SIGHD_W16_01:SIGHD_W16_17) %>% names()

nac_imp2 <- pmm_impute(vars = W16_sighd, max_missing = 8, data = nac_imp2)

# W20

W20_sighd <- nac_imp2 %>% select(SIGHD_W20_01:SIGHD_W20_17) %>% names()

nac_imp2 <- pmm_impute(vars = W20_sighd, max_missing = 8, data = nac_imp2)

## recalculate total scores

nac_imp2$SIGHD_BL_total2 <- rowSums(nac_imp2[,bl_sighd])

nac_imp2$SIGHD_W4_total2 <- rowSums(nac_imp2[,W4_sighd])

nac_imp2$SIGHD_W8_total2 <- rowSums(nac_imp2[,W8_sighd])

nac_imp2$SIGHD_W12_total2 <- rowSums(nac_imp2[,W12_sighd])

nac_imp2$SIGHD_W16_total2 <- rowSums(nac_imp2[,W16_sighd])

nac_imp2$SIGHD_W20_total2 <- rowSums(nac_imp2[,W20_sighd])

# raw means

nac_imp2 %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(SIGHD_BL_total2, SIGHD_W20_total2),
                   list(mean = ~mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

## long format 

sig_vars <- c("SIGHD_W4_total2", "SIGHD_W8_total2", "SIGHD_W12_total2", "SIGHD_W16_total2", "SIGHD_W20_total2")

nac_long_imp2 <- nac_imp2 %>% 
  pivot_longer(cols = all_of(sig_vars), 
               names_to = "time", 
               values_to = "SIGHD")

# check SIGHD variable for outliers

hist(nac_long_imp2$SIGHD, breaks = 30)

# time variable

nac_long_imp2$time <- factor(nac_long_imp2$time, levels = unique(nac_long_imp2$time))

nac_long_imp2$time <- as.numeric(nac_long_imp2$time) - 1

# model 

nac_long_imp2$zSIGHD_BL_total2 <- scale(nac_long_imp2$SIGHD_BL_total2)

sighd_m1 <- lmer(SIGHD ~ zSIGHD_BL_total2 + time*Group_allocation + zAge + Site +
                 Gender + zBL_standard_drinks_week + (1 | Participant_ID), 
                 data = nac_long_imp2)

summary(sighd_m1)

emmeans(sighd_m1, "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

# Cohen's d

c(-0.807, -2.64, 1.03) / sd(nac_imp2$SIGHD_BL_total2, na.rm=T)


#### BAI ####

# remove -3 (missing)

nac <- nac %>% mutate(across(contains("BAI"), ~na_if(., -3)))

# check histograms for outliers

nac %>% 
  select(contains("BAI")) %>% 
  select(contains("BL"), contains("W4")) %>% 
  multi.hist(global = F)

nac %>% 
  select(contains("BAI")) %>% 
  select(contains("W8"), contains("W12")) %>% 
  multi.hist(global = F)

nac %>% 
  select(contains("BAI")) %>% 
  select(contains("W16"), contains("W20")) %>% 
  multi.hist(global = F)

# missed items at BL, W4, W8, W20
# impute missing items

nac_imp3 <- nac %>% 
  mutate(across(contains("BAI"), haven::zap_labels))

# BL

bl_BAI <- nac_imp3 %>% select(BAI_BL_01:BAI_BL_21) %>% names()

nac_imp3 <- pmm_impute(vars = bl_BAI, max_missing = 10, data = nac_imp3)

# W4

W4_BAI <- nac_imp3 %>% select(BAI_W4_01:BAI_W4_21) %>% names()

nac_imp3 <- pmm_impute(vars = W4_BAI, max_missing = 10, data = nac_imp3)

# W8

W8_BAI <- nac_imp3 %>% select(BAI_W8_01:BAI_W8_21) %>% names()

nac_imp3 <- pmm_impute(vars = W8_BAI, max_missing = 10, data = nac_imp3)

# W12

W12_BAI <- nac_imp3 %>% select(BAI_W12_01:BAI_W12_21) %>% names()

nac_imp3 <- pmm_impute(vars = W12_BAI, max_missing = 10, data = nac_imp3)

# W16

W16_BAI <- nac_imp3 %>% select(BAI_W16_01:BAI_W16_21) %>% names()

nac_imp3 <- pmm_impute(vars = W16_BAI, max_missing = 10, data = nac_imp3)

# W20

W20_BAI <- nac_imp3 %>% select(BAI_W20_01:BAI_W20_21) %>% names()

nac_imp3 <- pmm_impute(vars = W20_BAI, max_missing = 10, data = nac_imp3)

## recalculate total scores

nac_imp3$BAI_BL_total2 <- rowSums(nac_imp3[,bl_BAI])

nac_imp3$BAI_W4_total2 <- rowSums(nac_imp3[,W4_BAI])

nac_imp3$BAI_W8_total2 <- rowSums(nac_imp3[,W8_BAI])

nac_imp3$BAI_W12_total2 <- rowSums(nac_imp3[,W12_BAI])

nac_imp3$BAI_W16_total2 <- rowSums(nac_imp3[,W16_BAI])

nac_imp3$BAI_W20_total2 <- rowSums(nac_imp3[,W20_BAI])

# raw means

nac_imp3 %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(BAI_BL_total2, BAI_W20_total2),
                   list(mean = ~mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

## long format 

bai_vars <- c("BAI_W4_total2", "BAI_W8_total2", "BAI_W12_total2", "BAI_W16_total2", "BAI_W20_total2")

nac_long_imp3 <- nac_imp3 %>% 
  pivot_longer(cols = all_of(bai_vars), 
               names_to = "time", 
               values_to = "BAI")

# check BAI variable for outliers

hist(nac_long_imp3$BAI, breaks = 30)

# time variable

nac_long_imp3$time <- factor(nac_long_imp3$time, levels = unique(nac_long_imp3$time))

nac_long_imp3$time <- as.numeric(nac_long_imp3$time) - 1

# model 

nac_long_imp3$zBAI_BL_total2 <- scale(nac_long_imp3$BAI_BL_total2)

BAI_m1 <- lmer(BAI ~ zBAI_BL_total2 + time*Group_allocation + zAge + Site +
                 Gender + zBL_standard_drinks_week + (1 + time | Participant_ID), 
               data = nac_long_imp3)

summary(BAI_m1)

emmeans(BAI_m1, "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

# Cohen's d

c(-1.78, -5.61, 2.04) / sd(nac_imp3$BAI_BL_total2, na.rm=T)


#### SDS ####

# calculate global impairment score

nac <- nac %>%
  rowwise() %>% 
  mutate(SDS_global_BL = sum(c_across(SDS_BL_work_school:SDS_BL_family_life_home)),
         SDS_global_W4 = sum(c_across(SDS_W4_work_school:SDS_W4_family_life_home)),
         SDS_global_W8 = sum(c_across(SDS_W8_work_school:SDS_W8_family_life_home)),
         SDS_global_W12 = sum(c_across(SDS_WK12_work_school:SDS_WK12_family_life_home)),
         SDS_global_W16 = sum(c_across(SDS_WK16_work_school:SDS_WK16_family_life_home)),
         SDS_global_W20 = sum(c_across(SDS_WK20_work_school:SDS_WK20_family_life_home))) %>% 
  ungroup()

# check histograms for outliers

nac %>% 
  select(starts_with("SDS_global")) %>% 
  multi.hist()

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(SDS_global_BL, SDS_global_W20),
                   list(mean = ~mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

## long format 

sds_vars <- nac %>% 
  select(starts_with("SDS_global")) %>% 
  select(-SDS_global_BL) %>% 
  names()

nac_long_sds <- nac %>% 
  pivot_longer(cols = all_of(sds_vars), 
               names_to = "time", 
               values_to = "SDS")

# time variable

nac_long_sds$time <- factor(nac_long_sds$time, levels = unique(nac_long_sds$time))

nac_long_sds$time <- as.numeric(nac_long_sds$time) - 1

# model 

nac_long_sds$zSDS_global_BL <- scale(nac_long_sds$SDS_global_BL)

SDS_m1 <- lmer(SDS ~ zSDS_global_BL*time + time*Group_allocation + zAge + Site +
                Gender + zBL_standard_drinks_week + (1 + time | Participant_ID), 
               data = nac_long_sds)

summary(SDS_m1)

emmeans(SDS_m1, "Group_allocation", 
        at = list(time = 4),
        weights = "proportional") %>% 
  pairs(reverse = T, infer = T)

# Cohen's d

c(-1.04, -4.05, 1.97) / sd(nac$SDS_global_BL, na.rm=T)

#### Quality of life ####

# Remove -3 

nac <- nac %>% 
  mutate(across(contains("WHOQOL"), ~ na_if(., -3)))

# check histograms for outliers

nac %>% 
  select(contains("WHOQOL")) %>% 
  select(ends_with("01")) %>% 
  multi.hist(global = F)

# raw proportions

nac %>% 
  filter(!is.na(WHOQOL_BL_01)) %>% 
  group_by(Group_allocation, WHOQOL_BL_01) %>% 
  tally() %>%
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

nac %>% 
  filter(!is.na(WHOQOL_WK20_01)) %>% 
  group_by(Group_allocation, WHOQOL_WK20_01) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(WHOQOL_BL_01, WHOQOL_WK20_01),
                   list(mean = ~mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

# data to long format 

nac_long3 <- nac %>% 
  pivot_longer(c(WHOQOL_WK12_01, WHOQOL_WK20_01),
               names_to = "time",
               values_to = "WHOQOL")

# center time variable at week 20

nac_long3$time <- as.numeric(factor(nac_long3$time,levels = unique(nac_long3$time))) -2

# WHOQOL variable is ordinal 

nac_long3$WHOQOL <- as.ordered(nac_long3$WHOQOL)

## Mixed effects ordinal regression model ##

# We fit the model using Bayesian proportional odds regression 
# from R package 'brms'. The first model uses weakly skeptical priors on the treatment effect.
# The normal(0,1) priors on treatment coefficients (group and group:time) give low (<5%)
# prior probability to odds ratios smaller than 0.17 or larger than 7.4. This shrinks very large ORs
# which are contextually very implausible towards zero. 
# We also perform a sensitivity analysis with wide (uninformative) priors to 
# demonstrate that results are robust to choice of prior

who_weakskep <- 
  brm(WHOQOL ~ zWHOQOL_BL_01 + time*Group_allocation + zAge + Site +
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long3,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,1), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,1), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, 
      cores = 4,
      file = "fits/whoqol_weakskep",
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      backend = "cmdstanr")


# posterior predictive check - does the model reproduce the data well?

pp_check(who_weakskep, 
         type = "bars_grouped", 
         group = "Group_allocation", 
         ndraws = 30) # good fit 

# results

summary(who_weakskep)

fixef(who_weakskep)

# Test hypothesis of NAC treatment improves QoL (coefficient > 0)

hyp <- hypothesis(who_weakskep, "Group_allocationNAC > 0")

hyp

# sensitivity analysis with less informative priors.
# normal(0, SD = 5) priors compatible with odds ratios
# from 0.007 to 148

who_sens <- 
  brm(WHOQOL ~ zWHOQOL_BL_01 + time*Group_allocation + zAge + Site +
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long3,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,5), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,5), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, 
      cores = 4,
      file = "fits/whoqol_sens",
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      backend = "cmdstanr")

summary(who_sens)

fixef(who_sens)

# hypothesis test

hyp <- hypothesis(who_sens, "Group_allocationNAC > 0")

hyp

#### CGI-S ####

# Remove -3 

nac <- nac %>% 
  mutate(across(contains("CGI"), ~ na_if(., -3)))

# check histograms for outliers

nac %>% 
  select(contains("CGI_S")) %>% 
  multi.hist(global = F)

# raw proportions

nac %>% 
  group_by(Group_allocation, CGI_S_BL) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

nac %>% 
  filter(!is.na(CGI_S_W20)) %>% 
  group_by(Group_allocation, CGI_S_W20) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(CGI_S_BL, CGI_S_W20),
                   list(mean = ~mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))

# data to long format 

nac_long4 <- nac %>% 
  pivot_longer(c(CGI_S_W4, CGI_S_W8, CGI_S_WK12, CGI_S_W16, CGI_S_W20),
               names_to = "time",
               values_to = "CGI_S")

hist(nac_long4$CGI_S)

# center time variable at week 20

nac_long4$time <- as.numeric(factor(nac_long4$time,levels = unique(nac_long4$time))) - 5

# WHOQOL variable is ordinal 

nac_long4$CGI_S <- as.ordered(nac_long4$CGI_S)

## Mixed effects ordinal regression model

# standardize continuous predictors to Z scores

nac_long4$zCGI_S_BL <- scale(nac_long4$CGI_S_BL)

# fit model

cgis_weakskep <- 
  brm(CGI_S ~ zCGI_S_BL + time*Group_allocation + zAge + Site +
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long4,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,1), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,1), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, cores = 4,
      file = "fits/cgis_weakskep",
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

# posterior predictive check

pp_check(cgis_weakskep, 
         type = "bars_grouped", 
         group = "Group_allocation", 
         ndraws = 30) # good fit

# results

summary(cgis_weakskep)
fixef(cgis_weakskep)

# Test hypothesis of NAC treatment benefit (coefficient < 0)

hyp <- hypothesis(cgis_weakskep, "Group_allocationNAC < 0")

hyp

# sensitivity model with weak priors 

cgis_sens <- 
  brm(CGI_S ~ zCGI_S_BL + time*Group_allocation + zAge + Site +
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long4,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,5), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,5), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, cores = 4,
      file = "fits/cgis_sens",
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

fixef(cgis_sens)

# hypothesis test

hyp <- hypothesis(cgis_sens, "Group_allocationNAC < 0")

hyp


#### PGI-S ####

# remove missing data 

nac <- nac %>% 
  mutate(across(contains("PGI"), ~ na_if(., -3)))

nac <- nac %>% 
  mutate(PGI_S_W20 = as.numeric(PGI_S_W20)) %>% 
  mutate(PGI_S_W20 = if_else(PGI_S_W20 == 0, 3, PGI_S_W20))


# check histograms for outliers

nac %>% 
  select(contains("PGI_S")) %>% 
  multi.hist(global = F)


# raw proportions

nac %>% 
  filter(!is.na(PGI_S_BL)) %>% 
  group_by(Group_allocation, PGI_S_BL) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

nac %>% 
  filter(!is.na(PGI_S_W20)) %>% 
  group_by(Group_allocation, PGI_S_W20) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(across(c(PGI_S_BL, PGI_S_W20),
                   list(mean = ~mean(.,na.rm=T),
                        sd = ~ sd(., na.rm=T))))


# data to long format 

nac_long5 <- nac %>% 
  pivot_longer(c(PGI_S_W4, PGI_S_W8, PGI_S_WK12, PGI_S_W16, PGI_S_W20),
               names_to = "time",
               values_to = "PGI_S")

hist(nac_long5$PGI_S)

# center time variable at week 20

nac_long5$time <- as.numeric(factor(nac_long5$time,levels = unique(nac_long5$time))) - 5

# WHOQOL variable is ordinal 

nac_long5$PGI_S <- as.ordered(nac_long5$PGI_S)

## Mixed effects ordinal regression model

# standardize continuous predictors to Z scores

nac_long5$zPGI_S_BL <- scale(nac_long5$PGI_S_BL)

# fit model

pgis_weakskep <- brm(PGI_S ~ zPGI_S_BL + time*Group_allocation + zAge + Site +
                       Gender + zBL_standard_drinks_week + (1 | Participant_ID),
                     data = nac_long5,
                     family = cumulative(link = "probit", threshold = "flexible"),
                     prior = c(prior(normal(0,1), class = "b", coef = "Group_allocationNAC"),
                               prior(normal(0,1), class = "b", coef = "time:Group_allocationNAC"),
                               prior(normal(0,10), class = "b"),
                               prior(normal(0,10), class = "Intercept"),
                               prior(normal(0,2), class = "sd")),
                     iter = 5000,
                     warmup = 1500, 
                     chains = 4, cores = 4,
                     file = "fits/pgis_weakskep",
                     control = list(adapt_delta = 0.95, max_treedepth = 15),
                     backend = "cmdstanr")

# posterior predictive check

pp_check(pgis_weakskep, 
         type = "bars_grouped", 
         group = "Group_allocation", 
         ndraws = 30) # good fit

# results

summary(pgis_weakskep)

fixef(pgis_weakskep)

# Test hypothesis of NAC treatment benefit (coefficient < 0)

hyp <- hypothesis(pgis_weakskep, "Group_allocationNAC < 0")

hyp

# sensitivity model with weak priors 

pgis_sens <- brm(PGI_S ~ zPGI_S_BL + time*Group_allocation + zAge + Site +
                    Gender + zBL_standard_drinks_week + (1 | Participant_ID),
                     data = nac_long5,
                     family = cumulative(link = "probit", threshold = "flexible"),
                     prior = c(prior(normal(0,5), class = "b", coef = "Group_allocationNAC"),
                               prior(normal(0,5), class = "b", coef = "time:Group_allocationNAC"),
                               prior(normal(0,10), class = "b"),
                               prior(normal(0,10), class = "Intercept"),
                               prior(normal(0,2), class = "sd")),
                     iter = 5000,
                     warmup = 1500, 
                     chains = 4, cores = 4,
                     file = "fits/pgis_sens",
                     control = list(adapt_delta = 0.95, max_treedepth = 15),
                     backend = "cmdstanr")

# results

summary(pgis_sens)

fixef(pgis_sens)

# Test hypothesis of NAC treatment benefit (coefficient < 0)

hyp <- hypothesis(pgis_sens, "Group_allocationNAC < 0")

hyp

#### CGI-I ####

# remove -3 (SPSS missing)

nac <- nac %>% 
  mutate(across(contains("CGI"), ~ na_if(., -3)))

# check histograms for outliers

nac %>% 
  select(contains("CGI_I")) %>% 
  multi.hist(global = F)

# raw proportions

nac %>% 
  filter(!is.na(CGI_I_W20)) %>% 
  group_by(Group_allocation, CGI_I_W20) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(mean = mean(CGI_I_W20, na.rm=T),
            sd = sd(CGI_I_W20, na.rm=T))

# data to long format 

nac_long6 <- nac %>% 
  pivot_longer(c(CGI_I_W4, CGI_I_W8, CGI_I_WK12, CGI_I_W16, CGI_I_W20),
               names_to = "time",
               values_to = "CGI_I")

# center time variable at week 20

nac_long6$time <- as.numeric(factor(nac_long6$time,levels = unique(nac_long6$time))) - 5

# WHOQOL variable is ordinal 

nac_long6$CGI_I <- as.ordered(nac_long6$CGI_I)

## Mixed effects ordinal regression model

cgii_weakskep <- 
  brm(CGI_I ~ time*Group_allocation + zAge + Site +
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long6,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,1), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,1), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, cores = 4,
      file = "fits/cgii_weakskep",
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

# posterior predictive check

pp_check(cgii_weakskep, 
         type = "bars_grouped",
         group = "Group_allocation",
         ndraws = 30) # poor fit

# results

summary(cgii_weakskep)

fixef(cgii_weakskep)

# plot

conditional_effects(cgii_weakskep, categorical = T, effects = "Group_allocation")

# Test hypothesis of NAC treatment benefit (coefficient < 0)

hyp <- hypothesis(cgii_weakskep, "Group_allocationNAC < 0")

hyp

# sensitivity model with weak priors 

cgii_sens <- brm(CGI_I ~ time*Group_allocation + zAge + Site + 
                    Gender + zBL_standard_drinks_week + (1 | Participant_ID),
                  data = nac_long6,
                  family = cumulative(link = "probit", threshold = "flexible"),
                  prior = c(prior(normal(0,5), class = "b", coef = "Group_allocationNAC"),
                            prior(normal(0,5), class = "b", coef = "time:Group_allocationNAC"),
                            prior(normal(0,10), class = "b"),
                            prior(normal(0,10), class = "Intercept"),
                            prior(normal(0,2), class = "sd")),
                  iter = 5000,
                  warmup = 1500, 
                  chains = 4, cores = 4,
                  file = "fits/cgii_sens",
                  control = list(adapt_delta = 0.95, max_treedepth = 15),
                  backend = "cmdstanr")

fixef(cgii_sens)

hyp <- hypothesis(cgii_sens, "Group_allocationNAC < 0")

hyp

#### PGI-I ####

# remove -3 (SPSS missing)

nac <- nac %>% 
  mutate(across(contains("PGI"), ~ na_if(., -3)))

# check histograms for outliers

nac %>% 
  select(contains("PGI_I")) %>% 
  multi.hist(global = F)

# raw proportions

nac %>% 
  filter(!is.na(PGI_I_W20)) %>% 
  group_by(Group_allocation, PGI_I_W20) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(total = sum(n))

# raw means

nac %>% 
  group_by(Group_allocation) %>% 
  summarise(mean = mean(PGI_I_W20, na.rm=T),
            sd = sd(PGI_I_W20, na.rm=T))

# data to long format 

nac_long7 <- nac %>% 
  pivot_longer(c(PGI_I_W4, PGI_I_W8, PGI_I_WK12, PGI_I_W16, PGI_I_W20),
               names_to = "time",
               values_to = "PGI_I")

# center time variable at week 20

nac_long7$time <- as.numeric(factor(nac_long7$time,levels = unique(nac_long7$time))) - 5

# WHOQOL variable is ordinal 

nac_long7$PGI_I <- as.ordered(nac_long7$PGI_I)

## Mixed effects ordinal regression model

PGIi_weakskep <- 
  brm(PGI_I ~ time*Group_allocation + zAge + 
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long7,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,1), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,1), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, cores = 4,
      file = "fits/PGIi_weakskep",
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

# posterior predictive check

pp_check(PGIi_weakskep, 
         type = "bars_grouped",
         group = "Group_allocation",
         ndraws = 30) 

# results

summary(PGIi_weakskep)

fixef(PGIi_weakskep)

# Test hypothesis of NAC treatment benefit (coefficient < 0)

hyp <- hypothesis(PGIi_weakskep, "Group_allocationNAC < 0")

hyp

# sensitivity model with weak priors 

PGIi_sens <- 
  brm(PGI_I ~ time*Group_allocation + zAge + 
        Gender + zBL_standard_drinks_week + (1 | Participant_ID),
      data = nac_long7,
      family = cumulative(link = "probit", threshold = "flexible"),
      prior = c(prior(normal(0,5), class = "b", coef = "Group_allocationNAC"),
                prior(normal(0,5), class = "b", coef = "time:Group_allocationNAC"),
                prior(normal(0,10), class = "b"),
                prior(normal(0,10), class = "Intercept"),
                prior(normal(0,2), class = "sd")),
      iter = 5000,
      warmup = 1500, 
      chains = 4, cores = 4,
      file = "fits/PGIi_sens",
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

fixef(PGIi_sens)

hyp <- hypothesis(PGIi_sens, "Group_allocationNAC < 0")

hyp


#### Blood pressure ####

### systolic ###

# check histograms

nac %>% 
  filter(!is.na(YBOCS_W4_total2)) %>% 
  group_by(Group_allocation) %>% 
  tally()

nac %>% 
  select(starts_with("systolic")) %>% 
  multi.hist(global = F)

# raw means

nac %>% 
  summarise(across(starts_with("systolic"), list(mean = ~ mean(., na.rm = T),
                                                 sd = ~ sd(., na.rm=T))))

## long format data ##

bp_long <- nac %>% 
  pivot_longer(c(Systolic_BP_W4, Systolic_BP_W8,
                 Systolic_BP_W12, Systolic_BP_W16, Systolic_BP_W20),
               names_to = "time",
               values_to = "systolic",
               names_prefix = "Systolic_BP_")

bp_long$time <- as.numeric(factor(bp_long$time, levels = unique(bp_long$time))) - 1

## model

m_bp <- lme(systolic ~ Systolic_BP_baseline + time * Group_allocation + Gender + Age, 
            random = ~ 1 | Participant_ID,
            data = bp_long,
            na.action = "na.omit")

summary(m_bp)

## Plot

preds <- ggemmeans(m_bp, 
                   terms = c("time","Group_allocation"))

preds

# add baseline data to predictions


bl <- data.frame(x = c(-1,-1), 
                 predicted = c(120.55, 120.55), # baseline mean
                 std.error = c(NA,NA),
                 conf.low = c(NA,NA),
                 conf.high = c(NA,NA),
                 group = c("Placebo","NAC"))

preds <- rbind(bl, preds)

# plot 

ggplot(preds, aes(x = x, y = predicted, colour = group)) +
  geom_line() + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0.1),
                  fatten = 2) +
  scale_x_continuous(labels = c("BL","W4","W8","W12","W16","W20")) +
  xlab("Visit") + ylab("Systolic BP, 95% CI") + 
  ylim(c(110,130)) +
  ggthemes::scale_colour_colorblind() +
  ggthemes::theme_stata()

# save plot to png

ggsave("NAC systolic.png", device = "png", path = "C:/Users/lachy/OneDrive - The University of Melbourne/Unimelb/NAC placebo/Plots")



### Diastolic ###

nac <- nac %>% 
  mutate(across(contains("Diastolic"), ~ na_if(., 3)))

nac %>% 
  select(starts_with("diastolic")) %>% 
  multi.hist(global = F)

# raw means

nac %>% 
  summarise(across(starts_with("diastolic"), list(mean = ~ mean(., na.rm = T),
                                                 sd = ~ sd(., na.rm=T))))

## long format data ##

bp_long <- nac %>% 
  pivot_longer(c(Diastolic_BP_W4, Diastolic_BP_W8,
                 Diastolic_BP_W12, Diastolic_BP_W16, Diastolic_BP_W20),
               names_to = "time",
               values_to = "diastolic",
               names_prefix = "Diastolic_BP_")

bp_long$time <- as.numeric(factor(bp_long$time, levels = unique(bp_long$time))) - 1

## model

m_bp <- lme(diastolic ~ Diastolic_BP_baseline + time * Group_allocation + Gender + Age, 
            random = ~ 1 | Participant_ID,
            data = bp_long,
            na.action = "na.omit")

summary(m_bp)

## Plot

preds <- ggemmeans(m_bp, 
                   terms = c("time","Group_allocation"))

preds

# add baseline data to predictions

bl <- data.frame(x = c(-1,-1), 
                 predicted = c(79.81, 79.81), # baseline mean
                 std.error = c(NA,NA),
                 conf.low = c(NA,NA),
                 conf.high = c(NA,NA),
                 group = c("Placebo","NAC"))

preds <- rbind(bl, preds)

# plot 

ggplot(preds, aes(x = x, y = predicted, colour = group)) +
  geom_line() + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0.1),
                  fatten = 2) +
  scale_x_continuous(labels = c("BL","W4","W8","W12","W16","W20")) +
  xlab("Visit") + ylab("Diastolic BP, 95% CI") + 
  ylim(c(70,85)) +
  ggthemes::scale_colour_colorblind() +
  ggthemes::theme_stata()

# save plot to png

ggsave("NAC diastolic.png", device = "png", path = "C:/Users/lachy/OneDrive - The University of Melbourne/Unimelb/NAC placebo/Plots")

  