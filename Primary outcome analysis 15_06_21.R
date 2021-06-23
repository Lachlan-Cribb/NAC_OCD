library(haven)
library(tidyverse)
library(nlme)
library(sjPlot)
library(sjlabelled)
library(ggeffects)
library(extrafont)
library(brms)
library(emmeans)
library(mice)
loadfonts(device = "win")

nac <- read_sav("C:\\Users\\lachy\\OneDrive - The University of Melbourne\\Unimelb\\NAC placebo\\datafiles\\COGNAC data file (2) MERGED 28-Apr-2021.sav")


# remove -9 and -999 SPSS missing data indicators 

nac <- nac %>% 
  mutate(across(where(is.numeric), ~ na_if(., -9))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -99)))

# drop unneeded variable with problematic formatting

nac <- nac[,-890]

#### Primary outcome ####

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

nac$Group_allocation <- sjlabelled::as_label(nac$Group_allocation)

nac$Group_allocation <- relevel(nac$Group_allocation, "B") # placebo reference group


### data to long format ###


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



#pattern mixture model 

pm1 <- lme(YBOCS ~ YBOCS_BL2_total2 + time*Group_allocation*drop_any,
              random = ~ time | Participant_ID,
              na.action = na.omit, data = nac_long)

summary(pm1)

# average over dropout pattern
em <- emmeans(pm1,
              pairwise ~ Group_allocation | time, 
              at = list(time = 4),
              weights = "proportional")

em


#### Per protocol analysis ####



#### DOCS ####

# check for missing data (repeat for each visit) #

print(nac %>% 
        rowwise() %>% 
        mutate(total_na = sum(is.na(c_across(DOCS_W20_contamination_01:DOCS_W20_total)))) %>% 
        dplyr::select(total_na) %>% 
        arrange(desc(total_na)), n = 100)

## Impute missing DOCS items ##

# remove labels for MICE imputation

nac_imp <- nac %>% 
  mutate(across(contains("DOCS"), haven::zap_labels))

# imputation function

docs_impute <- function(vars,max_missing){
  
  skip <- which(rowSums(is.na(nac[,vars]))>max_missing) # maximum missing for imputation
  
  imp <- mice(data = nac_imp[-skip,vars], m = 1, method = "pmm") # impute
  
  imp <- complete(imp) # save as dataframe 
  
  nac_imp[-skip,vars] <- imp # replace data with imputed data 
  
  nac_imp
}


# BL imputation - skipping participants with more than 50% missing #

bl_docs <- names(nac %>% 
                   select(DOCS_BL_contamination_01:DOCS_BL_symmetry_05) %>% 
                   select(-DOCS_BL_contamination_total, -DOCS_BL_harm_total, -DOCS_BL_taboo_total))


nac_imp <- docs_impute(vars = bl_docs, max_missing = 12)


# week 4

W4_docs <- names(nac %>% 
                   select(DOCS_W4_contamination_01:DOCS_W4_symmetry_05) %>% 
                   select(-DOCS_W4_contamination_total, -DOCS_W4_harm_total, -DOCS_W4_taboo_total))

nac_imp <- docs_impute(vars = W4_docs, max_missing = 12)


# week 8 

W8_docs <- names(nac %>% 
                   select(DOCS_W8_contamination_01:DOCS_W8_symmetry_05) %>% 
                   select(-DOCS_W8_contamination_total, -DOCS_W8_harm_total, -DOCS_W8_taboo_total))

nac_imp <- docs_impute(vars = W8_docs, max_missing = 12)


# week 12

W12_docs <- names(nac %>% 
                   select(DOCS_W12_contamination_01:DOCS_W12_symmetry_05) %>% 
                   select(-DOCS_W12_contamination_total, -DOCS_W12_harm_total, -DOCS_W12_taboo_total))

nac_imp <- docs_impute(vars = W12_docs, max_missing = 12)

# Week 16 

W16_docs <- nac %>% 
  select(DOCS_W16_contamination_01:DOCS_W16_symmetry_05) %>% 
  select(-DOCS_W16_contamination_total, -DOCS_W16_harm_total, -DOCS_W16_taboo_total) %>% 
  names()

nac_imp <- docs_impute(vars = W16_docs, max_missing = 12)

# week 20

W20_docs <- nac %>% 
  select(DOCS_W20_contamination_01:DOCS_W20_symmetry_05) %>% 
  select(-DOCS_W20_contamination_total, -DOCS_W20_harm_total, -DOCS_W20_taboo_total) %>% 
  names()

nac_imp <- docs_impute(vars = W20_docs, max_missing = 12)


## recalculate total scores

nac_imp$DOCS_BL_total2 <- rowSums(nac_imp[,bl_docs])

nac_imp$DOCS_W4_total2 <- rowSums(nac_imp[,W4_docs])

nac_imp$DOCS_W8_total2 <- rowSums(nac_imp[,W8_docs])

nac_imp$DOCS_W12_total2 <- rowSums(nac_imp[,W12_docs])

nac_imp$DOCS_W16_total2 <- rowSums(nac_imp[,W16_docs])

nac_imp$DOCS_W20_total2 <- rowSums(nac_imp[,W20_docs])


## DOCS mixed model 

# data to wide format

docs_vars <- c("DOCS_W4_total2", "DOCS_W8_total2", "DOCS_W12_total2", "DOCS_W16_total2", "DOCS_W20_total2")

nac_long_imp <- nac_imp %>% 
  pivot_longer(cols = all_of(docs_vars), names_to = "Visit", values_to = "DOCS")

# continuous time variable

nac_long_imp$Visit <- factor(nac_long_imp$Visit, levels = unique(nac_long_imp$Visit))

nac_long_imp$Visit <- as.numeric(nac_long_imp$Visit) - 1


# model 

docs_m1 <- lme(DOCS ~ Visit * Group_allocation,
               random = ~ 1 | Participant_ID,
               na.action = na.omit, data = nac_long_imp)

summary(docs_m1)

em <- emmeans(docs_m1,
              pairwise ~ Group_allocation | Visit, 
              at = list(Visit = 4))

em

