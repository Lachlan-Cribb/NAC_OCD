
# Load packages

library(tidyverse)
library(ggthemes)
library(bayesplot)
library(cowplot)

# read data

ae <- read.csv("PATH")

# format data

ae <- ae %>% 
  rename('Adverse_Event' = 'X') %>% 
  filter(!Adverse_Event == "")

ae <- ae %>% mutate(Adverse_Event = fct_recode(Adverse_Event,
                  'Ischaemic Colitis (Severe)' = 'Iscaemic Colitis (Severe)'))

ae$Adverse_Event <- factor(ae$Adverse_Event, levels = unique(ae$Adverse_Event))


# Calculate risk reduction (risk in NAC - risk in placebo)

ae <- ae %>% 
  mutate(total = A + B) %>% 
  filter(total > 2 | str_detect(Adverse_Event,"Severe")) %>% 
  rowwise() %>% 
  mutate(across(c("A","B"), ~ifelse(. < 1, 0.5, .), .names = "{.col}_2")) %>% 
  mutate(pA = A / 44, pB = B / 45) %>%
  mutate(pA2 = A_2 / 44, pB2 = B_2 / 45) %>% 
  mutate(RR = pA2 - pB2) %>% 
  mutate(lower = RR - 1.96 * sqrt(((pA2 * (1 - pA2)) / 44) + (((pB2 * (1 - pB2)) / 45)))) %>% 
  mutate(upper = RR + 1.96 * sqrt(((pA * (1 - pA2)) / 44) + (((pB2 * (1 - pB2)) / 45)))) %>% 
  select(Adverse_Event, A, B, pA, pB, A_2, B_2, RR, lower, upper) %>% 
  arrange(RR)



## plot

p1 <- ae %>% 
  ggplot(aes(x = RR, y = Adverse_Event)) + 
  geom_vline(xintercept = 0) + 
  geom_pointrange(aes(xmin = lower,
                     xmax = upper)) +
  labs(x = "Risk difference with 95% CI", y = "") +
  bayesplot_theme_set()
  
p1


## raw proportions in each group

p2 <- ae %>% pivot_longer(cols = c("pA","pB"), values_to = "prop",
                    names_to = "Group") %>% 
  ggplot(aes(x = prop, y = Adverse_Event, colour = Group, shape = Group)) +
  geom_point(size = 2) +
  bayesplot_theme_set() +
  labs(x = "Proportion", y = "") +
  theme(axis.text.y = element_blank()) +
  scale_colour_colorblind(labels = c("NAC","Placebo")) +
  scale_shape_discrete(labels = c("NAC","Placebo")) +
  geom_hline(aes(yintercept = Adverse_Event), size = 0.3)
  
p2


# Combine plots

cowplot::plot_grid(p1, p2, rel_widths = c(1.25,1))


# Save output 

ggsave("ae plot.png", device = "png", dpi = 450)


