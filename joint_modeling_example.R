### create dataset for joint modeling

rm(list = ls())

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(bayesCT)
library(survival)
library(survminer)
library("JM")

#####----- biomarker dataset generation -----#####

### generate 1000 patients and biomarker values for them,
### which has linear relationship with time: biomarker = intercept + slope * t.
### Parameters of the relationship for each patients will be sampled from normal distribution

### number of patients and times for which biomarkers values will be modeled
n = 1000
t = 1:35
t_last <- max(t)

### population means and SDs for parameters of biomarker model
slope_m = 0.2
slope_s = 0.3
intercept_m = 0
intercept_s = 1
error_s = 0.3

### dataset with biomarker profiles; biom_err - observed biomarker value

set.seed(1)
dataset_biom <- tibble(id = 1:n, intercept = rnorm(n, intercept_m, intercept_s), 
                      slope = rnorm(n, slope_m, slope_s)) %>% 
    left_join(tibble(time = t), by = character()) %>% 
    mutate(biom = intercept + slope * time,
           biom_obs = biom + rnorm(n * max(t), 0, error_s)) %>% 
    arrange(id, time)

dataset_biom_BL <- dataset_biom %>% 
    dplyr::select(id, intercept) %>% 
    unique() %>% 
    rename(biom = intercept) %>% 
    mutate(time = 0, biom_obs = biom + rnorm(n, 0, error_s))

### adding baseline value (at t = 0)

dataset_biom_mod <- dataset_biom %>% 
    dplyr::select(id, time, biom, biom_obs) %>% 
    bind_rows(dataset_biom_BL) %>% 
    arrange(id, time)

### summary stats for biomarker dynamics

dataset_biom_m <- dataset_biom_mod %>% 
    group_by(time) %>% 
    summarise_at(vars(biom_obs), list(mean = ~mean(., na.rm = T),
                                     sd = ~sd(., na.rm = T),
                                     se = ~sd(., na.rm = T)/n())) %>% 
    ungroup()

### plot observed biomarker profiles if no drop-out/deaths occured

ggplot(data = dataset_biom_mod) +
    geom_line(aes(x = time, y = biom_obs, group = id), size = 0.1, alpha = .3) +
    geom_point(data = dataset_biom_m, aes(x = time, y = mean), size = 1, col = "yellow") +
    geom_line(data = dataset_biom_m, aes(x = time, y = mean), size = 0.5, col = "yellow")


#####----- time to death generation -----#####
### time-to-event generated with function pw_exp_sim, which generate time-to-event,
### based on vector of hazards and vector of times. hazard value i from hazard vector corresponds
### to time period [i-1, i) from times vector. vector of hazards is generated using the expsression:
### h_t = h_0 * exp(a_m * m_t), where h_0 - default hazard function, a_m - slope of relationship h_t and m_t,
### m_t biomarker's value at time t.

### parameters of hazard function
def_hazard = 0.05
alpha_m = 0.3

### function to generate times-to-event for a patient
tte_gen <- function(h_0 = def_hazard, a_m = alpha_m, a, b, t_set = seq(0, t_last, 0.1)){
    
    m_t <- a + b * t_set
    h_t <- h_0 * exp(alpha_m * m_t)
    time <- pw_exp_sim(h_t, n = 1, cutpoint = t_set[1:(length(t_set)-1)])
    return(time)
}

### generate times of event for patients from dataset_biom_mod

dataset_tte <- dataset_biom %>% 
    dplyr::select(id, intercept, slope) %>% 
    unique()

for (i in dataset_tte$id){
    
    tte_i <- tte_gen(a = dataset_tte[dataset_tte$id == i,]$intercept, 
                     b = dataset_tte[dataset_tte$id == i,]$slope)
    
    dataset_tte[dataset_tte$id == i,4] <- tte_i$time
}

names(dataset_tte)[4] <- "tte"

### plotting survival function for generated times of event versus
### survival function from simple exponential distribution with hazard = def_hazard

dataset_km_jm <- dataset_tte %>% 
    dplyr::select(id, tte) %>% 
    mutate(event = if_else(tte > 25, 0, 1), 
           tte = if_else(tte > 25, 25, tte),
           group = "jm")

dataset_km_exp <- tibble(id = 1:1000, tte = rexp(1000, def_hazard), event = 1) %>% 
    mutate(event = if_else(tte > 25, 0, 1), 
           tte = if_else(tte > 25, 25, tte),
           group = "exp")

dataset_km <- bind_rows(dataset_km_jm, dataset_km_exp)

fit <- survfit(Surv(tte, event) ~ group,
               data = dataset_km)

ggsurvplot(fit, data = dataset_km, risk.table = TRUE)

### plotting observed biomarker's profiles as some subject withdrew due to event occured

dataset_biom_mod_e <- dataset_biom_mod %>% 
    left_join(dataset_km_jm %>% dplyr::select(id, tte, event), by = "id") %>% 
    filter(time <= tte) %>% 
    mutate(event_char = if_else(event == 1, "event occured", "no event"))

dataset_biom_mod_e_m <- dataset_biom_mod_e %>% 
    group_by(time) %>% 
    summarise_at(vars(biom_obs), list(mean = ~mean(., na.rm = T),
                                      sd = ~sd(., na.rm = T),
                                      se = ~sd(., na.rm = T)/n())) %>% 
    ungroup()

ggplot(data = dataset_biom_mod_e) +
    geom_line(aes(x = time, y = biom_obs, group = id, col = event_char), size = 0.1, alpha = .3) +
    geom_point(data = dataset_biom_mod_e_m, aes(x = time, y = mean), size = 1, col = "yellow") +
    geom_line(data = dataset_biom_mod_e_m, aes(x = time, y = mean), size = 0.5, col = "yellow") +
    scale_color_manual(values = c("red", "blue"))
    labs(color = "Event")

#####----- Modeling -----#####

### linear model for biomarker with random effects for intercept and slope

### model for dataset without withdrawn subjects
    
lmemodel_full_dt <- lme(biom_obs ~ time, #formula for fixed effects
                        random = ~ time |id, #formula for random effects
                        data = dataset_biom_mod)
    
### model for dataset with withdrawn subjects (observed dataset)

lmemodel_e_dt <- lme(biom_obs ~ time, #formula for fixed effects
                        random = ~ time |id, #formula for random effects
                        data = dataset_biom_mod_e)

summary(lmemodel_full_dt)
summary(lmemodel_e_dt)

### Cox submodel for baseline biomarker's value

dataset_km_jm_mod <- dataset_km_jm %>% 
    left_join(dataset_biom_BL %>% dplyr::select(id, biom_obs), by = "id")
    
cox <- coxph(Surv(tte,event) ~ biom_obs, data = dataset_km_jm_mod, x = TRUE)
summary(cox)

### joint model

joint_model <- jointModel(lmemodel_e_dt, cox, timeVar = "time",
                          method = "weibull-PH-aGH" # a time-dependent relative risk model
)

### summaries for joint model show alpha_m parameter (Assoct) closer to the value used for dataset generation than
### from Cox model (coef)

summary(joint_model)
