#####################################################################
##                        Program set up                           ##
#####################################################################
options(mc.cores = parallel::detectCores())

nc <- switch(tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_")), 
             "true" =, "warn" = 2, 
             parallel::detectCores())
options(mc.cores = nc)

setwd("D:\\MPH-Project")
memory.limit(size=40000) # set memory

#####################################################################
##          Redistribute the outcome to increase                   ##
##          treatment effects in younger population                ##  
#####################################################################

## checking the missing data
full_ipd <- plaque_psoriasis_ipd %>%
  mutate(   # Check complete cases for covariates of interest
    complete = complete.cases(durnpso, prevsys, bsa, weight, psa, age, male))

#
s.comp <- sum(!full_ipd$complete)
mean.comp <- mean(!full_ipd$complete)

# remove the missing data
full_ipd <- filter(full_ipd, complete)

full_agd <- plaque_psoriasis_agd

## initial count by studyc, age group, pasi75 group
full.crs <- full_ipd %>%
  mutate(sub_age = ifelse(age>=47, "gt or eq 47", "lt 47"),
         cat_pas = ifelse(pasi75 == 1, "yes", "no")) %>%
  count(studyc, trtc, cat_pas, sub_age)

## histogram by studyc before redistribution
ggplot(full.crs) +
  geom_bar(aes(x=sub_age, y=n, fill = cat_pas), stat = "identity") +
  facet_wrap(vars(studyc, trtc)) +
  labs(title = "Distribution before simulation by trials, studies")+
  guides(fill=guide_legend(title="PASI75")) +
  xlab("Age group") +
  ylab("Frequency")

## nest data
full_ipd.muta <- full_ipd %>%
  group_by(studyc) %>%
  nest()

## create model and then predict the f(pasi75)
set.seed(1234)

full_ipd.muta$data <- map(full_ipd.muta$data, function(df){
  a <- glm(pasi75 ~ age*trtc + bmi*trtc + male*trtc + bsa*trtc + weight*trtc + durnpso*trtc, family = "binomial", data = df, na.action = "na.exclude")
  df %>%
    mutate(lp = predict(a, type = "link"))
})

# Ensures reproducible
full_ipd.muta$data <- map(full_ipd.muta$data, function(df){
  df %>%
    mutate(lp2 = if_else(age < 47 & trtc %in% c("PBO"),  lp - 0.2, lp),
           lp2 = if_else(age < 47 & !trtc %in% c("PBO"), lp2 + 0.2, lp2),
           risk = plogis(lp2),
           ## Next line is stochastic
           pasi75_sim = rbinom(n = length(lp2), size = 1, prob = risk))
})

full_ipd <- full_ipd.muta %>%
  unnest(data) %>%
  filter(!is.na(pasi75_sim))  ## 6 records with missing PASI75 deleted from data

full.crs2 <- full_ipd %>%
  mutate(cat_pas = ifelse(pasi75_sim == 1, "yes", "no"),
         sub_age = ifelse(age>=47, "gt or eq 47", "lt 47")) %>%
  count(studyc, trtc, cat_pas, sub_age)

## histogram by studyc after redistribution
ggplot(full.crs2) +
  geom_bar(aes(x=sub_age, y=n, fill = cat_pas), stat = "identity") +
  facet_wrap(vars(studyc, trtc)) +
  labs(title = "Distribution after simulation by trials, studies") +
  guides(fill=guide_legend(title="PASI75"))


#####################################################################
##                          Data preparation                       ##
#####################################################################

## Data transformation in IPD
sub_ipd <- full_ipd %>% 
  mutate(# Variable transformations
    bsa = bsa / 100,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    weight = weight / 10,
    durnpso = durnpso / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6, 7) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker"),
    sub_age = ifelse(age>=47, "gt or eq 47", "lt 47"),
    sub_sex = ifelse(male == "TRUE", "male", "female"),
    male = ifelse(male == "TRUE", 1, 0))

## subgroup result generated from IPD
## five covariate durnpso, prevsys, bsa, weight, psa, (age/male)
## include studyc, trt, trtclass, trtn, trtc_long into by vars

sub_eff_func <- function(group_vars, name_) {
  
  sub_ipd %>% 
    group_by(!!!group_vars) %>%
    summarise(age_mean = mean(age),
              age_sd = sd(age),
              bmi_mean = mean(bmi, na.rm = TRUE),
              bmi_sd = sd(bmi, na.rm = TRUE),
              weight_mean = mean(bmi, na.rm = TRUE),
              weight_sd = sd(bmi, na.rm = TRUE),
              durnpso_mean = mean(durnpso),
              durnpso_sd = sd(durnpso),
              pasi75_n = n(),
              pasi75_r = sum(pasi75),
              pasi90_n = n(),
              pasi90_r = sum(pasi90),
              pasi100_n = n(),
              pasi100_r = sum(pasi100), 
              psa = mean(psa), 
              prevsys = mean(prevsys), 
              male = mean(male),
              bsa_mean = mean(bsa), 
              bsa_sd = sd(bsa)) %>%
    mutate(studyc = paste(studyc, name_)) 
  
}

## convert IPD to AgD
ipd_agd <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long), "IPD")

## aggregrate data age-specific subgroup
agd_sub_age <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long, sub_age), "SUB")

## aggregrate data sex-specific subgroup
agd_sub_sex <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long, sub_sex), "SUB")

## aggregrate data factorial combination of age and sex subgroup
agd_sub_agesex <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long, sub_age, sub_sex), "SUB")


## original AgD mutation
sub_agd <- full_agd %>% 
  mutate(
    # Variable transformations
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    prevsys = prevsys / 100,
    psa = psa / 100,
    male = male / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6, 7) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker")
  )



## combine AGD with study level/treatment level subgroup
sub_agd1 <- bind_rows(sub_agd,
                      agd_sub_age, 
                      ipd_agd)
## sex
# sub_agd2 <- bind_rows(sub_agd,
#                       agd_sub_sex, 
#                       ipd_agd)
# 
# age and sex
# sub_agd3 <- bind_rows(sub_agd,
#                       agd_sub_agesex, 
#                       ipd_agd)

#####################################################################
##        Generate a self-defined function to run the model        ##
##        under different scenarios preparation                    ##
#####################################################################

#####################################################################
##             specification of function
##
## define a function to run the multiple ML-NMRs
## return a list contains nmr model
## agd.dt: the AgD
## col_name: default value is studyc, which helps to select 
## corresponding IPD and AgD
## value1: filter criteria in IPD
## value2: filter criteria in AgD
#####################################################################

func.mult <- function(agd.dt, col_name, value1, value2, sub_var) {
  
  filter_criteria1 <- interp(~y %in% x, .values=list(y = as.name(col_name), x = value1))
  
  filter_criteria2 <- interp(~y %in% x, .values=list(y = as.name(col_name), x = value2))
  
  net_org <- combine_network(
    set_ipd(sub_ipd %>% filter_(filter_criteria1),
            # set_ipd(sub_ipd,
            study = studyc, 
            trt = trtc, 
            r = pasi75_sim,
            trt_class = trtclass, 
            trt_ref = "PBO"),
    
    set_agd_arm(agd.dt %>% filter_(filter_criteria2), 
                study = studyc, 
                trt = trtc, 
                r = pasi75_r, 
                n = pasi75_n,
                trt_class = trtclass,
                trt_ref = "PBO")
  )
  
  # net_org
  # 
  netplot <- plot(net_org, weight_nodes = T, weight_edges = T, show_trt_class = T) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  
  
  ## Add integration points of to the AgD studies in the network
  if (sub_var == "age") {
    net_org <- add_integration(net_org,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               age = distr(qgamma, mean = age_mean, sd = age_sd),
                               n_int = 1000)
    ## fixed effect NL-NMR
    FE_net_org <- nma(net_org, 
                      trt_effects = "fixed",
                      link = "probit", 
                      likelihood = "bernoulli2",
                      regression = ~(durnpso + prevsys + bsa + weight + psa + age)*.trt,
                      class_interactions = "common",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_reg = normal(scale = 10),
                      init_r = 0.1,
                      QR = TRUE)
  }
  else if (sub_var == "male") {
    net_org <- add_integration(net_org,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               male = distr(qbern, prob = male),
                               n_int = 1000)
    FE_net_org <- nma(net_org, 
                      trt_effects = "fixed",
                      link = "probit", 
                      likelihood = "bernoulli2",
                      regression = ~(durnpso + prevsys + bsa + weight + psa + age)*.trt,
                      class_interactions = "common",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_reg = normal(scale = 10),
                      init_r = 0.1,
                      QR = TRUE)
    
  } 
  
  else if (sub_var == "male+age") {
    net_org <- add_integration(net_org,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               male = distr(qbern, prob = male),
                               age = distr(qgamma, mean = age_mean, sd = age_sd),
                               n_int = 1000)
    
    FE_net_org <- nma(net_org, 
                      trt_effects = "fixed",
                      link = "probit", 
                      likelihood = "",
                      regression = ~(durnpso + prevsys + bsa + weight + psa + age + male)*.trt,
                      class_interactions = "common",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_reg = normal(scale = 10),
                      init_r = 0.1,
                      QR = TRUE)
    
  }
  else {
    net_org <- add_integration(net_org,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               n_int = 1000)
    
    FE_net_org <- nma(net_org, 
                      trt_effects = "fixed",
                      link = "probit", 
                      likelihood = "bernoulli2",
                      regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
                      class_interactions = "common",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_reg = normal(scale = 10),
                      init_r = 0.1,
                      QR = TRUE)
  }
  
  return(FE_net_org)
  
}

#####################################################################
##                        Original scenario                        ##
##                                  vs                             ##
##                        3 IPD + 1 agg IPD + 5 agd                ##
##                                  vs                             ##
##                        3 IPD + 1 sub IPD + 5 agd                ##
#####################################################################

#####################################################################
##                        Original IPD+AGD                         ##
#####################################################################

## 4 IPD + 5 AgD
FE_net_org_df <- func.mult(sub_agd, "studyc", 
                           c("IXORA", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3"), 
                           c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                             "JUNCTURE"), 
                           "age")

#####################################################################
##                  Scenario 1：                                   ##
##                  3 IPD + 1 aggregated IPD + 5 AGD               ##
#####################################################################

#####################################################################
##                        selected IPD is IXORA                    ##
#####################################################################

## 3 IPD + 1 IPD aggregated + 5 AgD
FE_net_age_agd.sc1.1 <- func.mult(sub_agd1, "studyc", 
                                  c("UNCOVER-1", "UNCOVER-2", "UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "IXORA IPD"), 
                                  "age")

#####################################################################
##                        selected IPD is UNCOVER-1                ##
#####################################################################

## 3 IPD + 1 IPD aggregated + 5 AgD
FE_net_age_agd.sc2.1 <- func.mult(sub_agd1, "studyc", 
                                  c("IXORA", "UNCOVER-2", "UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "UNCOVER-1 IPD"), 
                                  "age")

#####################################################################
##                        selected IPD is UNCOVER-2                ##
#####################################################################

## 3 IPD + 1 IPD aggregated + 5 AgD
FE_net_age_agd.sc3.1 <- func.mult(sub_agd1, "studyc", 
                                  c("IXORA", "UNCOVER-1", "UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "UNCOVER-2 IPD"), 
                                  "age")

#####################################################################
##                        selected IPD is UNCOVER-3                ##
#####################################################################

## 3 IPD + 1 IPD aggregated + 5 AgD
FE_net_age_agd.sc4.1 <- func.mult(sub_agd1, "studyc", 
                                  c("IXORA", "UNCOVER-1", "UNCOVER-2"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "UNCOVER-3 IPD"), 
                                  "age")

#####################################################################
##                  Scenario 2：                                   ##
##                  3 IPD + 1 subgrouped IPD + 5 AGD               ##
#####################################################################

#####################################################################
##                        selected IPD is IXORA                    ##
#####################################################################

## 3 IPD + 1 IPD subgrouped + 5 AgD
FE_net_age_agd.sc1.2 <- func.mult(sub_agd1, "studyc", 
                                  c("UNCOVER-1", "UNCOVER-2", "UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", "JUNCTURE", "IXORA SUB"), "age")

#####################################################################
##                        selected IPD is UNCOVER-1                ##
#####################################################################

## 3 IPD + 1 IPD subgroup + 5 AgD
FE_net_age_agd.sc2.2 <- func.mult(sub_agd1, "studyc", 
                                  c("IXORA", "UNCOVER-2", "UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "UNCOVER-1 SUB"), 
                                  "age")

#####################################################################
##                        selected IPD is UNCOVER-2                ##
#####################################################################

## 3 IPD + 1 IPD subgroup + 5 AgD
FE_net_age_agd.sc3.2 <- func.mult(sub_agd1, "studyc", 
                                  c("IXORA", "UNCOVER-1", "UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "UNCOVER-2 SUB"), 
                                  "age")

#####################################################################
##                        selected IPD is UNCOVER-3                ##
#####################################################################

## 3 IPD + 1 IPD subgroup + 5 AgD
FE_net_age_agd.sc4.2 <- func.mult(sub_agd1, "studyc", 
                                  c("IXORA", "UNCOVER-1", "UNCOVER-2"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "UNCOVER-3 SUB"), 
                                  "age")


#####################################################################
##                        Original scenario                        ##
##                                  vs                             ##
##                        1 IPD + 3 agg IPD + 5 agd                ##
##                                  vs                             ##
##                        1 IPD + 3 sub IPD + 5 agd                ##
#####################################################################

#####################################################################
##    1 IPD + 3 agg IPD + 5 agd: selected IPD is UNCOVER-3         ##
#####################################################################

## 1 IPD + 3 IPD aggregated + 5 AgD
FE_net_age_agd.sc4.3 <- func.mult(sub_agd1, "studyc", 
                                  c("UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "IXORA IPD", "UNCOVER-1 IPD",
                                    "UNCOVER-2 IPD"), 
                                  "age")

#####################################################################
##    1 IPD + 3 sub IPD + 5 agd: selected IPD is UNCOVER-3         ##
#####################################################################

## 1 IPD + 3 subgrouped IPD + 5 AgD
FE_net_age_agd.sc4.4 <- func.mult(sub_agd1, "studyc", c("UNCOVER-3"), 
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE", 
                                    "JUNCTURE", "IXORA SUB", "UNCOVER-1 SUB", 
                                    "UNCOVER-2 SUB"), 
                                  "age")

#####################################################################
##                      Save the model result                      ##
#####################################################################

rlist::list.save(FE_net_org_df, 'D:\\MPH-Project\\FE_net_org_df.rdata')
rlist::list.save(FE_net_age_agd.sc1.1, 'D:\\MPH-Project\\FE_net_age_agd.sc1.1.rdata')
rlist::list.save(FE_net_age_agd.sc1.2, 'D:\\MPH-Project\\FE_net_age_agd.sc1.2.rdata')
rlist::list.save(FE_net_age_agd.sc2.1, 'D:\\MPH-Project\\FE_net_age_agd.sc2.1.rdata')
rlist::list.save(FE_net_age_agd.sc2.2, 'D:\\MPH-Project\\FE_net_age_agd.sc2.2.rdata')
rlist::list.save(FE_net_age_agd.sc3.1, 'D:\\MPH-Project\\FE_net_age_agd.sc3.1.rdata')
rlist::list.save(FE_net_age_agd.sc3.2, 'D:\\MPH-Project\\FE_net_age_agd.sc3.2.rdata')
rlist::list.save(FE_net_age_agd.sc4.1, 'D:\\MPH-Project\\FE_net_age_agd.sc4.1.rdata')
rlist::list.save(FE_net_age_agd.sc4.2, 'D:\\MPH-Project\\FE_net_age_agd_sc4.2.rdata')
rlist::list.save(FE_net_age_agd.sc4.3, 'D:\\MPH-Project\\FE_net_age_agd.sc4.3.rdata')
rlist::list.save(FE_net_age_agd.sc4.4, 'D:\\MPH-Project\\FE_net_age_agd_sc4.4.rdata')
