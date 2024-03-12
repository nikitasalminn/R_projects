install.packages("here", version='1.0.1')
install.packages("tidyverse", version='2.0.0')
install.packages("fixest", version='0.11.1')
install.packages("huxtable", version='5.5.2')
install.packages("flextable", version='0.8.6')



##############################################################
# Data preparation
##############################################################

# SETUP -------------------------------------------------------------------

# Load packages
library(tidyverse)
library(here) #easy file referencing, root is set at project root
library(fixest) # perform estimations with multiple fixed-effects
library(flextable) # Output tables to word
library(huxtable) # Output tables
library(dplyr)

#  Gravity data ------------------------------------------------------------

#Load Data 
data <- readRDS(here("input","Gravity.rds"))

# Select necessary data
data = data %>%
  filter(year %in% seq(1970, 2015))
  
#Save the dataset to .csv format
write_csv(data,file=here("output","Gravity.csv"))

##############################################################
# Data preparation
##############################################################

# SETUP -------------------------------------------------------------------
# Load data ---------------------------------------------------------------

d = read_csv(here("output","Gravity.csv"))

## Select necessary variables

x <- d %>%
  select(year, iso3_o, iso3_d, iso3num_o, iso3num_d, dist, comlang_off, comcol, contig, 
         pop_d, pop_o, gdp_d, gdp_o, eu_d, eu_o, rta_coverage,
         rta_type, tradeflow_comtrade_d, tradeflow_comtrade_o)

## remove rows with no data for trade flows and distance. Remove suspicious trade values from 0.00001 to 0.1

x_app1 <- x %>% 
  filter(!is.na(tradeflow_comtrade_d) & !is.na(tradeflow_comtrade_o)) %>%
  filter(!is.na(dist)) %>%
  filter(!(0.00001 <= tradeflow_comtrade_d & tradeflow_comtrade_d <= 0.1))%>%
  filter(!(0.00001 <= tradeflow_comtrade_o & tradeflow_comtrade_o <= 0.1))

#Since partner countries tend to report different values for same trade flow, we generally use 
#the mean of reported values when possible. If either country fails
#to report a value, we use the non-missing value.

x_app1 <- x_app1 %>%
  mutate(tradeflow = coalesce((tradeflow_comtrade_d + tradeflow_comtrade_o)/2, 
                              tradeflow_comtrade_d, 
                              tradeflow_comtrade_o))

## Construct symmetric pair id's
x_app1 = x_app1 %>%
  mutate(pair = paste(pmin(iso3_o, iso3_d),pmax(iso3_o, iso3_d),sep = "_")) %>%
  group_by(pair) %>%
  mutate(pair_id = cur_group_id())


## Construct assymmetric pair id's
x_app1 = x_app1 %>%
  mutate(pair_assym = paste(iso3_o, iso3_d,sep = "_")) %>%
  group_by(pair_assym) %>%
  mutate(pair_id_assym = cur_group_id())


## construct exporter output and importer expenditure
x_app1 = x_app1 %>%
  group_by(iso3_o,year) %>%
  mutate(Y_it = sum(tradeflow_comtrade_d)) %>%
  group_by(iso3_d,year) %>%
  mutate(E_jt = sum(tradeflow_comtrade_d))

## construct data for intranational trade (it creates measurment error)
#
#x_app1 = x_app1 %>%
#  group_by(iso3_o, year) %>%
#  mutate(intra = gdp_o - sum(tradeflow_comtrade_o))

## removing negative and missing intranational trade flows

#x_app1 = x_app1 %>% 
#  filter(!is.na(intra))

## calculate logs
x_app1 = x_app1 %>%
  mutate(across(c(tradeflow,Y_it,E_jt,dist),~log(.x),.names="ln_{.col}"))


#4. Running naive gravity models on our sample to check whether the main pattern remains the same


## Estimation: OLS, dep. var.: log trade (ln_tradeflow) ; explan. vars.: common border (contig) ; log distance (ln_dist) ; common language (comlang_off) ; common colony (comcol) ; log out put exporter (ln_Y_it) ; log expenditure importer (ln_E_jt)
fit_ols = feols(ln_tradeflow ~ ln_dist + contig + comlang_off + comcol + ln_Y_it + ln_E_jt,
                data = x_app1 %>%
                  filter(ln_tradeflow > 0 & iso3_o != iso3_d), 
                vcov = cluster ~ pair_id)

summary(fit_ols)

#results are in line with the theory and other emphirical research.

# 2. Proxy for multilateral resistance - despite we know that FE will suit better, we want to compare all possible ways of obtaining
# estimations for traditional gravity

## Calculate remoteness indices
x_app1 = x_app1 %>%
  group_by(year) %>%
  mutate(Y_t = sum(Y_it), E_t = sum(E_jt)) %>%
  group_by(iso3_o, year) %>%
  mutate(remoteness_exp = sum(dist /(E_jt / E_t))) %>%
  group_by(iso3_d, year) %>%
  mutate(remoteness_imp = sum(dist / (Y_it / Y_t))) %>%
  mutate(ln_remoteness_exp = log(remoteness_exp), 
         ln_remoteness_imp = log(remoteness_imp))


## Estimation with remotness
fit_remoteness = feols(ln_tradeflow ~ ln_dist + contig + comlang_off + comcol + ln_Y_it + 
                         ln_E_jt + ln_remoteness_exp + ln_remoteness_imp,
                       data = x_app1 %>% filter(tradeflow > 0),
                       vcov = cluster ~ pair_id)
summary(fit_remoteness)


# Coefficients of distance, contiguity and common languages are different relative to
#Naive Gravity
#- Not accounting for multilateral resistance biases the estimates!
#  - Remoteness coefficients are positive for export and neegative for import
#, Meanwhile, they are quite small and highly significant
#- Theory predicts that all else equal, more isolated/remote regions tend to trade
# more with each other.


# 3. Creating Fixed Effects -----------------------------------------------------------

## Create importer- exporter-year and country-pair specific fixed effects 
x_app1 = x_app1 %>%
  unite("fe_exp_year",c(iso3_o,year),sep="_",remove=FALSE) %>%
  unite("fe_imp_year",c(iso3_d,year),sep="_",remove=FALSE)%>%
  unite("fe_imp_exp",c(iso3_o,iso3_d),sep="_",remove=FALSE)
  

## OLS Estimation with FE
fit_fixedeffects = feols(ln_tradeflow ~ ln_dist + contig + comlang_off + comcol |
                           fe_exp_year + fe_imp_year,
                         data = x_app1,
                         vcov = cluster ~ pair_id)
summary(fit_fixedeffects)

# 4. Fit gravity with PPML and fixed effects -----------------------------
fit_poisson = fepois(tradeflow ~ ln_dist + contig + comlang_off + comcol |
                       fe_exp_year + fe_imp_year,
                     data = x_app1,
                     vcov = cluster ~ pair_id)
### should be added %>% filter(exporter != importer)
summary(fit_poisson)


# Overview  of traditional gravity results----------------------------------------------------------------
tab_traditional_gravity =  huxreg(" " = fit_ols,
                                  "Remoteness" = fit_remoteness,
                                  "Fixed Effects" = fit_fixedeffects,
                                  "Fixed Effects " = fit_poisson,
                                  coefs = c("Intercept" = "(Intercept)",
                                            "Log Distance" = "ln_dist",
                                            "Contiguity" = "contig",
                                            "Common language" = "comlang_off",
                                            "Colony" = "comcol",
                                            "Log output" = "ln_Y_it",
                                            "Log expenditure" = "ln_E_jt",
                                            "Exporter remoteness" = "ln_remoteness_exp",
                                            "Importer remoteness" = "ln_remoteness_imp"),
                                  note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1970-2015. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row("","(1) OLS", "(2) OLS", "(3) OLS", " (4) PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "No", "No", "Yes", "Yes"),
             after = 24) %>%
  insert_row(c("Importer-time fixed effects", "No", "No", "Yes", "Yes"),
             after = 25) %>%
  set_number_format(23:24, everywhere, 0) %>%
  set_tb_borders(24:25,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Traditional Gravity Estimates") %>%
  set_label("tab_traditional_gravity") 

width(tab_traditional_gravity) = 1 #Set relative table width for use in documents

## Export table to word
tab_traditional_gravity_docx = as_flextable(tab_traditional_gravity)
save_as_docx(tab_traditional_gravity_docx, path = here("output","tables","tab_traditional_gravity.docx"))

##################
#Analysis of average RTAs effects, including of lags. Splitting RTAs into sub-categories

x_app1$rta_type <- ifelse(is.na(x_app1$rta_type), 0, x_app1$rta_type)

x_app1 <- x_app1 %>%
  mutate(CU = ifelse(rta_type == "CU", 1, 0))
x_app1 <- x_app1 %>%
  mutate(FTA_EIA = ifelse(rta_type == "FTA & EIA", 1, 0))
x_app1 <- x_app1 %>%
  mutate(CU_EIA = ifelse(rta_type == "CU & EIA", 1, 0))
x_app1 <- x_app1 %>%
  mutate(FTA = ifelse(rta_type == "FTA", 1, 0))
x_app1 <- x_app1 %>%
  mutate(PSA_EIA = ifelse(rta_type == "PSA & EIA", 1, 0))
x_app1 <- x_app1 %>%
  mutate(RTA = ifelse(rta_type == "CU" | rta_type == "CU&EIA" | rta_type == "FTA & EIA" | rta_type == "FTA" | rta_type == "PSA & EIA" , 1, 0))


rta_ols1 = feols(ln_tradeflow ~  RTA|
                  fe_exp_year + fe_imp_year + pair_id, 
                data = x_app1 %>%
                  filter(iso3_o != iso3_d), 
                vcov = cluster ~ pair_id)
summary(rta_ols1)

# Getting lagged values, using function written by Ruben Dewitte
tlag <- function(x, n = 1L, along_with, default = NA) { 
  if (!is.numeric(n) | (length(n)>1)) stop("n must be a numeric of length one")
  index <- match(along_with - n, along_with, incomparables = NA)
  out <- x[index]
  if (!is.na(default)) out[which(is.na(index))] <- default
  out
}

x_app1 = x_app1 %>%
  group_by(iso3_o, iso3_d) %>%
  mutate(rta_lag4 = tlag(RTA,n=4,along_with = year),
         rta_lag8 = tlag(RTA,n=8,along_with = year),
         rta_lag12 = tlag(RTA,n=12,along_with = year))
# As R doesn't know values for years before 1970, it can't create 4,8, and 12 years lags up until 1974, 1978 and 1982. 
# we assume that most of the trade agreements are signed once (e.g it can't be RTA in 1967, if there are no
# agreement in 1970). So, we replace all NAs with 0s.
x_app1$rta_lag4 <- ifelse(is.na(x_app1$rta_lag4), 0, x_app1$rta_lag4)
x_app1$rta_lag8 <- ifelse(is.na(x_app1$rta_lag8), 0, x_app1$rta_lag8)
x_app1$rta_lag12 <- ifelse(is.na(x_app1$rta_lag12), 0, x_app1$rta_lag12)

### same for EEA

#x_app1 = x_app1 %>%
#  group_by(iso3_o, iso3_d) %>%
#  mutate(FTA_EIA_lag4 = tlag(FTA_EIA,n=4,along_with = year),
#         FTA_EIA_lag8 = tlag(FTA_EIA,n=8,along_with = year),
#         FTA_EIA_lag12 = tlag(FTA_EIA,n=12,along_with = year))
# As R doesn't know values for years before 1970, it can't create 4,8, and 12 years lags up until 1974, 1978 and 1982. 
# we assume that most of the trade agreements are signed once (e.g it can't be RTA in 1967, if there are no
# agreement in 1970). So, we replace all NAs with 0s.
#x_app1$FTA_EIA_lag4 <- ifelse(is.na(x_app1$FTA_EIA_lag4), 0, x_app1$FTA_EIA_lag4)
#x_app1$FTA_EIA_lag8 <- ifelse(is.na(x_app1$FTA_EIA_lag8), 0, x_app1$FTA_EIA_lag8)
#x_app1$FTA_EIA_lag12 <- ifelse(is.na(x_app1$FTA_EIA_lag12), 0, x_app1$FTA_EIA_lag12)



#Estimation of RTAs impact including "phase-in" effect. All controls are excluded because country-specific variations
#are captured by fixed-effects
rta_ols = feols(ln_tradeflow ~  RTA + rta_lag4 + rta_lag8 + rta_lag12|
                  fe_exp_year + fe_imp_year, 
                data = x_app1 %>%
                  filter(iso3_o != iso3_d), 
                vcov = cluster ~ pair_id)
summary(rta_ols)


rta_ols_ppml = fepois(tradeflow ~  RTA + rta_lag4 + rta_lag8 + rta_lag12|
                  fe_exp_year + fe_imp_year, 
                data = x_app1 %>%
                  filter(iso3_o != iso3_d), 
                vcov = cluster ~ pair_id)
summary(rta_ols_ppml)

#Estimation of FTA_EIA impact including "phase-in" effect. All controls are excluded because country-specific variations
#are captured by fixed-effects
#FTA_EIA_ols = feols(ln_tradeflow ~  FTA_EIA + FTA_EIA_lag4 + FTA_EIA_lag8 + FTA_EIA_lag12|
#                  fe_exp_year + fe_imp_year, 
#                data = x_app1 %>%
#                  filter(iso3_o != iso3_d), 
#                vcov = cluster ~ pair_id)
#summary(FTA_EIA_ols)
# results show that the joint value of RTA, 4-, 8-, and 12-year lags is around 2.3, while 
# the joint value of RTA, and 4-year lag is around 2.4

# Allowing for heterogeneity between trade agreements
types_rta_ols = feols(ln_tradeflow ~  CU + CU_EIA + FTA_EIA + FTA + PSA_EIA|
                  fe_exp_year + fe_imp_year, 
                data = x_app1 %>%
                  filter(iso3_o != iso3_d), 
                vcov = cluster ~ pair_id)
summary(types_rta_ols)

# all types of RTAs are positively affecting trade with CU being twice as efficient as other RTAs 
#In line with works with other authors

# 4. Testing for potential "reverse causality" --------

# To get lead values in function of CONSECUTIVE observations (here years)
tlead <- function(x, n = 1L, along_with, default = NA) { 
  if (!is.numeric(n) | (length(n)>1)) stop("n must be a numeric of length one")
  index <- match(along_with + n, along_with, incomparables = NA)
  out <- x[index]
  if (!is.na(default)) out[which(is.na(index))] <- default
  out
}

# ## Identify future RTAs

x_app1 = x_app1 %>%
  group_by(iso3_o, iso3_d) %>%
  mutate(rta_lead4 = tlead(RTA,n=4,along_with = year))
x_app1$rta_lead4 <- ifelse(is.na(x_app1$rta_lead4), 0, x_app1$rta_lead4)
## Estimation
rta_lead = fepois(tradeflow ~  RTA + rta_lead4 |
                    fe_exp_year + fe_imp_year + pair_id,
                  data = x_app1, 
                  vcov = cluster ~ pair_id)

summary(rta_lead)

# no reverse causality

# ## Identify future FTA_EIA

#x_app1 = x_app1 %>%
#  group_by(iso3_o, iso3_d) %>%
#  mutate(FTA_EIA_lead4 = tlead(FTA_EIA,n=4,along_with = year))
#x_app1$FTA_EIA_lead4 <- ifelse(is.na(x_app1$FTA_EIA_lead4), 0, x_app1$FTA_EIA_lead4)
## Estimation
#FTA_EIA_lead = fepois(tradeflow ~  FTA_EIA + FTA_EIA_lead4 |
#                    fe_exp_year + fe_imp_year + pair_id,
#                  data = x_app1, 
#                  vcov = cluster ~ pair_id)

#summary(FTA_EIA_lead)

# 6. Addressing globalization effects -------------------------------------

## Create time-specific international border variables - no data for intranational trade
#x_app1 = x_app1 %>%
#  mutate(D_inter = ifelse(iso3_o != iso3_d, 1, 0),
#         intl_brdr_year = paste0("intl_brdr_", year)) %>%
#  pivot_wider(names_from="intl_brdr_year",
#              values_from="D_inter",
#              values_fill = 0)


#rta_glob = fepois(tradeflow ~  RTA + rta_lag4 + rta_lag8 + rta_lag12 +
#                    intl_brdr_1990 + intl_brdr_1994 +
#                    intl_brdr_1998 + intl_brdr_2002 + intl_brdr_2006 |
#                    fe_exp_year + fe_imp_year + pair_id,
#                  data = x_app1, 
#                  vcov = cluster ~ pair_id)


#####################################################################################
#####################################################################################

#summary(rta_glob)

# Overview ------------------------------------------------------------------

## Overview
tab_RTA_leads_lags =  huxreg("(1) RTA OLS" = rta_ols1,
                           "(2) RTA phase-in" = rta_ols,
                           "(3) RTA phase-in PPML" = rta_ols_ppml,
                           "(4) RTA reverse causality" = rta_lead,
                           "(5) RTAs types" = types_rta_ols,
                           coefs = c("RTA" = "RTA",
                                     "RTA(t+4)" = "rta_lead4",
                                     "RTA(t-4)" = "rta_lag4",
                                     "RTA(t-8)" = "rta_lag8",
                                     "RTA(t-12)" = "rta_lag12",
                             "Custom Union" = "CU",
                                     "Custom Union & EIA" = "CU_EIA",
                                     "FTA & EIA" = "FTA_EIA",
                                     "FTA" = "FTA",
                                     "PSA & EIA" = "PSA_EIA"
                                     ),
                           note = "Notes: All estimates are obtained with data for the years 1970-2015 and use exporter-time and importer-time fixed effects. The estimates of the fixed effects are omitted for brevity. Standard errors are clustered by country pair and are reported in parentheses.") %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes", "Yes", "Yes", "Yes"),
             after = 20) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes", "Yes", "Yes", "Yes"),
             after = 21) %>%
  insert_row(c("Pair-id fixed effects", "Yes", "No", "No", "No", "No"),
             after = 22) %>%
  set_number_format(15:16, everywhere, 0) %>%
  set_tb_borders(16:18,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Estimating the Effects of Regional Trade Agreements") %>%
  set_label("tab_RTA_leads_lags")

width(tab_RTA_leads_lags) = 0.95 #Set relative table width for use in documents


## Export table to word
tab_RTA_leads_lags_docx = as_flextable(tab_RTA_leads_lags)
save_as_docx(tab_RTA_leads_lags_docx, path = here("output","tables","tab_RTA_leads_lags.docx"))




#B. Measuring the impact of EEA (European Economic Area) agreement
#via gravity estimations

#a. Dummy creation
#------------------------------------------------------------------------------------------------------------------

#Data set up for EEA specific dummies
unique(x_app1$iso3_o)
x <- unique(x_app1$iso3_o)
vector_ex_EU <- x[!x %in% c("AUT", "BEL", "DNK", "FIN", "FRA", "DEU",
                            "GBR","GRC", "IRL", "ITA", "LUX", "NLD", "PRT", 
                            "ESP", "SWE")]
vector_ex_EFTA <- x[!x %in% c("ISL", "LIE", "NOR")]
vector_ex_EUandEFTA <- x[!x %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                   "DEU", "GBR", "GRC", "IRL", "ITA", "LUX",
                                   "NLD", "PRT", "ESP", "SWE", "ISL",
                                   "LIE", "NOR")]

## Dummies for symmetric country-pair effects

x_app1$AUT_NOR_S <- ifelse(x_app1$pair == "AUT_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$BEL_NOR_S <- ifelse(x_app1$pair == "BEL_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$DNK_NOR_S <- ifelse(x_app1$pair == "DNK_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$FIN_NOR_S <- ifelse(x_app1$pair == "FIN_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$FRA_NOR_S <- ifelse(x_app1$pair == "FRA_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$DEU_NOR_S <- ifelse(x_app1$pair == "DEU_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$GRC_NOR_S <- ifelse(x_app1$pair == "GRC_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$IRL_NOR_S <- ifelse(x_app1$pair == "IRL_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$ITA_NOR_S <- ifelse(x_app1$pair == "ITA_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$LUX_NOR_S <- ifelse(x_app1$pair == "LUX_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$NLD_NOR_S <- ifelse(x_app1$pair == "NLD_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$PRT_NOR_S <- ifelse(x_app1$pair == "NOR_POR" & x_app1$year >= 1994, 1, 0)
x_app1$ESP_NOR_S <- ifelse(x_app1$pair == "ESP_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$SWE_NOR_S <- ifelse(x_app1$pair == "NOR_SWE" & x_app1$year >= 1994, 1, 0)
x_app1$GBR_NOR_S <- ifelse(x_app1$pair == "GBR_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$ISL_NOR_S <- ifelse(x_app1$pair == "ISL_NOR" & x_app1$year >= 1994, 1, 0)


x_app1$AUT_ISL_S <- ifelse(x_app1$pair == "AUT_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$BEL_ISL_S <- ifelse(x_app1$pair == "BEL_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$DNK_ISL_S <- ifelse(x_app1$pair == "DNK_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$FIN_ISL_S <- ifelse(x_app1$pair == "FIN_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$FRA_ISL_S <- ifelse(x_app1$pair == "FRA_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$DEU_ISL_S <- ifelse(x_app1$pair == "DEU_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$GRC_ISL_S <- ifelse(x_app1$pair == "GRC_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$IRL_ISL_S <- ifelse(x_app1$pair == "IRL_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$ITA_ISL_S <- ifelse(x_app1$pair == "ITA_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$LUX_ISL_S <- ifelse(x_app1$pair == "ISL_LUX" & x_app1$year >= 1994, 1, 0)
x_app1$NLD_ISL_S <- ifelse(x_app1$pair == "ISL_NED" & x_app1$year >= 1994, 1, 0)
x_app1$PRT_ISL_S <- ifelse(x_app1$pair == "ISL_PRT" & x_app1$year >= 1994, 1, 0)
x_app1$ESP_ISL_S <- ifelse(x_app1$pair == "ESP_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$SWE_ISL_S <- ifelse(x_app1$pair == "ISL_SWE" & x_app1$year >= 1994, 1, 0)
x_app1$GBR_ISL_S <- ifelse(x_app1$pair == "GBR_ISL" & x_app1$year >= 1994, 1, 0)


## Dummies for assymetric country-pair effects

x_app1$AUT_NOR <- ifelse(x_app1$pair_assym == "AUT_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$BEL_NOR <- ifelse(x_app1$pair_assym == "BEL_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$DNK_NOR <- ifelse(x_app1$pair_assym == "DNK_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$FIN_NOR <- ifelse(x_app1$pair_assym == "FIN_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$FRA_NOR <- ifelse(x_app1$pair_assym == "FRA_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$DEU_NOR <- ifelse(x_app1$pair_assym == "DEU_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$GRC_NOR <- ifelse(x_app1$pair_assym == "GRC_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$IRL_NOR <- ifelse(x_app1$pair_assym == "IRL_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$ITA_NOR <- ifelse(x_app1$pair_assym == "ITA_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$LUX_NOR <- ifelse(x_app1$pair_assym == "LUX_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$NLD_NOR <- ifelse(x_app1$pair_assym == "NLD_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$PRT_NOR <- ifelse(x_app1$pair_assym == "PRT_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$ESP_NOR <- ifelse(x_app1$pair_assym == "ESP_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$SWE_NOR <- ifelse(x_app1$pair_assym == "SWE_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$GBR_NOR <- ifelse(x_app1$pair_assym == "GBR_NOR" & x_app1$year >= 1994, 1, 0)
x_app1$ISL_NOR <- ifelse(x_app1$pair_assym == "ISL_NOR" & x_app1$year >= 1994, 1, 0)


x_app1$AUT_ISL <- ifelse(x_app1$pair_assym == "AUT_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$BEL_ISL <- ifelse(x_app1$pair_assym == "BEL_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$DNK_ISL <- ifelse(x_app1$pair_assym == "DNK_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$FIN_ISL <- ifelse(x_app1$pair_assym == "FIN_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$FRA_ISL <- ifelse(x_app1$pair_assym == "FRA_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$DEU_ISL <- ifelse(x_app1$pair_assym == "DEU_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$GRC_ISL <- ifelse(x_app1$pair_assym == "GRC_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$IRL_ISL <- ifelse(x_app1$pair_assym == "IRL_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$ITA_ISL <- ifelse(x_app1$pair_assym == "ITA_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$LUX_ISL <- ifelse(x_app1$pair_assym == "LUX_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$NLD_ISL <- ifelse(x_app1$pair_assym == "NLD_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$PRT_ISL <- ifelse(x_app1$pair_assym == "PRT_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$ESP_ISL <- ifelse(x_app1$pair_assym == "ESP_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$SWE_ISL <- ifelse(x_app1$pair_assym == "SWE_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$GBR_ISL <- ifelse(x_app1$pair_assym == "GBR_ISL" & x_app1$year >= 1994, 1, 0)
x_app1$ISL_NOR <- ifelse(x_app1$pair_assym == "NOR_ISL" & x_app1$year >= 1994, 1, 0)


x_app1$ISL_AUT <- ifelse(x_app1$pair_assym == "ISL_AUT"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_BEL <- ifelse(x_app1$pair_assym == "ISL_BEL"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_DNK <- ifelse(x_app1$pair_assym == "ISL_DNK"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_FIN <- ifelse(x_app1$pair_assym == "ISL_FIN"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_FRA <- ifelse(x_app1$pair_assym == "ISL_FRA"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_DEU <- ifelse(x_app1$pair_assym == "ISL_DEU"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_GRC <- ifelse(x_app1$pair_assym == "ISL_GRC"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_IRL <- ifelse(x_app1$pair_assym == "ISL_IRL"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_ITA <- ifelse(x_app1$pair_assym == "ITA_ISL"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_LUX <- ifelse(x_app1$pair_assym == "ISL_LUX"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_NLD <- ifelse(x_app1$pair_assym == "ISL_NLD"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_PRT <- ifelse(x_app1$pair_assym == "ISL_PRT"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_ESP <- ifelse(x_app1$pair_assym == "ISL_ESP"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_SWE <- ifelse(x_app1$pair_assym == "ISL_SWE"  & x_app1$year >= 1994, 1, 0)
x_app1$ISL_GBR <- ifelse(x_app1$pair_assym == "ISL_GBR"  & x_app1$year >= 1994, 1, 0)


x_app1$NOR_AUT <- ifelse(x_app1$pair_assym == "NOR_AUT" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_BEL <- ifelse(x_app1$pair_assym == "NOR_BEL" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_DNK <- ifelse(x_app1$pair_assym == "NOR_DNK" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_FIN <- ifelse(x_app1$pair_assym == "NOR_FIN" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_FRA <- ifelse(x_app1$pair_assym == "NOR_FRA" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_DEU <- ifelse(x_app1$pair_assym == "NOR_DEU" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_GRC <- ifelse(x_app1$pair_assym == "NOR_GRC" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_IRL <- ifelse(x_app1$pair_assym == "NOR_IRL" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_ITA <- ifelse(x_app1$pair_assym == "NOR_ITA" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_LUX <- ifelse(x_app1$pair_assym == "NOR_LUX" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_NLD <- ifelse(x_app1$pair_assym == "NOR_NLD" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_PRT <- ifelse(x_app1$pair_assym == "NOR_PRT" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_ESP <- ifelse(x_app1$pair_assym == "NOR_ESP" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_SWE <- ifelse(x_app1$pair_assym == "NOR_SWE" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_GBR <- ifelse(x_app1$pair_assym == "NOR_GBR" & x_app1$year >= 1994, 1, 0)
x_app1$NOR_ISL <- ifelse(x_app1$pair_assym == "NOR_ISL" & x_app1$year >= 1994, 1, 0)

x_app1$intraEU <- ifelse(x_app1$iso3_o %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                              "DEU", "GRC", "IRL", "ITA", "LUX",
                                              "NLD", "PRT", "ESP", "SWE", "GBR")
                         & x_app1$iso3_d %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                "DEU", "GRC", "IRL", "ITA", "LUX",
                                                "NLD", "PRT", "ESP", "SWE", "GBR")
                         & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                              "1999", "2000", "2001", "2002", "2003",
                                              "2004", "2005", "2006", "2007", "2008",
                                              "2009", "2010", "2011", "2012", "2013",
                                              "2014", "2015"),1,0)

# There are no trade observations for Liechtenstein neither in IMF nor in COMTRADE
#------------------------------------------------------------------------------------------------------------------


#Dummy 1: EEA specific estimates

x_app1$EEA = ifelse(x_app1$iso3_o %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                      "DEU", "GRC", "IRL", "ITA", "LUX",
                                                      "NLD", "PRT", "ESP", "SWE", "GBR", "ISL",
                                                      "LIE", "NOR")
                                 & x_app1$iso3_d %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                          "DEU", "GRC", "IRL", "ITA", "LUX",
                                                          "NLD", "PRT", "ESP", "SWE", "GBR", "ISL",
                                                          "LIE", "NOR")
                                 & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                        "1999", "2000", "2001", "2002", "2003",
                                                        "2004", "2005", "2006", "2007", "2008",
                                                        "2009", "2010", "2011", "2012", "2013",
                                                        "2014", "2015"),1,0)

x_app1$RTAexceptEEA = ifelse((x_app1$EEA == 0 & x_app1$RTA == 1),1,0)


#Estimations including the EEA dummies

EEA_effect = feols(ln_tradeflow ~ RTAexceptEEA + EEA |
                             fe_exp_year + fe_imp_year + pair_id,
                           data = x_app1 %>% filter(iso3_o != iso3_d),
                           vcov = cluster ~ pair_id)

summary(EEA_effect)



#Dummy 2: Norway and Iceland specific estimates
x_app1$NOR_EU_trade = ifelse((x_app1$iso3_o %in% "NOR"
                            & x_app1$iso3_d %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                   "DEU", "GRC", "IRL", "ITA", "LUX",
                                                   "NLD", "PRT", "ESP", "SWE", "GBR", "ISL")
                            & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                 "1999", "2000", "2001", "2002", "2003",
                                                 "2004", "2005", "2006", "2007", "2008",
                                                 "2009", "2010", "2011", "2012", "2013",
                                                 "2014", "2015")) | (x_app1$iso3_d %in% "NOR"
                                                                     & x_app1$iso3_o %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                                                            "DEU", "GRC", "IRL", "ITA", "LUX",
                                                                                            "NLD", "PRT", "ESP", "SWE", "GBR", "ISL")
                                                                     & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                                                          "1999", "2000", "2001", "2002", "2003",
                                                                                          "2004", "2005", "2006", "2007", "2008",
                                                                                          "2009", "2010", "2011", "2012", "2013",
                                                                                          "2014", "2015"))  ,1,0)

x_app1$ICL_EU_trade = ifelse((x_app1$iso3_o %in% "ISL"
                              & x_app1$iso3_d %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                     "DEU", "GRC", "IRL", "ITA", "LUX",
                                                     "NLD", "PRT", "ESP", "SWE", "GBR", "ISL")
                              & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                   "1999", "2000", "2001", "2002", "2003",
                                                   "2004", "2005", "2006", "2007", "2008",
                                                   "2009", "2010", "2011", "2012", "2013",
                                                   "2014", "2015")) | (x_app1$iso3_d %in% "ISL"
                                                                       & x_app1$iso3_o %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                                                              "DEU", "GRC", "IRL", "ITA", "LUX",
                                                                                              "NLD", "PRT", "ESP", "SWE", "GBR", "NOE", "ISL")
                                                                       & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                                                            "1999", "2000", "2001", "2002", "2003",
                                                                                            "2004", "2005", "2006", "2007", "2008",
                                                                                            "2009", "2010", "2011", "2012", "2013",
                                                                                            "2014", "2015"))  ,1,0)


ICE_NORtradewithEU_ols = feols(ln_tradeflow ~ RTAexceptEEA +  intraEU + ICL_EU_trade + NOR_EU_trade |
                             fe_exp_year + fe_imp_year + pair_id,
                           data = x_app1 %>% filter(iso3_o != iso3_d),
                           vcov = cluster ~ pair_id)

summary(ICE_NORtradewithEU_ols)


##Dummy 3: 
#Assymetric country-specific results


x_app1$NORtoEU_EXP = ifelse(x_app1$iso3_o %in% "NOR"
                                      & x_app1$iso3_d %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                             "DEU", "GRC", "IRL", "ITA", "LUX",
                                                             "NLD", "PRT", "ESP", "SWE", "GBR", "ISL", "NOR")
                                      & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                 "1999", "2000", "2001", "2002", "2003",
                                                 "2004", "2005", "2006", "2007", "2008",
                                                 "2009", "2010", "2011", "2012", "2013",
                                                 "2014", "2015"),1,0)

x_app1$ICEtoEU_EXP = ifelse(x_app1$iso3_o %in% "ISL"
                            & x_app1$iso3_d %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                   "DEU", "GRC", "IRL", "ITA", "LUX",
                                                   "NLD", "PRT", "ESP", "SWE", "GBR", "NOR", "ISL")
                            & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                 "1999", "2000", "2001", "2002", "2003",
                                                 "2004", "2005", "2006", "2007", "2008",
                                                 "2009", "2010", "2011", "2012", "2013",
                                                 "2014", "2015"),1,0)


x_app1$NORfromEU_IMP = ifelse(x_app1$iso3_d %in% "NOR"
                            & x_app1$iso3_o %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                   "DEU", "GRC", "IRL", "ITA", "LUX",
                                                   "NLD", "PRT", "ESP", "SWE", "GBR", "NOR", "ISL")
                            & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                 "1999", "2000", "2001", "2002", "2003",
                                                 "2004", "2005", "2006", "2007", "2008",
                                                 "2009", "2010", "2011", "2012", "2013",
                                                 "2014", "2015"),1,0)

x_app1$ICEfromEU_IMP = ifelse(x_app1$iso3_d %in% "ISL"
                            & x_app1$iso3_o %in% c("AUT", "BEL", "DNK", "FIN", "FRA",
                                                   "DEU", "GRC", "IRL", "ITA", "LUX",
                                                   "NLD", "PRT", "ESP", "SWE", "GBR", "NOR", "ISL")
                            & x_app1$year %in% c("1994", "1995", "1996", "1997", "1998",
                                                 "1999", "2000", "2001", "2002", "2003",
                                                 "2004", "2005", "2006", "2007", "2008",
                                                 "2009", "2010", "2011", "2012", "2013",
                                                 "2014", "2015"),1,0)


NORICE_effect = feols(ln_tradeflow ~ RTAexceptEEA + intraEU + NORtoEU_EXP + NORfromEU_IMP + ICEtoEU_EXP + ICEfromEU_IMP|
                     fe_exp_year + fe_imp_year + pair_id,
                   data = x_app1 %>% filter(iso3_o != iso3_d),
                   vcov = cluster ~ pair_id)
summary(NORICE_effect)



##Dummy 4: 
#Symetric countrypair-specific results

EEA_countrypair_sym = feols(ln_tradeflow ~ RTAexceptEEA + intraEU + AUT_ISL_S + BEL_ISL_S + DNK_ISL_S + FIN_ISL_S +
                             FRA_ISL_S + DEU_ISL_S + GBR_ISL_S + GRC_ISL_S + IRL_ISL_S + ITA_ISL_S +
                             LUX_ISL_S + NLD_ISL_S + PRT_ISL_S + ESP_ISL_S + SWE_ISL_S + AUT_NOR_S + BEL_NOR_S + DNK_NOR_S + FIN_NOR_S +
                             FRA_NOR_S + DEU_NOR_S + GBR_NOR_S + GRC_NOR_S + IRL_NOR_S + ITA_NOR_S +
                             LUX_NOR_S + NLD_NOR_S + PRT_NOR_S + ESP_NOR_S + SWE_NOR_S + ISL_NOR_S |
                             fe_exp_year + fe_imp_year + pair_id,
                           data = x_app1 %>% filter(iso3_o != iso3_d),
                           vcov = cluster ~ pair_id)
summary(EEA_countrypair_sym)

##Dummy 4: 
#Asymmetric country-pair-specific results


EEA_countrypair_assym = feols(ln_tradeflow ~ RTAexceptEEA + intraEU + ISL_AUT + ISL_BEL + ISL_DNK + ISL_FIN +
                              ISL_FRA + ISL_DEU + ISL_GBR + ISL_GRC + ISL_IRL + ISL_ITA + ISL_LUX + ISL_NLD +
                              ISL_PRT + ISL_ESP + ISL_SWE + NOR_AUT + NOR_BEL + NOR_DNK + NOR_FIN + NOR_FRA +
                              NOR_DEU + NOR_GBR + NOR_GRC + NOR_IRL + NOR_ITA + FRA_ISL + DEU_ISL + GBR_ISL +
                              NOR_LUX + NOR_NLD + NOR_PRT + NOR_ESP + NOR_SWE + AUT_ISL+ BEL_ISL + DNK_ISL + 
                              FIN_ISL + GRC_ISL + IRL_ISL + ITA_ISL + LUX_ISL + NLD_ISL + PRT_ISL + ESP_ISL +
                              SWE_ISL + AUT_NOR + BEL_NOR + DNK_NOR + FIN_NOR + FRA_NOR + DEU_NOR + GBR_NOR +
                              GRC_NOR + IRL_NOR + ITA_NOR +  LUX_NOR + NLD_NOR + PRT_NOR + ESP_NOR + SWE_NOR +
                              NOR_ISL + ISL_NOR| fe_exp_year + fe_imp_year + pair_id,
                            data = x_app1 %>% filter(iso3_o != iso3_d),
                            vcov = cluster ~ pair_id)
summary(EEA_countrypair_assym)

### country-pair results are visualized manually as they contain up to 65 variables

tab_EEA =  huxreg("(1) EEA average effect (OLS)" = EEA_effect,
                             "(2) Symmetric effects for Norway and Iceland (OLS)" = ICE_NORtradewithEU_ols,
                             "(3) Assymmetric effects for Norway and Iceland (OLS)" = NORICE_effect,
                             coefs = c("RTA excpet for EEA" = "RTAexceptEEA",
                                       "EEA average effect" = "EEA",
                                       "Intra EU-15 trade" = "intraEU",
                                       "Iceland-EEA trade" = "ICL_EU_trade",
                                       "Norway-EEA trade" = "NOR_EU_trade",
                                       "Norway export to EEA" = "NORtoEU_EXP",
                                       "Norway import from EEA" = "NORfromEU_IMP",
                                       "Iceland export to EEA " = "ICEtoEU_EXP",
                                       "Iceland import from EEA" = "ICEfromEU_IMP"),
                             note = "Notes: All estimates are obtained with data for the years 1970-2015 and use exporter-time and importer-time fixed effects. The estimates of the fixed effects are omitted for brevity. Standard errors are clustered by country pair and are reported in parentheses.") %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  set_number_format(15:16, everywhere, 0) %>%
  set_tb_borders(16:18,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Estimating the Effects of EEA") %>%
  set_label("tab_EEA")

width(tab_EEA) = 0.95 #Set relative table width for use in documents


## Export table to word
tab_EEA_docx = as_flextable(tab_EEA)
save_as_docx(tab_EEA_docx, path = here("output","tables","tab_EEA.docx"))



########################################## visualization in Tableau


x_app2 <- subset(d, (iso3_d == "AUT" | iso3_d == "BEL"|iso3_d == "DNK"|iso3_d == "FIN"|iso3_d == "FRA"|
                   iso3_d == "DEU"|iso3_d == "GRC"| iso3_d =="IRL"| iso3_d =="ITA"|iso3_d == "LUX"|
                   iso3_d == "NLD"| iso3_d =="PRT"| iso3_d =="ESP"| iso3_d == "SWE") &  (iso3_o == "ISL"|
                   iso3_o == "LIE"| iso3_o == "NOR") & (year >= 1980 & year <= 2010))


## remove rows with no data for trade flows and distance. Remove suspicious trade values from 0.00001 to 0.1

x_app2 <- x_app2 %>% 
  filter(!is.na(tradeflow_comtrade_d) & !is.na(tradeflow_comtrade_o)) %>%
  filter(!is.na(dist)) %>%
  filter(!(0.00001 <= tradeflow_comtrade_d & tradeflow_comtrade_d <= 0.1))%>%
  filter(!(0.00001 <= tradeflow_comtrade_o & tradeflow_comtrade_o <= 0.1))

#Since partner countries tend to report different values for same trade flow, we generally use 
#the mean of reported values when possible. If either country fails
#to report a value, we use the non-missing value.

x_app2 <- x_app2 %>%
  mutate(tradeflow = coalesce((tradeflow_comtrade_d + tradeflow_comtrade_o)/2, 
                              tradeflow_comtrade_d, 
                              tradeflow_comtrade_o))

x_app2<- x_app2 %>%
  select(year, iso3_o, iso3_d, tradeflow)


save(x_app2,file=here("output","trade-flows.Rdata"))

