########################################################################################################################

#Measuring the AKFTA's effect on the AKFTA member states' trade as a bloc as also on Korean trade

########################################################################################################################


#preparing the libraries to run our models

install.packages("dplyr")

library(dplyr)
library(readxl)
library(stargazer)
library(tidyverse)
library(here)
library(fixest) # perform estimations with multiple fixed-effects
library(flextable) # Output tables to word
library(huxtable)
library(haven)


#inserting our dataset, retrieved from the CEPPI database and cleaned for our purposes 

Cleaned_RIA_Dataset_Tradeflows <- read_dta("Cleaned RIA Dataset + Tradeflows.dta")
View(Cleaned_RIA_Dataset_Tradeflows)


# construction of a dataframe to tailor the dataset even more to our project

d<-Cleaned_RIA_Dataset_Tradeflows

data_1<- data.frame(d$year, d$iso3_o, d$iso3_d, d$country_exists_o, d$country_exists_d, d$dist,
                    d$comlang_off, d$comcol, d$col45, d$pop_d, d$pop_o, d$gdp_d, d$gdp_o, 
                    d$gatt_d, d$gatt_o, d$wto_d, d$wto_o, d$eu_d, d$eu_o, d$rta, d$rta_coverage,
                    d$rta_type, d$tradeflow_comtrade_d, d$tradeflow_comtrade_o, d$tradeflow_baci,
                    d$manuf_tradeflow_baci, d$tradeflow_imf_d, d$tradeflow_imf_o, d$contig)

View(data_1)


##### exclusion of missing data in order to run our models

data_1<-data_1[-which(is.na(data_1$d.dist)),]
data_1<-data_1[-which(is.na(data_1$d.tradeflow_comtrade_o)),]


#------------------------------------------------------------------------------------------------------------------

# A. preparation of our data in order to run the gravity models and baseline regressions

#We first want to note that our trade variable is :tradeflow_comtrade_o

# 1. Construction of symmetric pair id's
data_1 = data_1 %>%
  mutate(pair = paste(pmin(d.iso3_o,d.iso3_d),pmax(d.iso3_o,d.iso3_d),sep = "_")) %>%
  group_by(pair) %>%
  mutate(pair_id = cur_group_id())

# 2. construction of exporter (d.iso3_o) output (Y_it) and importer (d.iso3_d) expenditure (E_jt)
data_1 = data_1 %>%
  group_by(d.iso3_o,d.year) %>%
  mutate(Y_it = sum(d.tradeflow_comtrade_o)) %>%
  group_by(d.iso3_d,d.year) %>%
  mutate(E_jt = sum(d.tradeflow_comtrade_o))


# 3. calculation of logs for the variables d.tradeflow_comtrade_o ; Y_it ; E_jt ; d.dist
data_1 = data_1 %>%
  mutate(across(c(d.tradeflow_comtrade_o,Y_it,E_jt,d.dist),~log(.x),.names="ln_{.col}"))


#4. testing our data amendments to run gravity models

## Estimation: OLS, dep. var.: log trade (ln_d.tradeflow_comtrade_o) ; explan. vars.: common border (d.contig) ; log distance (ln_d.dist) ; common language (d.comlang_off) ; common colony (d.comcol) ; log out put exporter (ln_Y_it) ; log expenditure importer (ln_E_jt)
fit_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + ln_Y_it + ln_E_jt,
                data = data_1 %>%
                  filter(ln_d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d), 
                vcov = cluster ~ pair_id)

summary(fit_ols)



# 5. Proxy for multilateral resistance

# Despite applying a set of fixed effects for our final regressions on a later stage, for completeness, we consider it useful to construct multilateral resistence terms.

#a. Calculate remoteness indices

data_1 = data_1 %>%
  group_by(d.year) %>%
  mutate(Y_t = sum(Y_it, na.rm = T), E_t = sum(E_jt, na.rm = T)) %>%
  group_by(d.iso3_o, d.year) %>%
  mutate(remoteness_exp = sum(d.dist /(E_jt / E_t), na.rm = T)) %>%
  group_by(d.iso3_d, d.year) %>%
  mutate(remoteness_imp = sum(d.dist / (Y_it / Y_t), na.rm = T)) %>%
  mutate(ln_remoteness_exp = log(remoteness_exp), 
         ln_remoteness_imp = log(remoteness_imp))


## b. OLS-Estimation incl. multilateral resistance proxies
#dep. var.: log trade (ln_d.tradeflow_comtrade_o) ; explan. vars.: common border (d.contig) ; log distance (ln_d.dist) ; common language (d.comlang_off) ; common colony (d.comcol) ; log out put exporter (ln_Y_it) ; log expenditure importer (ln_E_jt) ; log multil. resistence exp. (ln_remoteness_exp) ; log multil. resistence imp. (ln_remoteness_imp)

fit_remoteness = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + ln_Y_it + 
                         ln_E_jt + ln_remoteness_exp + ln_remoteness_imp,
                       data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                       vcov = cluster ~ pair_id)
summary(fit_remoteness)

# 6. Fixed Effects

#a. Construction of fixed effects

data_1 = data_1 %>%
  unite("fe_exp_year",c(d.iso3_o,d.year),sep="_",remove=FALSE) %>%
  unite("fe_imp_year",c(d.iso3_d,d.year),sep="_",remove=FALSE)

#b.  OLS-Estimation excl. multilateral resistance proxies and instead using fixed effects exp. year and imp. year
#dep. var.: log trade (ln_d.tradeflow_comtrade_o) ; explan. vars.: common border (d.contig) ; log distance (ln_d.dist) ; common language (d.comlang_off) ; common colony (d.comcol) ; log out put exporter (ln_Y_it) ; log expenditure importer (ln_E_jt)

fit_fixedeffects = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol |
                           fe_exp_year + fe_imp_year,
                         data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                         vcov = cluster ~ pair_id)

summary(fit_fixedeffects)



# c. OLS-Estimation excl. multilateral resistance proxies and instead using fixed effects exp. year and imp. year
#dep. var.: log trade (ln_d.tradeflow_comtrade_o) ; explan. vars.: common border (d.contig) ; log distance (ln_d.dist) ; common language (d.comlang_off) ; common colony (d.comcol) ; log out put exporter (ln_Y_it) ; log expenditure importer (ln_E_jt)

fit_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol |
                       fe_exp_year + fe_imp_year,
                     data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                     vcov = cluster ~ pair_id)
summary(fit_poisson)

#exporting the results into the output depository in .docx formation

tab_traditional_gravity =  huxreg("OLS" = fit_ols,
                                  "Remoteness" = fit_remoteness,
                                  "Fixed Effects" = fit_fixedeffects,
                                  "Fixed Effects " = fit_poisson,
                                  coefs = c("Intercept" = "(Intercept)",
                                            "Contiguity" = "d.contig",
                                            "Log Distance" = "ln_d.dist",
                                            "Common language" = "d.comlang_off",
                                            "Colony" = "d.comcol",
                                            "Log output" = "ln_Y_it",
                                            "Log expenditure" = "ln_E_jt",
                                            "Exporter remoteness" = "ln_remoteness_exp",
                                            "Importer remoteness" = "ln_remoteness_imp"),
                                  note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
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
  set_caption("Table 1. Traditional Gravity Estimates") %>%
  set_label("tab_traditional_gravity") 
width(tab_traditional_gravity) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_traditional_gravity)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_traditional_gravity.docx")



#-------------------------------------------------------------------------------------------------------------------

#B. Measuring the impact of the AKFTA via gravity estimations - Inclusion of RTA dummies

#a. Dummy creation
#------------------------------------------------------------------------------------------------------------------

#Data set up for KOREAN specific dummies

x <- unique(data_1$d.iso3_o)
vector_ex_AKFTA_incl_KOR <- x[!x %in% c("KOR","IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")]
vector_ex_AKFTA_excl_KOR <- x[!x %in% c("IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")]

#------------------------------------------------------------------------------------------------------------------


#Dummy 1: Intra-bloc trade creation
#Dummy AKFTA_intra_bloc takes the value 1 if both partners belong to AKFTA & 0 otherwise
#If the dummy coefficient is found to be positive and statistically significant we can assume trade creation effects and that the AKFTA has promoted intra-bloc trade.


data_1$AKFTA_intra_bloc = ifelse(data_1$d.iso3_o %in% c("KOR", "IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")
                                 & data_1$d.iso3_d %in% c("KOR", "IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")
                                 & data_1$d.year %in% c("2007", "2008","2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"),1,0)

#Dummy 2: Development of AKFTA-bloc exports towards non-AKFTA member states/ third countries
#Dummy AKFTA_bloc_exp_CR_3C takes the value 1 if the exporter belongs to the AKFTA-bloc & the importer is a non-AKFTA member state/ third country
#If the dummy coefficient is found to be positive and statistically significant, we can assume export trade creation effects and that the AKFTA has promoted AKFTA-bloc exports to third countries.


data_1$AKFTA_bloc_exp_CR_3C = ifelse(data_1$d.iso3_o %in% c("KOR", "IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")
                                     & data_1$d.iso3_d %in% vector_ex_AKFTA_incl_KOR
                                     & data_1$d.year %in% c("2007", "2008","2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", 
                                                            "2017", "2018", "2019")
                                     ,1,0)

#Dummy 3: Development of AKFTA-bloc imports from non-AKFTA member states/ third countries
#Dummy AKFTA_bloc_exp_CR_3C takes the value 1 if the importer belongs to the AKFTA-bloc & the exporter is a non-AKFTA member state/ third country
#If the dummy coefficient is found to be positive and statistically significant, we can assume import trade creation effects and that the AKFTA has promoted AKFTA-bloc imports from third countries.


data_1$AKFTA_bloc_imp_CR_3C = ifelse(data_1$d.iso3_o %in% vector_ex_AKFTA_incl_KOR
                                     & data_1$d.iso3_d %in% c("KOR", "IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")
                                     & data_1$d.year %in% c("2007", "2008","2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", 
                                                            "2017", "2018", "2019")
                                     ,1,0)

##Dummy 4: Development of Korean exports to the remaining AKFTA members
#Dummy KOR_ASEAN_bloc_EXP_CR takes the value 1 if KOREA is the exporter and an AKFTA member is the importing country
#If the dummy coefficient is found to be positive and statistically significant we can assume that KOREA has increased its export to AKFTA-members due to the AKFTA

data_1$KOR_ASEAN_bloc_EXP_CR = ifelse(data_1$d.iso3_o %in% "KOR"
                                      & data_1$d.iso3_d %in% c("IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")
                                      & data_1$d.year %in% c("2007", "2008","2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", 
                                                             "2017", "2018", "2019")
                                      ,1,0)

##Dummy 5: Development of Korean imports from the remaining AKFTA members
#Dummy KOR_ASEAN_bloc_IMP_CR takes the value 1 if KOREA is the importer and an AKFTA member is the exporting country
#If the dummy coefficient is found to be positive and statistically significant we can assume that KOREA has increased its imports from AKFTA-members due to the AKFTA


data_1$KOR_ASEAN_bloc_IMP_CR = ifelse(data_1$d.iso3_d %in% "KOR"
                                      & data_1$d.iso3_o %in% c("IDN", "MYS", "PHL", "SGP", "THA", "BRN", "VNM", "LAO", "MMR", "KHM")
                                      & data_1$d.year %in% c("2007", "2008","2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", 
                                                             "2017", "2018", "2019")
                                      ,1,0)


##Dummy 6: Development of Korean imports from non-AKFTA members/ third countries
#Dummy AKFTA_KOR_CR_exp takes the value 1 if KOREA is the importer and non-AKFTA members/ third countries are the exporting countries
#If the dummy coefficient is found to be positive and statistically significant we can assume that KOREA has increased its imports from non AKFTA-members/third countries due to the AKFTA. The latter finding would imply import trade diversion.


data_1 = data_1 %>% mutate(AKFTA_KOR_CR_imp = ifelse((d.iso3_o %in% vector_ex_AKFTA_incl_KOR 
                                                      & d.iso3_d %in% "KOR" 
                                                      & d.year >= 2007), 1, 0))

#Dummy 7: Development of Korean exports to non-AKFTA members/ third countries
#Dummy AKFTA_KOR_CR_exp takes the value 1 if KOREA is the exporter and non-AKFTA members/ third countries are the importing countries
#If the dummy coefficient is found to be positive and statistically significant we can assume that KOREA has increased its exports to non AKFTA-members/third countries due to the AKFTA. The latter finding would imply export trade diversion.

data_1 = data_1 %>% mutate(AKFTA_KOR_CR_exp = ifelse((d.iso3_o %in% "KOR" 
                                                      & d.iso3_d %in% vector_ex_AKFTA_incl_KOR 
                                                      & d.year >= 2007), 1, 0))


#-----------------------------------------------------------------------------------------------------------------

#b. Estimations including the AKFTA dummies



#aa. (1) Estimating intra-bloc trade effects of the AKFTA, inclusion of dummy AKFTA_intra_bloc, exluding pair_id fixed effects

#(1) OLS

AKFTA_intra_bloc_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol+ AKFTA_intra_bloc |
                               fe_exp_year + fe_imp_year,
                             data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                             vcov = cluster ~ pair_id)

summary(AKFTA_intra_bloc_ols)


# The ols estimation indicates that the AKFTA has decreased intra-bloc trade by -0.358742. The results are statistically significant at a 0.1 level.


# (2) PPML

AKFTA_intra_bloc_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + AKFTA_intra_bloc |
                                    fe_exp_year + fe_imp_year,
                                  data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                  vcov = cluster ~ pair_id)

summary(AKFTA_intra_bloc_poisson)



# The PPML estimation indicates that the AKFTA has increased intra-bloc trade by 0.274339. The results are not statistically significant.


#aa. (2) Estimating intra-bloc trade effects of the AKFTA, inclusion of dummy AKFTA_intra_bloc, including pair_id fixed effects

#(1) OLS

AKFTA_intra_bloc_pair_id_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol +AKFTA_intra_bloc |
                                       fe_exp_year + fe_imp_year + pair_id,
                                     data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                                     vcov = cluster ~ pair_id)

summary(AKFTA_intra_bloc_pair_id_ols)

# The ols estimation including pair_id as fixed effects indicates that the AKFTA has decreased intra-bloc trade by -0.315477. The results are statistically significant at a 0.001 level.


# (2) PPML

AKFTA_intra_bloc_pair_id_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + AKFTA_intra_bloc |
                                            fe_exp_year + fe_imp_year + pair_id,
                                          data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                          vcov = cluster ~ pair_id)

summary(AKFTA_intra_bloc_pair_id_poisson)

#exporting the results into the output depository in .docx formation

tab_rta_1 =  huxreg("Fixed Effects (1)" = AKFTA_intra_bloc_ols,
                    "Fixed Effects (2)" = AKFTA_intra_bloc_poisson,
                    "Fixed Effects (3)" = AKFTA_intra_bloc_pair_id_ols,
                    "Fixed Effects (4)" = AKFTA_intra_bloc_pair_id_poisson,
                    coefs = c("Contiguity" = "d.contig",
                              "Log Distance" = "ln_d.dist",
                              "Common language" = "d.comlang_off",
                              "Colony" = "d.comcol",
                              "FTA intra-block" = "AKFTA_intra_bloc"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row(""," OLS", " PPML", " OLS", "  PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 16) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 17) %>%
  insert_row(c("Pair-id fixed effects", "No", "No", "Yes", "Yes"),
             after = 18) %>%
  set_number_format(15:16, everywhere, 0) %>%
  set_tb_borders(16:18,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Table 2. Estimating the Effects of Free Trade Agreement") %>%
  set_label("tab_traditional_gravity") 
width(tab_rta_1) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_1)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_1.docx")



#---------------------------------------------------------------------------------------------------------------------

#bb. (1) Estimating the AKFTA's effect for overall trade levels of the AKFTA members (bloc), inclusion of dummy AKFTA_intra_bloc & AKFTA_bloc_exp_CR_3C & AKFTA_bloc_imp_CR_3C ; exluding pair_id fixed effects

#(1) OLS

AKFTA_bloc_TRADE_total_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + AKFTA_intra_bloc + AKFTA_bloc_exp_CR_3C + AKFTA_bloc_imp_CR_3C |
                                     fe_exp_year + fe_imp_year,
                                   data = data_1 %>% filter(d.tradeflow_comtrade_o > 0),
                                   vcov = cluster ~ pair_id)

# The variable 'AKFTA_bloc_imp_CR_3C' has been removed because of collinearity.

summary(AKFTA_bloc_TRADE_total_ols)

# The ols estimation indicates that the AKFTA has increased intra-bloc trade by 0.635389. The results are statistically significant at a 0.01 level.
# Additionally, the ols estimation indicates that the AKFTA has increased AKFTA-bloc exports to third countries by 1.066354. The results are statistically significant at a 0 level.


# (2) PPML

AKFTA_bloc_TRADE_total_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + AKFTA_intra_bloc + AKFTA_bloc_exp_CR_3C + AKFTA_bloc_imp_CR_3C |
                                          fe_exp_year + fe_imp_year,
                                        data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                        vcov = cluster ~ pair_id)

# The variable 'AKFTA_bloc_imp_CR_3C' has been removed because of collinearity.

summary(AKFTA_bloc_TRADE_total_poisson)


# The poisson estimation indicates that the AKFTA has decreased intra-bloc trade by -1.211362. The results are statistically significant at a 0.001 level.
# Additionally, the poisson estimation indicates that the AKFTA has decreased AKFTA-bloc exports to third countries by -1.496762. The results are statistically significant at a 0 level.

#bb. (2) Estimating the AKFTA's effect for overall trade levels of the AKFTA members (bloc), inclusion of dummy AKFTA_intra_bloc & AKFTA_bloc_exp_CR_3C & AKFTA_bloc_imp_CR_3C ; including pair_id fixed effects

#(1) OLS

AKFTA_bloc_TRADE_total_pair_id_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + AKFTA_intra_bloc + AKFTA_bloc_exp_CR_3C + AKFTA_bloc_imp_CR_3C |
                                             fe_exp_year + fe_imp_year + pair_id,
                                           data = data_1 %>% filter(d.tradeflow_comtrade_o > 0),
                                           vcov = cluster ~ pair_id)

# Variables 'd.contig', 'ln_d.dist' and 3 others have been removed because of collinearity.

summary(AKFTA_bloc_TRADE_total_pair_id_ols)


# The ols estimation including pair_id as fixed effects indicates that the AKFTA has decreased intra-bloc trade by -0.255489. The results are not statistically significant.
# Additionally, the ols estimation including pair_id as fixed effects indicates that the AKFTA has increased AKFTA-bloc exports to third countries by 0.063819. The results are not statistically significant.


# (2) PPML

AKFTA_bloc_TRADE_total_pair_id_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + AKFTA_intra_bloc + AKFTA_bloc_exp_CR_3C + AKFTA_bloc_imp_CR_3C |
                                                  fe_exp_year + fe_imp_year + pair_id,
                                                data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                                vcov = cluster ~ pair_id)

# Variables 'ln_d.dist', 'd.comlang_off', 'd.comcol' and 'AKFTA_bloc_imp_CR_3C' have been removed because of collinearity.

summary(AKFTA_bloc_TRADE_total_pair_id_poisson)



# The poisson estimation including pair_id as fixed effects indicates that the AKFTA has decreased intra-bloc trade by -0.508626. The results are not statistically significant.
# Additionally, the poisson estimation including pair_id as fixed effects indicates that the AKFTA has decreased AKFTA-bloc exports to third countries by -0.465542. The results are not statistically significant.

#exporting the results into the output depository in .docx formation

tab_rta_2 =  huxreg("Fixed Effects (1)" = AKFTA_bloc_TRADE_total_ols,
                    "Fixed Effects (2)" = AKFTA_bloc_TRADE_total_poisson,
                    "Fixed Effects (3)" = AKFTA_bloc_TRADE_total_pair_id_ols,
                    "Fixed Effects (4)" = AKFTA_bloc_TRADE_total_pair_id_poisson,
                    coefs = c("Contiguity" = "d.contig",
                              "Log Distance" = "ln_d.dist",
                              "Common language" = "d.comlang_off",
                              "Colony" = "d.comcol",
                              "FTA intra-block" = "AKFTA_intra_bloc",
                              "AKFTA bloc export to RW" = "AKFTA_bloc_exp_CR_3C",
                              "RW export to AKFTA block" = "AKFTA_bloc_exp_CR_3C"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row(""," OLS", " PPML", " OLS", "  PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 18) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 19) %>%
  insert_row(c("Pair-id fixed effects", "No", "No", "Yes", "Yes"),
             after = 20) %>%
  set_number_format(17:18, everywhere, 0) %>%
  set_tb_borders(18:20,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Table 3. Estimating the Effects of Free Trade Agreement") %>%
  set_label("tab_traditional_gravity") 
width(tab_rta_2) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_2)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_2.docx")

#---------------------------------------------------------------------------------------------------------------------

#cc. (1) Estimating the AKFTA's effect on KOREAN trade levels with AKFTA members (bloc), inclusion of dummies KOR_ASEAN_bloc_EXP_CR & KOR_ASEAN_bloc_IMP_CR ; exluding pair_id fixed effects

#(1) OLS

KOR_RAKFTA_EXPIMP_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_IMP_CR |
                                fe_exp_year + fe_imp_year,
                              data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                              vcov = cluster ~ pair_id)

summary(KOR_RAKFTA_EXPIMP_ols)


# The ols estimation indicates that the AKFTA has increased KOREAN exports to the remaining AKFTA members by 0.180177. The results are not statistically significant.
# Additionally, the ols estimation indicates that the AKFTA has increased KOREAN imports from the remaining AKFTA members by 0.463406. The results are not statistically significant.


# (2) PPML

KOR_RAKFTA_EXPIMP_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_IMP_CR |
                                     fe_exp_year + fe_imp_year,
                                   data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                   vcov = cluster ~ pair_id)


summary(KOR_RAKFTA_EXPIMP_poisson)


# The poisson estimation indicates that the AKFTA has increased KOREAN exports to the remaining AKFTA members by 0.776633. The results are statistically significant at a 0 level.
# Additionally, the poisson estimation indicates that the AKFTA has increased KOREAN imports from the remaining AKFTA members by 0.578234. The results are statistically significant at a 0.01 level.


#cc. (2) Estimating the AKFTA's effect on KOREAN trade levels with AKFTA members (bloc), inclusion of dummies KOR_ASEAN_bloc_EXP_CR & KOR_ASEAN_bloc_IMP_CR ; including pair_id fixed effects

#(1) OLS

KOR_RAKFTA_EXPIMP_pair_id_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_IMP_CR |
                                        fe_exp_year + fe_imp_year + pair_id,
                                      data = data_1 %>% filter(d.tradeflow_comtrade_o > 0),
                                      vcov = cluster ~ pair_id)

# Variables 'd.contig', 'ln_d.dist' and 2 others have been removed because of collinearity.

summary(KOR_RAKFTA_EXPIMP_pair_id_ols)


# The ols estimation including pair_id as fixed effects indicates the AKFTA has decreased KOREAN exports to the remaining AKFTA members by -0.343903. The results are not statistically significant.
# Additionally, the ols estimation including pair_id as fixed effects indicates that the AKFTA has decreased KOREAN imports from the remaining AKFTA members by -0.002443. The results are not statistically significant.


# (2) PPML

KOR_RAKFTA_EXPIMP_pair_id_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_IMP_CR |
                                             fe_exp_year + fe_imp_year + pair_id,
                                           data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                           vcov = cluster ~ pair_id)

# Variables 'ln_d.dist', 'd.comlang_off' and 'd.comcol' have been removed because of collinearity.

summary(KOR_RAKFTA_EXPIMP_pair_id_poisson)


# The poisson estimation including pair_id as fixed effects indicates that the AKFTA has increased KOREAN exports to the remaining AKFTA members by 0.181685. The results are not statistically significant.
# Additionally, the poisson estimation including pair_id as fixed effects indicates that the AKFTA has decreased KOREAN imports from the remaining AKFTA members by -0.045247. The results are not statistically significant.


tab_rta_3 =  huxreg("Fixed Effects (1)" = KOR_RAKFTA_EXPIMP_ols,
                    "Fixed Effects (2)" = KOR_RAKFTA_EXPIMP_poisson,
                    "Fixed Effects (3)" = KOR_RAKFTA_EXPIMP_pair_id_ols,
                    "Fixed Effects (4)" = KOR_RAKFTA_EXPIMP_pair_id_poisson,
                    coefs = c("Contiguity" = "d.contig",
                              "Log Distance" = "ln_d.dist",
                              "Common language" = "d.comlang_off",
                              "Colony" = "d.comcol",
                              "Korean export to bloc AKFTA" = "KOR_ASEAN_bloc_EXP_CR",
                              "AKFTA bloc export to South Korea" = "KOR_ASEAN_bloc_IMP_CR"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row(""," OLS", " PPML", " OLS", "  PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 17) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 18) %>%
  insert_row(c("Pair-id fixed effects", "No", "No", "Yes", "Yes"),
             after = 19) %>%
  set_number_format(16:17, everywhere, 0) %>%
  set_tb_borders(17:19,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Table 4. Estimating the Effects of Free Trade Agreement") %>%
  set_label("tab_traditional_gravity") 
width(tab_rta_3) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_3)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_3.docx")
#---------------------------------------------------------------------------------------------------------------------

#dd. (1) Estimating the AKFTA's effect on KOREAN export levels towards the AKFTA members (bloc) & KOREAN exports towards and imports from third countries /trade diversion export and import , inclusion of dummies KOR_ASEAN_bloc_EXP_CR & AKFTA_KOR_CR_imp & AKFTA_KOR_CR_exp ; excluding pair_id fixed effects

#(1) OLS

KOR_RAKFTA_3CEXP_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol  + KOR_ASEAN_bloc_EXP_CR + AKFTA_KOR_CR_imp + AKFTA_KOR_CR_exp |
                               fe_exp_year + fe_imp_year,
                             data = data_1 %>% filter(d.tradeflow_comtrade_o > 0 & d.iso3_o != d.iso3_d),
                             vcov = cluster ~ pair_id)

summary(KOR_RAKFTA_3CEXP_ols)


# The ols estimation indicates that the AKFTA has increased KOREAN exports to the remaining AKFTA members by 1.076322. The results are statistically significant at a 0.1 level.
# Additionally, the ols estimation indicates that the AKFTA has decreased KOREAN imports from non-AKFTA members by -0.463362. The results are not statistically significant.
#Furthermore, the ols estimation indicates that the AKFTA has increased KOREAN exports towards non-AKFTA members by 0.984859. The results are not statistically significant.

# (2) PPML

KOR_RAKFTA_3CEXP_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + AKFTA_KOR_CR_imp + AKFTA_KOR_CR_exp |
                                    fe_exp_year + fe_imp_year,
                                  data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                  vcov = cluster ~ pair_id)

summary(KOR_RAKFTA_3CEXP_poisson)


# The poisson estimation indicates that the AKFTA has decreased KOREAN exports to the remaining AKFTA members by -1.420851. The results are statistically significant at a 0 level.
# Additionally, the poisson estimation indicates that the AKFTA has decreased KOREAN imports from non-AKFTA members by -0.579211. The results are statistically significant at a 0.01 level.
#Furthermore, the poisson estimation indicates that the AKFTA has decreased KOREAN exports towards non-AKFTA members by -2.219558. The results are statistically significant at a 0 level.


#dd. (2) Estimating the AKFTA's effect on KOREAN export levels with AKFTA members (bloc),& KOREAN exports towards and imports from third countries /trade diversion export and import , inclusion of dummies KOR_ASEAN_bloc_EXP_CR & AKFTA_KOR_CR_imp & AKFTA_KOR_CR_exp ; including pair_id fixed effects

#(1) OLS

KOR_RAKFTA_3CEXP_pair_id_ols = feols(ln_d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + AKFTA_KOR_CR_imp + AKFTA_KOR_CR_exp |
                                       fe_exp_year + fe_imp_year + pair_id,
                                     data = data_1 %>% filter(d.tradeflow_comtrade_o > 0),
                                     vcov = cluster ~ pair_id)

# Variables 'd.contig', 'ln_d.dist' and 2 others have been removed because of collinearity.

summary(KOR_RAKFTA_3CEXP_pair_id_ols)


# The ols estimation including pair_id as fixed effects indicates the AKFTA has decreased KOREAN exports to the remaining AKFTA members by -0.347584. The results are not statistically significant.
# Additionally, the ols estimation including pair_id as fixed effects indicates that the AKFTA has increased KOREAN imports from non-AKFTA members by 0.002326. The results are not statistically significant.
#Furthermore, the ols estimation including pair_id as fixed effects indicates that the AKFTA has decreased KOREAN exports towards non-AKFTA members by -0.003952. The results are not statistically significant.


# (2) PPML

KOR_RAKFTA_3CEXP_pair_id_poisson = fepois(d.tradeflow_comtrade_o ~ d.contig + ln_d.dist + d.comlang_off + d.comcol + KOR_ASEAN_bloc_EXP_CR + AKFTA_KOR_CR_imp + AKFTA_KOR_CR_exp |
                                            fe_exp_year + fe_imp_year + pair_id,
                                          data = data_1 %>% filter(d.iso3_o != d.iso3_d),
                                          vcov = cluster ~ pair_id)

# Variables 'ln_d.dist', 'd.comlang_off' and 'd.comcol' have been removed because of collinearity.

summary(KOR_RAKFTA_3CEXP_pair_id_poisson)


# The poisson estimation including pair_id as fixed effects indicates that the AKFTA has increased KOREAN exports to the remaining AKFTA members by 0.134743. The results are not statistically significant.
# Additionally, the poisson estimation including pair_id as fixed effects indicates that the AKFTA has increased KOREAN imports from non-AKFTA members by 0.045063. The results are not statistically significant.
#Furthermore, the poisson estimation including pair_id as fixed effects indicates that the AKFTA has decreased KOREAN exports towards non-AKFTA members by -0.047247. The results are not statistically significant.



# Overview, tabs ------------------------------------------------------------------

## Overview
tab_rta_4 =  huxreg("Fixed Effects (1)" = KOR_RAKFTA_3CEXP_ols,
                    "Fixed Effects (2)" = KOR_RAKFTA_3CEXP_poisson,
                    "Fixed Effects (3)" = KOR_RAKFTA_3CEXP_pair_id_ols,
                    "Fixed Effects (4)" = KOR_RAKFTA_3CEXP_pair_id_poisson,
                    coefs = c("Contiguity" = "d.contig",
                              "Log Distance" = "ln_d.dist",
                              "Common language" = "d.comlang_off",
                              "Colony" = "d.comcol",
                              "Export of the RW to Korea" = "AKFTA_KOR_CR_exp",
                              "Korean export to bloc AKFTA" = "KOR_ASEAN_bloc_EXP_CR",
                              "Korean export to RW" = "AKFTA_KOR_CR_imp"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row(""," OLS", " PPML", " OLS", "  PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 18) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 19) %>%
  insert_row(c("Pair-id fixed effects", "No", "No", "Yes", "Yes"),
             after = 20) %>%
  set_number_format(17:18, everywhere, 0) %>%
  set_tb_borders(18:20,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Table 5. Estimating the Effects of Free Trade Agreement") %>%
  set_label("tab_traditional_gravity") 


width(tab_rta_4) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_4)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_4.docx")


#---------------------------------------------------------------------------------------------------------------------

#C. Measuring the impact of the AKFTA via gravity estimations - Addressing potential endogeneity of the AKFTA-dummies

#1.The AKFTA's effect for overall trade levels of the AKFTA members (bloc), inclusion of dummy AKFTA_intra_bloc & AKFTA_bloc_exp_CR_3C & AKFTA_bloc_imp_CR_3C

# Estimation 
AKFTA_bloc_TRADE_total_endo = fepois(d.tradeflow_comtrade_o ~  AKFTA_intra_bloc + AKFTA_bloc_exp_CR_3C |
                                       fe_exp_year + fe_imp_year + pair_id,
                                     data = data_1,
                                     vcov = cluster ~ pair_id)

# The variable 'AKFTA_bloc_imp_CR_3C' has been removed because of collinearity.

summary(AKFTA_bloc_TRADE_total_endo)

# The endogeneity estimation indicates that the AKFTA has decreased intra-bloc trade by -0.508716. The results are not statistically significant.
# Additionally, the endogeneity estimation indicates that the AKFTA has decreased AKFTA-bloc exports towards third countries by -0.465638. The results are not statistically significant.


#2. Estimating the AKFTA's effect on KOREAN export levels with AKFTA members (bloc),& KOREAN exports towards and imports from third countries /trade diversion export and import , inclusion of dummies KOR_ASEAN_bloc_EXP_CR & AKFTA_KOR_CR_imp & AKFTA_KOR_CR_exp

KOR_RAKFTA_3CEXP_endo = fepois(d.tradeflow_comtrade_o ~  KOR_ASEAN_bloc_EXP_CR + AKFTA_KOR_CR_imp + AKFTA_KOR_CR_exp |
                                 fe_exp_year + fe_imp_year + pair_id,
                               data = data_1,
                               vcov = cluster ~ pair_id)

summary(KOR_RAKFTA_3CEXP_endo)

# The endogeneity estimation indicates that the AKFTA has increased KOREAN exports to the remaining AKFTA members by 0.13465. The results are not statistically significant.
# Additionally, the endogeneity estimation indicates that the AKFTA has increased KOREAN imports from non-AKFTA members by 0.04499. The results are not statistically significant.
#Furthermore, the endogeneity estimation indicates that the AKFTA has decreased KOREAN exports towards non-AKFTA members by -0.04728. The results are not statistically significant.


tab_rta_4 =  huxreg("Fixed Effects (1)" = AKFTA_bloc_TRADE_total_endo,
                    "Fixed Effects (2)" = KOR_RAKFTA_3CEXP_endo,
                    coefs = c("Export of the RW to Korea" = "AKFTA_KOR_CR_exp",
                              "Korean export to bloc AKFTA" = "KOR_ASEAN_bloc_EXP_CR",
                              "Korean export to RW" = "AKFTA_KOR_CR_imp"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row("", " PPML", "  PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes"),
             after = 14) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes"),
             after = 15) %>%
  insert_row(c("Pair-id fixed effects", "Yes", "Yes"),
             after = 16) %>%
  set_number_format(13:14, everywhere, 0) %>%
  set_tb_borders(14:16,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Table 7. Addressing potential endogeneity of the AKFTA-dummies") %>%
  set_label("tab_traditional_gravity") 


width(tab_rta_4) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_4)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_4.docx")





#---------------------------------------------------------------------------------------------------------------------

#D. Measuring the impact of the AKFTA via gravity estimations - Testing for potential "reverse causality" of the AKFTA-dummies

# 1. Identification of future RTAs tailored towards the intra-bloc creation dummy AKFTA_intra_bloc

data_1 = data_1 %>%
  group_by(d.iso3_o, d.iso3_d) %>%
  mutate(AKFTA_intra_bloc_lead4 = tlead(AKFTA_intra_bloc,n=4,along_with = d.year))

# Estimation
AKFTA_intra_bloc_lead = fepois(d.tradeflow_comtrade_o ~  AKFTA_intra_bloc + AKFTA_intra_bloc_lead4 |
                                 fe_exp_year + fe_imp_year + pair_id,
                               data = data_1, 
                               vcov = cluster ~ pair_id)

summary(AKFTA_intra_bloc_lead)


# The reverse causality estimation indicates that the AKFTA has decreased intra-bloc trade by -0.039372. The results are not statistically significant.
# Additionally, the reverse causality estimation indicates that the AKFTA will increase future intra-bloc trade (i.e. + 4 years) by 0.056054. The results are not statistically significant.



# 2. Identification of future RTAs tailored towards the KOREAN-export-towards-the-remaining-AKFTA-members dummy KOR_ASEAN_bloc_EXP_CR

data_1 = data_1 %>%
  group_by(d.iso3_o, d.iso3_d) %>%
  mutate(KOR_ASEAN_bloc_EXP_CR_lead4 = tlead(KOR_ASEAN_bloc_EXP_CR,n=4,along_with = d.year))

## Estimation
KOR_ASEAN_bloc_EXP_CR_lead = fepois(d.tradeflow_comtrade_o ~  KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_EXP_CR_lead4 |
                                      fe_exp_year + fe_imp_year + pair_id,
                                    data = data_1, 
                                    vcov = cluster ~ pair_id)

summary(KOR_ASEAN_bloc_EXP_CR_lead)


# The reverse causality estimation indicates that the AKFTA has increased KOREAN exports towards remaining AKFTA members by 0.189303. The results are statistically significant at a 0.1 level.
# Additionally, the reverse causality estimation indicates that the AKFTA will decrease future KOREAN exports towards remaining AKFTA members (i.e. + 4 years) by -0.011842. The results are not statistically significant.


tab_rta_6 =  huxreg("Fixed Effects: Exporter/importer-time fixed effects +pair-id" = AKFTA_bloc_TRADE_total_endo,
                    "Fixed Effects: Exporter/importer-time fixed effects +pair-id)" = KOR_RAKFTA_3CEXP_endo,
                    coefs = c("FTA intra-block" = "AKFTA_intra_bloc",
                              "FTA intra-block (t+4)" = "AKFTA_intra_bloc",
                              "Korean export to bloc AKFTA (t+4)" = "KOR_ASEAN_bloc_EXP_CR",
                              "Korean export to bloc AKFTA" = "KOR_ASEAN_bloc_EXP_CR"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  set_caption("Table 7. Estimating the Effects of Free Trade Agreement") %>%
  set_label("tab_traditional_gravity") 
width(tab_rta_6) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_6)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_6.docx")




#---------------------------------------------------------------------------------------------------------------------

#E. Measuring the impact of the AKFTA via gravity estimations - Allowing for potential non-linear and phasing-in effects of the AKFTA

# 1. Identifying past AKFTAs tailored towards the intra-bloc creation dummy AKFTA_intra_bloc

data_1 = data_1 %>%
  group_by(d.iso3_o,d.iso3_d) %>%
  mutate(AKFTA_intra_bloc_lag4 = lag(AKFTA_intra_bloc,n=4,along_with = d.year),
         AKFTA_intra_bloc_lag8 = lag(AKFTA_intra_bloc,n=8,along_with = d.year))

## Estimation
AKFTA_intra_bloc_lags = fepois(d.tradeflow_comtrade_o ~  AKFTA_intra_bloc + AKFTA_intra_bloc_lag4 + AKFTA_intra_bloc_lag8 |
                                 fe_exp_year + fe_imp_year + pair_id,
                               data = data_1, 
                               vcov = cluster ~ pair_id)

summary(AKFTA_intra_bloc_lags)

# The potential non-linear and phasing-in effects estimation indicates that the AKFTA has decreased intra-bloc trade by -0.015607. The results are not statistically significant.
# Additionally, the potential non-linear and phasing-in effects estimation indicates that the AKFTA has increased past intra-bloc trade (- 4 years) by 0.017722. The results are not statistically significant.
# Furthermore, the potential non-linear and phasing-in effects estimation indicates that the AKFTA has decreased past intra-bloc trade (- 8 years) by -0.076635. The results are not statistically significant.



## 1. Identifying past AKFTAs tailored towards the KOREAN-export-towards-the-remaining-AKFTA-members dummy KOR_ASEAN_bloc_EXP_CR

data_1 = data_1 %>%
  group_by(d.iso3_o,d.iso3_d) %>%
  mutate(KOR_ASEAN_bloc_EXP_CR_lag4 = lag(KOR_ASEAN_bloc_EXP_CR,n=4,along_with = d.year),
         KOR_ASEAN_bloc_EXP_CR_lag8 = lag(KOR_ASEAN_bloc_EXP_CR,n=8,along_with = d.year))

## Estimation
KOR_ASEAN_bloc_EXP_CR_lags = fepois(d.tradeflow_comtrade_o ~  KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_EXP_CR_lag4 + KOR_ASEAN_bloc_EXP_CR_lag8 |
                                      fe_exp_year + fe_imp_year + pair_id,
                                    data = data_1, 
                                    vcov = cluster ~ pair_id)

summary(KOR_ASEAN_bloc_EXP_CR_lags)

#exporting the results into the output depository in .docx formation

tab_KOR_ASEAN_bloc_EXP_CR_lags_RIA =  huxreg(" " = KOR_ASEAN_bloc_EXP_CR_lags, statistics = c(N = "nobs"))
tab_KOR_ASEAN_bloc_EXP_CR_lags_RIA_docx = as_flextable(tab_KOR_ASEAN_bloc_EXP_CR_lags_RIA)
save_as_docx(tab_KOR_ASEAN_bloc_EXP_CR_lags_RIA_docx, path = "/Users/lotts/Desktop/AKFTA_gravity_est_May_2022/output/tab_KOR_ASEAN_bloc_EXP_CR_lags_RIA.docx")

# The potential non-linear and phasing-in effects estimation indicates that the AKFTA has increased KOREAN exports towards remaining AKFTA members by 0.119949. The results are not statistically significant.
# Additionally, the potential non-linear and phasing-in effects estimation indicates that the AKFTA has increased past KOREAN exports towards remaining AKFTA members (- 4 years) by 0.107498. The results are statistically significant at a 0.1 level.
# Furthermore, the potential non-linear and phasing-in effects estimation indicates that the AKFTA has increased past KOREAN exports towards remaining AKFTA members (- 8 years) by 0.051552. The results are not statistically significant.


tab_rta_5 =  huxreg("Fixed Effects (1)" = AKFTA_bloc_TRADE_total_endo,
                    "Fixed Effects (2)" = KOR_RAKFTA_3CEXP_endo,
                    'Fixed Effects (3)' = AKFTA_intra_bloc_lags,
                    'Fixed Effects (4)' = KOR_ASEAN_bloc_EXP_CR_lags,
                    coefs = c("Export of the RW to Korea" = "AKFTA_KOR_CR_exp",
                              "Korean export to RW" = "AKFTA_KOR_CR_imp",
                              "FTA intra-block" = "AKFTA_intra_bloc",
                              "AKFTA bloc export to RW" = "AKFTA_bloc_exp_CR_3C",
                              "Korean export to bloc AKFTA" = "KOR_ASEAN_bloc_EXP_CR"),
                    note = "Notes: Statistics based on author's calculations. All estimates are obtained with data for the years 1986, 1990, 1994, 1998, 2002, and 2006. Columns (1)-(3) use the OLS estimator. Column (1) does not control for the multilateral resistances. Column (2) uses remoteness indexes to control for multilateral resistances. Column (3) uses importer-time and exporter-time fixed effects, whose estimates are omitted for brevity, to control for multilateral resistances. Finally, column (4) employs the PPML estimator. Standard errors are clustered by country pair and are reported in parentheses."
) %>%
  insert_row(""," OLS", " PPML", " OLS", "  PPML", after = 0) %>%
  set_top_border(1,everywhere,1) %>%
  set_align(1, everywhere, "center") %>%
  insert_row(c("Exporter-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 16) %>%
  insert_row(c("Importer-time fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 17) %>%
  insert_row(c("Pair-id fixed effects", "Yes", "Yes", "Yes", "Yes"),
             after = 18) %>%
  set_number_format(15:16, everywhere, 0) %>%
  set_tb_borders(16:18,everywhere,0) %>%
  set_tb_padding(0) %>%
  set_col_width(c(0.3,rep(0.7/4,4))) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Table 6. Identifying past AKFTAs tailored towards the KOREAN-export-towards-the-remaining-AKFTA-members") %>%
  set_label("tab_traditional_gravity") 
width(tab_rta_5) = 1 #Set relative table width for use in documents
tab_fit_ols_RIA_docx = as_flextable(tab_rta_5)
save_as_docx(tab_fit_ols_RIA_docx, path = "C:/Users/User/Desktop/hurt/tab_rta_5.docx")


#---------------------------------------------------------------------------------------------------------------------

#F. Measuring the impact of the AKFTA via gravity estimations - Addressing globalization effects which could affect our estimations of the AKFTA's impact

## Create time-specific international border variables
data_1 = data_1 %>%
  mutate(D_inter = ifelse(d.iso3_d != d.iso3_o, 1, 0),
         intl_brdr_year = paste0("intl_brdr_",d.year)) %>%
  pivot_wider(names_from="intl_brdr_year",
              values_from="D_inter",
              values_fill = 0)

#view(data_1) --> when showing the data again, the r-software always crashes.

# Estimation

rta_glob = fepois(d.tradeflow_comtrade_o ~  KOR_ASEAN_bloc_EXP_CR + KOR_ASEAN_bloc_EXP_CR_lag4 + KOR_ASEAN_bloc_EXP_CR_lag8 + intl_brdr_1997 + intl_brdr_2001 + intl_brdr_2005 +
                    intl_brdr_2009 + intl_brdr_2014 + intl_brdr_2018 |
                    fe_exp_year + fe_imp_year + pair_id,
                  data = data_1, 
                  vcov = cluster ~ pair_id)

summary(rta_glob)

#Due to removal of the for this regression important dummies collinearity (intl_brdr_1997 ; intl_brdr_2001 ; intl_brdr_2005 ; intl_brdr_2009 ; intl_brdr_2014 ; intl_brdr_2018), we were unfortunately unable to estimate the potential effects of gloablisation on the effect of the AKFTA.


Y_2004 %>% #new variable, wealth of household member 
  mutate(Share_KOR = 100*Trade_Value/sum(Trade_Value))-> Y_2004
data_2004<-Y_2004 %>% filter(ISO_3  == "idn"| ISO_3  =="mys"| ISO_3  == "phl"| ISO_3  ==  "sgp"| ISO_3  ==  "tha"| ISO_3  ==  "brn"| ISO_3  ==  "vnm"| ISO_3  ==  "lao"| ISO_3  ==  "mmr"| ISO_3  ==  "khm")


Y_210 %>% #new variable, wealth of household member 
  mutate(Share_KOR = 100*Trade_Value/sum(Trade_Value))-> Y_210
data_2010<-Y_210 %>% filter(ISO_3  == "idn"| ISO_3  =="mys"| ISO_3  == "phl"| ISO_3  ==  "sgp"| ISO_3  ==  "tha"| ISO_3  ==  "brn"| ISO_3  ==  "vnm"| ISO_3  ==  "lao"| ISO_3  ==  "mmr"| ISO_3  ==  "khm")

library(ggplot2)
DATA_1 <- cbind.data.frame (data_2010$ISO_3, data_2004$Share_KOR, data_2010$Share_KOR)
colnames(DATA_1) <- c("Importer", "Y_2004", "Y_2010")
DATA_1<-arrange(DATA_1,desc(Y_2004))
value<-c(1:length(DATA_1$Importer))
DATA <- cbind.data.frame (DATA_1, value)
colnames(DATA) <- c("Importer", "Y_2004", "Y_2010", "value")





library(grid)
g.mid<-ggplot(DATA,aes(x=1,y=reorder(Importer, -value)))+geom_text(aes(label=Importer))+
  geom_segment(aes(x=0.94,xend=0.96,yend=Importer))+
  geom_segment(aes(x=1.04,xend=1.065,yend=Importer))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))

g1 <- ggplot(data = DATA, aes(x = reorder(Importer, -value), y =Y_2004)) +
  geom_bar(stat = "identity") + ggtitle("Share of SK's total export which went to AKFTA member-state in Year 2004") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = DATA, aes(x = reorder(Importer, -value), y = Y_2010)) +xlab(NULL)+
  geom_bar(stat = "identity") + ggtitle("Share of SK's total export which went to AKFTA member-state in Year 2010") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  coord_flip()

library(gridExtra)
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))




#nominal values

Y_2004 %>% #new variable, wealth of household member 
  mutate(Share_KOR = 100*Trade_Value/sum(Trade_Value))-> Y_2004
data_2004<-Y_2004 %>% filter(ISO_3  == "idn"| ISO_3  =="mys"| ISO_3  == "phl"| ISO_3  ==  "sgp"| ISO_3  ==  "tha"| ISO_3  ==  "brn"| ISO_3  ==  "vnm"| ISO_3  ==  "lao"| ISO_3  ==  "mmr"| ISO_3  ==  "khm")


Y_210 %>% #new variable, wealth of household member 
  mutate(Share_KOR = 100*Trade_Value/sum(Trade_Value))-> Y_210
data_2010<-Y_210 %>% filter(ISO_3  == "idn"| ISO_3  =="mys"| ISO_3  == "phl"| ISO_3  ==  "sgp"| ISO_3  ==  "tha"| ISO_3  ==  "brn"| ISO_3  ==  "vnm"| ISO_3  ==  "lao"| ISO_3  ==  "mmr"| ISO_3  ==  "khm")




library(ggplot2)
DATA_1 <- cbind.data.frame (data_2010$ISO_3, data_2004$Trade_Value/1000000, data_2010$Trade_Value/1000000)
colnames(DATA_1) <- c("Importer", "Y_2004", "Y_2010")
DATA_1<-arrange(DATA_1,desc(Y_2004))
value<-c(1:length(DATA_1$Importer))
DATA <- cbind.data.frame (DATA_1, value)
colnames(DATA) <- c("Importer", "Y_2004", "Y_2010", "value")
sum(DATA$Y_2010)/sum(DATA$Y_2004)

library(grid)
g.mid<-ggplot(DATA,aes(x=1,y=reorder(Importer, -value)))+geom_text(aes(label=Importer))+
  geom_segment(aes(x=0.94,xend=0.96,yend=Importer))+
  geom_segment(aes(x=1.04,xend=1.065,yend=Importer))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))

g1 <- ggplot(data = DATA, aes(x = reorder(Importer, -value), y =Y_2004)) +
  geom_bar(stat = "identity") + ggtitle("SK's total export which went to AKFTA member-state in Year 2004") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = DATA, aes(x = reorder(Importer, -value), y = Y_2010)) +xlab(NULL)+
  geom_bar(stat = "identity") + ggtitle("SK's total export which went to AKFTA member-state in Year 2010") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  coord_flip()

library(gridExtra)
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
