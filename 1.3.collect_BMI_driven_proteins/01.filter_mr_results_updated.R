library(vroom)
library(tidyverse)
library(magrittr)

system('mkdir -p output/')

########################
# curate Step 1 MR results (BMI to proteins)
########################
nsum <- vroom('/scratch/butler/bangli.cao/1.1.BMI_to_proteins/output/summary_results.txt', delim = '\t')
#nsum <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/0.BMI_no_MHC_proxy/output/summary_results.txt', delim = '\t')
#nrow(nsum)
#Mine:48167112
#Satoshi:60208890


# IVW
nsum2 <- nsum %>% filter(method == 'Inverse variance weighted') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
#nrow(nsum2)
# I have same results from Satoshi's. 4907.
nsum2$shortname <- unlist(lapply(strsplit(nsum2$protein, "\\."), function(x)x[1])) 
nsum2$shortname <- unlist(lapply(strsplit(nsum2$protein, "[.]"), function(x)x[1]))
# I checked the shortnames from mine and sato; they are the same.
nsum2 %<>% dplyr::select(protein, shortname, method, nsnp, b, se, pval, lo_ci, up_ci) 
# sum2 contains 4670 unique proteins. why not 4907? -> because multiple aptamers target the same proteins (e.g. 4670 unique proteins)
# Are you sure? I already checked 4907 rows and the duplicates removed already.

# weighted median
nsum2_median <- nsum %>% filter(method == 'Weighted median') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
nsum2_median %<>% dplyr::select(protein, b, se, pval, lo_ci, up_ci) %>% dplyr::rename(b.median = b, se.median = se, lo_ci.median = lo_ci, up_ci.median = up_ci, pval.median = pval)

# weighted mode
nsum2_mode <- nsum %>% filter(method == 'Weighted mode') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
nsum2_mode %<>% dplyr::select(protein, b, se, pval, lo_ci, up_ci) %>% dplyr::rename(b.mode = b, se.mode = se, lo_ci.mode = lo_ci, up_ci.mode = up_ci, pval.mode = pval)

# egger
nsum2_egger <- nsum %>% filter(method == 'MR Egger') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
nsum2_egger %<>% dplyr::select(protein, b, se, pval, lo_ci, up_ci) %>% dplyr::rename(b.egger = b, se.egger = se, lo_ci.egger = lo_ci, up_ci.egger = up_ci, pval.egger = pval)

# hetero results for IVW
# hetero <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/0.BMI_no_MHC_proxy/output/summary_hetero.txt', delim = '\t') 
hetero <- vroom('/scratch/butler/bangli.cao/1.1.BMI_to_proteins/output/summary_hetero.txt', delim = '\t') 
hetero2 <- hetero %>% filter(!is.na(Q_pval))
hetero2 %<>% dplyr::select(protein, Q, Q_df, Q_pval, isquared)

# pleio results for IVW
# pleio <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/0.BMI_no_MHC_proxy/output/summary_pleio.txt', delim = '\t')
pleio <- vroom('/scratch/butler/bangli.cao/1.1.BMI_to_proteins/output/summary_pleio.txt', delim = '\t')
pleio2 <- pleio %>% filter(!is.na(pval))
# pleio2$shortname <- unlist(lapply(strsplit(pleio2$protein, "[.]"), function(x)x[1]))
pleio2 %<>% dplyr::select(protein, egger_intercept, se, pval) %>% dplyr::rename(se.egger_intercept = se, pval.egger_intercept = pval)

# joining all IVW data first
njoin <- left_join(nsum2, hetero2, by='protein') %>% left_join(., pleio2, by = 'protein') # nested left njoin
#the dot (.) is used as a placeholder representing the current data being processed in a pipeline of commands

# then join median, mode, and egger results
njoin2 <- left_join(njoin, nsum2_median, by='protein') %>% left_join(., nsum2_mode, by = 'protein') %>% left_join(., nsum2_egger, by = 'protein') # nested left njoin
#4907

# remove duplicated (to ensure there is no duplicated protein)
njoin3 <- njoin2 %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein))
#4907, no duplicate

########################
# annotate pass or not
########################
# Significance
#(sig_threshold <- signif(0.05/4907, digits = 3))

njoin4 <- njoin3 %>% 
  mutate(Significance = case_when(
    pval < 1e-5  ~ 'Bonferroni',
    pval < 0.05 ~ 'Nominal',
    TRUE ~ 'No'
  ))
#4907

njoin4$Significance <- factor(njoin4$Significance, levels = c('No', 'Nominal', 'Bonferroni'))
table(njoin4$Significance)
# No    Nominal Bonferroni
# 1187       2412       1308


# heterogeneity
njoin4 %<>% mutate(Heterogeneity_test = case_when(isquared < 50 ~ 'Pass',
                                                  TRUE ~ 'Fail')) # 50%
#4907

# count
table(njoin4[c('Significance', 'Heterogeneity_test')])
# My Heterogeneity_test
# Significance Fail Pass
# No           14 1173
# Nominal      11 2401
# Bonferroni    0 1308



# pleiotropy
njoin4 <- njoin4 %>%
  mutate(Pleiotropy_test = case_when(
    (pval.egger_intercept >= 0.05) ~ 'Pass',
    TRUE ~ 'Fail'
  ))

# count
table(njoin4$Pleiotropy_test)
# Mine results
# Fail Pass
# 211 4696


 
# combine causal estimate, Heterogeneity_test, and Pleiotropy_test
njoin4 %>% filter(Significance == "Bonferroni") %>% filter(Heterogeneity_test == 'Pass') %>% nrow() #1304 # My results are 1308
njoin4 %>% filter(Significance == "Bonferroni") %>% filter(Heterogeneity_test == 'Pass') %>% filter(Pleiotropy_test == 'Pass') %>% nrow() #1229 # My results are 1234

# sensitivity test
# combine pass/fail results for Heterogeneity test and Pleiotropy test
njoin4 %<>% mutate(Sensitivity_test = case_when(Heterogeneity_test == 'Pass' & Pleiotropy_test == 'Pass' ~ 'Pass',
                                                TRUE ~ 'Fail'))

# add seqid
njoin5 <- njoin4 %>%
  separate(protein, into = c('shortname', 'seqid'), sep = "\\.", remove = FALSE) %>%
  select(protein, seqid, shortname, everything())
njoin5 <- njoin4 %>%
  separate(protein, into = c('shortname', 'seqid'), sep = "[.]", remove = FALSE) %>%
  select(protein, seqid, shortname, everything())
# njoin5 <- nsum2 %>% separate(protein, into = c('seqid1', 'seqid2', 'shortname'), sep = '\\.', remove = F) %>% unite(seqid, c('seqid1', 'seqid2'), sep = '_') %>% unite(protein, c('shortname', 'seqid'), sep = '.', remove = F)

# write.table(
#   njoin5,
#   file = '1.aptamer_summary_without_reverse_filter.tsv',
#   sep = '\t',
#   quote = F,
#   row.names = F
# )

########################
# To incorporate reverse causation results,
# curate reverse MR results (proteins to BMI)
########################

# collect reverse MR results (proteins to BMI)
# IVW or wald
# reverse <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.1.pQTL_to_BMI_opengwas/output/summary_results.txt', delim = '\t')
reverse <- vroom('/scratch/butler/bangli.cao/1.2.BMI_reverse/output/summary_results.txt', delim = '\t')
reverse2 <- reverse %>% filter(method == 'Inverse variance weighted' | method == 'Wald ratio') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
# Mine 
# 757 lines

# Satoshi's
# 357 lines

reverse2$shortname <- unlist(lapply(strsplit(reverse2$protein, "[.]"), function(x)x[1]))  
reverse2$seqid <- unlist(lapply(strsplit(reverse2$protein, "[.]"), function(x)x[2]))  
reverse2 %<>% dplyr::select(protein, seqid, shortname, method, nsnp, b, se, pval, lo_ci, up_ci) 

# get hetero results for IVW
# hetero <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.1.pQTL_to_BMI_opengwas/output/summary_hetero.txt', delim = '\t') 
hetero <- vroom('/scratch/butler/bangli.cao/1.2.BMI_reverse/output/summary_hetero.txt', delim = '\t')
hetero2 <- hetero %>% filter(!is.na(Q_pval))
hetero2 %<>% dplyr::select(protein, Q, Q_df, Q_pval, isquared) 

# get pleio results for IVW
# pleio <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.1.pQTL_to_BMI_opengwas/output/summary_pleio.txt', delim = '\t')
pleio <- vroom('/scratch/butler/bangli.cao/1.2.BMI_reverse/output/summary_pleio.txt', delim = '\t')
pleio2 <- pleio %>% filter(!is.na(pval))
pleio2 %<>% dplyr::select(protein, egger_intercept, se, pval) %>% dplyr::rename(se.egger_intercept = se, pval.egger_intercept = pval)

# combine results
# reverse_joining all IVW data first
reverse_join2 <- left_join(reverse2, hetero2, by='protein') %>% left_join(., pleio2, by = 'protein') # nested left reverse_join

# remove duplicated (to ensure there is no duplicated protein)
reverse_join3 <- reverse_join2 %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
# Mine still 757 lines
# Satoshi's still 757 lines

# annotate
# Significance
reverse_join4 <- reverse_join3 %>%  mutate(Significance = case_when(pval < 0.05/nrow(reverse_join3) ~ 'Bonferroni', 
                                                    TRUE ~ 'No'))
reverse_join4$Significance <- factor(reverse_join4$Significance, levels = c('Bonferroni', 'No'))


# sensitivity: heterogeneity
reverse_join4 %<>% mutate(Heterogeneity_test = case_when(
  isquared >= 0.5 & Q_pval < 0.05 ~ 'Fail',
  isquared < 0.5 | Q_pval >0.05 ~ 'Pass',
  is.na(isquared) ~ 'NA'))
table(reverse_join4[c('Significance', 'Heterogeneity_test')])
# My results
#             Heterogeneity_test
# Significance Fail  NA Pass
# Bonferroni    0  17   10
# No           64 319  347

# Satoshi's results
# Heterogeneity_test
# Significance Fail  NA Pass
# Bonferroni    0  10    3
# No           32 139  173


# sensitivity: pleiotropy
reverse_join4 %<>% mutate(Pleiotropy_test = case_when(pval.egger_intercept< 0.05  ~ 'Fail',
                                              pval.egger_intercept >= 0.5 ~ 'Pass',
                                             TRUE ~ 'NA'))

# Combine heterogeneity test and pleiotropy test results
reverse_join5 <- reverse_join4 %<>% mutate(Sensitivity_test = case_when((Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA') & (Pleiotropy_test == 'Pass' | Pleiotropy_test == 'NA') ~ 'Pass',
                                                       TRUE ~ 'Fail'))
# count
table(reverse_join5[c('Significance', 'Sensitivity_test')])
# Mine
#             Sensitivity_test
# Significance Fail Pass
# Bonferroni    0   27
# No           64  666

# Satoshi's
# Sensitivity_test
# Significance Fail Pass
# Bonferroni    0   13
# No           32  312



# finally, make a list of proteins with reverse causation
revese_protein <- reverse_join5 %>% filter(Significance == 'Bonferroni' & Sensitivity_test == 'Pass') %>% dplyr::select(protein) 
# Mine 27 lines
# Satoshi's 13 lines

########################
# incorporate reverse MR results into step 1 MR results
########################

# finally, for BMI-to-proteins, get proteins that are significant, non-heterogeneous, and non-pleiotropic  
sum5 <- njoin5 %>% mutate(Reverse_causation_test = case_when(protein %in% revese_protein$protein ~ 'reverse',
                                                                TRUE ~ 'Pass'))

# check 
table(sum5$Reverse_causation_test)
# Pass reverse
# 4881      26

# 26 lines from reverse are in the 4907 lines.

# sensitvivity test: consistency, heterogeneity, and reverse causation
sum5 %<>% mutate(Sensitivity_test = case_when(Heterogeneity_test == 'Pass' & Pleiotropy_test == 'Pass' & Reverse_causation_test== 'Pass' ~ 'Pass',
                                                      TRUE ~ 'Fail'))
sum5$Sensitivity_test <- factor(sum5$Sensitivity_test, levels = c('Pass', 'Fail'))

# significance + sensitivity
sum5 %<>% mutate(Significance_Sensitivity = case_when((Significance == 'Bonferroni') & Sensitivity_test  == 'Pass' ~ 'Pass',
                                                      TRUE ~ 'Fail'))

# count
table(sum5$Significance_Sensitivity)
# Mine results
# Fail Pass
# 3681 1226


table(sum5[c('Significance', 'Sensitivity_test')])
# Mine
# Sensitivity_test
# Significance Pass Fail
# No         1144   43
# Nominal    2278  134
# Bonferroni 1226   82



write_tsv(sum5, file = '/scratch/butler/bangli.cao/1.3.collect_BMI_driven_proteins/output/aptamer_summary.tsv')
# Satoshi's path:/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.1.pQTL_to_BMI_opengwas/output/aptamer_summary.tsv


