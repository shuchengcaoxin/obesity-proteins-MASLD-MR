library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)

#####################
# Reverse MR  (proteins to BMI)
#####################

wd <- "/scratch/butler/bangli.cao/1.2.BMI_reverse/" 
setwd(wd)
output_dir <- paste0(wd, "output/") 
# paste0() this function concatenates the input values in a single character string

# read the outcome GWAS
outcome_path <- '/scratch/richards/satoshi.yoshiji/09.proMR/01.exposure/BMI_opengwas/ukb-b-19953.tsv' # Use UKB BMI here 
outcome_GWAS <- vroom(outcome_path)
outcome_GWAS <- outcome_GWAS %>% dplyr::rename(SNP = ID) 
# rename(), new_name = old_name

args <- commandArgs(trailingOnly=TRUE) 
# Extract command-line arguments; trailingOnly = TRUE, only arguments after '--args' are returned
# If '--args' is not explicitly mentioned in your command, 'trailingOnly =TRUE' simply captures all arguments after the script name
exp_path <- args[1]
protname <- args[2]

#exposure
exp_dat <- read_exposure_data(filename = exp_path,
                             sep='\t', # specify delimeter in file
                             snp_col='variant', # required name of column with SNP rs IDs
                             beta_col = 'beta_unadj', # Required for MR. Name of column with effect sizes
                             se_col = 'se', # Required for MR. Name of column with standard errors
                             effect_allele_col = 'Amin', # Name of column with effect allele
                             other_allele_col = 'Amaj', # Name of column with non effect allele
                             eaf_col = 'MAF', # Name of column with effect allele frequency
                             pval_col = 'pval') # Required for enrichment tests. Name of column with p-value

# read the outcome
formatted_outcome <- format_data(outcome_GWAS, 
                                 snps = exp_dat$SNP, 
                                 type="outcome", 
                                 snp_col="SNP", 
                                 beta_col="ES", 
                                 se_col = "SE", 
                                 eaf_col = "AF", 
                                 effect_allele_col = "ALT", 
                                 other_allele_col = "REF", 
                                 pval_col = "pval", 
                                 chr_col = "seqnames", 
                                 pos_col = "start"
                                 )
# format_data(): read exposure or outcome data; checks and organises columns fro use with MR or
# enrichment tests. Infers p-values when possible from beta and se.

formatted_outcome$id.outcome <- 'outcome'

# For proxy search, you may use snappy v1.0 (https://gitlab.com/richards-lab/vince.forgetta/snappy/-/blob/master/snappy)

# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)

# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))
exp_data_outcome_name <- paste0(output_dir,"harmonized/", protname, ".harmonized.txt")
write_tsv(exp_dat_outcome, file=exp_data_outcome_name)

# mr result
mr_results <- mr(exp_dat_outcome) #mr() perform all Mendelian randomization tests
mr_name <- paste0(output_dir, "mr/", protname, ".mr.txt") 
write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write_tsv(OR, file=OR_name)

# # scatter plot
# pdf_name <- paste0(output_dir, "pdf/", protname, ".scatter.pdf")
# pdf(pdf_name)
# mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1)
# dev.off()

# horizontal pleiotropy
pleio_res <- mr_pleiotropy_test(exp_dat_outcome) # Performs MR Egger and returns intercept values.
pleio_name <- paste0(output_dir,"pleio/", protname, ".pleio.txt")
write_tsv(pleio_res, file=pleio_name)

# hetero test
tryCatch({
  res_single <- mr_singlesnp(exp_dat_outcome)
  res_single_beta <- res_single$b
  res_single_se <- res_single$se
  res_Isq <- Isq(res_single_beta, res_single_se)
  hetero_res <- mr_heterogeneity(exp_dat_outcome) # Get heterogeneity statistics
  hetero_res$i_squared <- res_Isq
  hetero_name <- paste0(output_dir,"hetero/", protname, ".hetero.txt")
  write_tsv(hetero_res, file=hetero_name)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# tryCatch() to handle potential errors during the execution of Mendelian randomization (MR) analysis
# `error=` defines what to do if an error occurs in the try block. In this case, it prints an error message
# using `cat`, displaying the condition message associated with the error.
# `e` represents the error object
# conditionMessage() is a function in R that extracts the message from a condition object. 
# In this case, it's used to get a human-readable message from the error object `e`.
# cat() concatenates and prints objects.
# "\n" is a newline character.

# mr_singlesnp(): perform 2 sample MR on each SNP individually
# Isq(): calculates the I-squared statistic, which measures heterogeneity among the genetic instrments
# used in the MR analysis

