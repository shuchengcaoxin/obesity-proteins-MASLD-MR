library(vroom)
library(dplyr)
library(magrittr)

# ---- Paths ----
mr_file <- "/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/2.1.pQTL_to_MASLD_GWAS_MVP_new_filtered/output/aptamer_summary.tsv"
protein_file <- "/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/1.3.collect_BMI_driven_proteins/output/all_unique_driven_proteins.tsv"
output_file <- "/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/2.1.pQTL_to_MASLD_GWAS_MVP_new_filtered/output/aptamer_summary_BMI_WHR_BFP.tsv"

# ---- Read data ----
mr_df <- vroom(mr_file, delim = "\t")
proteins <- vroom(protein_file, delim = "\t")$protein

# ---- Filter: keep only rows whose shortname is in the obesity-driven protein set ----
njoin3 <- mr_df %>% 
  filter(shortname %in% proteins)

cat("[INFO] Rows kept after protein filter:", nrow(njoin3), "\n")

# ---- Bonferroni threshold (based on the number of rows retained) and re-annotation of Significance ----
n_tests <- nrow(njoin3)
sig_threshold <- signif(0.05 / n_tests, digits = Inf)
cat("[INFO] Bonferroni threshold (0.05 /", n_tests, ") = ", format(sig_threshold, scientific = TRUE), "\n", sep = "")

# Assumes the p-value column is named pval; change this if your column name differs
njoin4 <- njoin3 %>% 
  mutate(Significance = case_when(
    pval < sig_threshold ~ "Bonferroni",
    pval < 0.05          ~ "Nominal",
    TRUE                 ~ "No"
  ))

# Set the factor order and inspect the counts
njoin4$Significance <- factor(njoin4$Significance, levels = c("No", "Nominal", "Bonferroni"))
print(table(njoin4$Significance))

# ---- Cross-tabulation and single-column checks ----
# Contingency table of Significance by Heterogeneity_test
print(table(njoin4[c("Significance", "Heterogeneity_test")]))

# Counts for Pleiotropy_test
print(table(njoin4$Pleiotropy_test))

# Number passing Bonferroni correction and the heterogeneity test
bonf_pass_het_n <- njoin4 %>% 
  filter(Significance == "Bonferroni", Heterogeneity_test == "Pass") %>% 
  nrow()
cat("[INFO] Bonferroni & Heterogeneity Pass count:", bonf_pass_het_n, "\n")

# Number that additionally pass the pleiotropy test
bonf_pass_both_n <- njoin4 %>% 
  filter(Significance == "Bonferroni", Heterogeneity_test == "Pass", Pleiotropy_test == "Pass") %>% 
  nrow()
cat("[INFO] Bonferroni & (Heterogeneity+Pleiotropy) Pass count:", bonf_pass_both_n, "\n")

# ---- Combined sensitivity test: Pass only if both tests pass ----
njoin4 %<>% mutate(
  Sensitivity_test = case_when(
    Heterogeneity_test == "Pass" & Pleiotropy_test == "Pass" ~ "Pass",
    TRUE ~ "Fail"
  )
)

# Optional: distribution of Sensitivity_test
print(table(njoin4$Sensitivity_test))

# ---- Write results ----
vroom::vroom_write(njoin4, output_file, delim = "\t")
cat("[INFO] Saved filtered & annotated MR results to:\n", output_file, "\n", sep = "")
