library(vroom)
library(dplyr)
library(magrittr)

# ---- 路径设定 ----
mr_file <- "/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/2.1.pQTL_to_MASLD_GWAS_MVP_new_filtered/output/aptamer_summary.tsv"
protein_file <- "/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/1.3.collect_BMI_driven_proteins/output/all_unique_driven_proteins.tsv"
output_file <- "/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/2.1.pQTL_to_MASLD_GWAS_MVP_new_filtered/output/aptamer_summary_BMI_WHR_BFP.tsv"

# ---- 读取数据 ----
mr_df <- vroom(mr_file, delim = "\t")
proteins <- vroom(protein_file, delim = "\t")$protein

# ---- 过滤：只保留 shortname 在并集蛋白内的行 ----
njoin3 <- mr_df %>% 
  filter(shortname %in% proteins)

cat("[INFO] Rows kept after protein filter:", nrow(njoin3), "\n")

# ---- 计算 Bonferroni 阈值（基于筛完后的行数）并重写 Significance ----
n_tests <- nrow(njoin3)
sig_threshold <- signif(0.05 / n_tests, digits = Inf)
cat("[INFO] Bonferroni threshold (0.05 /", n_tests, ") = ", format(sig_threshold, scientific = TRUE), "\n", sep = "")

# 假设 p 值列名为 pval（若不同请改为你的列名）
njoin4 <- njoin3 %>% 
  mutate(Significance = case_when(
    pval < sig_threshold ~ "Bonferroni",
    pval < 0.05          ~ "Nominal",
    TRUE                 ~ "No"
  ))

# 设定因子顺序并查看频数
njoin4$Significance <- factor(njoin4$Significance, levels = c("No", "Nominal", "Bonferroni"))
print(table(njoin4$Significance))

# ---- 联合/单列检视 ----
# 按 Significance × Heterogeneity_test 的列联表
print(table(njoin4[c("Significance", "Heterogeneity_test")]))

# Pleiotropy_test 的频数
print(table(njoin4$Pleiotropy_test))

# 满足 Bonferroni 且异质性通过的数量
bonf_pass_het_n <- njoin4 %>% 
  filter(Significance == "Bonferroni", Heterogeneity_test == "Pass") %>% 
  nrow()
cat("[INFO] Bonferroni & Heterogeneity Pass count:", bonf_pass_het_n, "\n")

# 同时 pleiotropy 也通过的数量
bonf_pass_both_n <- njoin4 %>% 
  filter(Significance == "Bonferroni", Heterogeneity_test == "Pass", Pleiotropy_test == "Pass") %>% 
  nrow()
cat("[INFO] Bonferroni & (Heterogeneity+Pleiotropy) Pass count:", bonf_pass_both_n, "\n")

# ---- 汇总敏感性测试：两者都Pass才算Pass ----
njoin4 %<>% mutate(
  Sensitivity_test = case_when(
    Heterogeneity_test == "Pass" & Pleiotropy_test == "Pass" ~ "Pass",
    TRUE ~ "Fail"
  )
)

# 可选：看一下 Sensitivity_test 的分布
print(table(njoin4$Sensitivity_test))

# ---- 写出结果 ----
vroom::vroom_write(njoin4, output_file, delim = "\t")
cat("✅ Saved filtered & annotated MR results to:\n", output_file, "\n", sep = "")
