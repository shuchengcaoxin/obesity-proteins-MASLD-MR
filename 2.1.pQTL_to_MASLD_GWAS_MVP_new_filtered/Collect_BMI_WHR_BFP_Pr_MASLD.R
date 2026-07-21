setwd("/scratch/butler/bangli.cao/BMI_driven_proteins_MASLD/3.colocalization/sharepro_NCAN_UKBPPP_eQTL/data/")

suffix <- "_NCAN_UKBPPP_19_19215154_SharePro.txt"
target_snp <- "chr19:19215154:G:T"
p_threshold <- 5e-8  # genome-wide significance

# 获取符合条件的文件名
all_files <- list.files(pattern = paste0(".*", suffix, "$"))
target_files <- all_files[!grepl("MASLD", all_files, fixed = FALSE)]

results <- data.frame(Gene = character(), Tissue = character(), stringsAsFactors = FALSE)

for (f in target_files) {
  cat("Processing:", f, "\n")
  
  if (file.info(f)$size == 0) {
    next
  }
  
  df <- tryCatch({
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    warning("Failed to read: ", f)
    return(NULL)
  })
  
  if (is.null(df) || nrow(df) == 0) next
  if (!("SNP" %in% colnames(df)) || !("P" %in% colnames(df))) {
    warning("Missing SNP or P column in: ", f)
    next
  }
  
  # 条件1: 文件中包含目标 SNP
  has_target_snp <- any(df$SNP == target_snp, na.rm = TRUE)
  
  # 条件2: 文件中存在至少一个 P < 5e-8
  has_significant_p <- any(df$P < p_threshold, na.rm = TRUE)
  
  if (has_target_snp && has_significant_p) {
    prefix <- sub(paste0(suffix, "$"), "", f)
    parts <- strsplit(prefix, "_", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      warning("Unexpected filename format: ", f)
      next
    }
    gene <- parts[1]
    tissue <- paste(parts[-1], collapse = "_")
    results <- rbind(results, data.frame(Gene = gene, Tissue = tissue, stringsAsFactors = FALSE))
  }
}

# 去重并输出
results <- unique(results)
output_file <- "Gene_Tissue_NCAN_UKBPPP.tsv"
write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done. Output written to:", output_file, "\n")
cat("Number of qualifying files:", nrow(results), "\n")