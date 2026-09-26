library(DESeq2)

raw_file <- file.choose()

counts <- as.matrix(
  read.csv(raw_file, row.names = 1, check.names = FALSE)
)

# Check that these are valid raw counts.
stopifnot(
  is.numeric(counts),
  !anyNA(counts),
  all(counts >= 0),
  all(counts == floor(counts)),
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  "SOX9-AS1" %in% rownames(counts)
)

storage.mode(counts) <- "integer"

dim(counts)







soft_file <- file.choose()
soft <- readLines(soft_file, warn = FALSE)

starts <- grep("^\\^SAMPLE = ", soft)
ends <- c(starts[-1] - 1L, length(soft))

get_field <- function(block, prefix) {
  value <- block[startsWith(block, prefix)]
  
  if (length(value) != 1L) {
    stop("Missing or duplicated metadata field: ", prefix)
  }
  
  substring(value, nchar(prefix) + 1L)
}

metadata <- do.call(rbind, lapply(seq_along(starts), function(i) {
  block <- soft[starts[i]:ends[i]]
  
  title <- get_field(block, "!Sample_title = ")
  
  data.frame(
    sample = sub(" \\(RNA-seq\\)$", "", title),
    patient = get_field(
      block, "!Sample_characteristics_ch1 = patient: "
    ),
    tissue = get_field(
      block, "!Sample_characteristics_ch1 = tissue: "
    ),
    stringsAsFactors = FALSE
  )
}))

# Match by sample name, never by position.
stopifnot(!anyDuplicated(metadata$sample))

sample_order <- match(colnames(counts), metadata$sample)
stopifnot(!anyNA(sample_order))

metadata <- metadata[sample_order, ]
rownames(metadata) <- metadata$sample

stopifnot(all(
  metadata$tissue %in% c("Lung adjacent normal", "Tumor")
))

metadata$tissue <- factor(
  ifelse(metadata$tissue == "Lung adjacent normal", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)

metadata$patient <- factor(metadata$patient)

# Every patient must have exactly one sample of each tissue.
pair_check <- table(metadata$patient, metadata$tissue)

stopifnot(
  identical(colnames(counts), rownames(metadata)),
  nrow(pair_check) == 123L,
  all(pair_check == 1L)
)

table(metadata$tissue)








dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ patient + tissue
)

# Remove only genes with zero counts in every sample.
dds <- dds[rowSums(counts(dds)) > 0, ]





# 1. Check dimensions first
dim(dds)

# 2. Remove very low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds2 <- dds[keep, ]

dim(dds2)

# 3. Run DESeq2 in parallel
library(BiocParallel)

param <- SnowParam(workers = 4)

dds2 <- DESeq(
  dds2,
  parallel = TRUE,
  BPPARAM = param
)





"SOX9-AS1" %in% rownames(dds2)

design(dds2)

resultsNames(dds2)



res <- results(
  dds2,
  contrast = c("tissue", "Tumor", "Normal"),
  alpha = 0.05
)

sox9_result <- as.data.frame(res["SOX9-AS1", ])

print(sox9_result)






mcols(dds2)["SOX9-AS1", "betaConv"]

warnings()




# Confirm the fitted model has result coefficients
resultsNames(dds2)

# Confirm SOX9-AS1 survived filtering
"SOX9-AS1" %in% rownames(dds2)

# Check convergence for this gene
mcols(dds2)["SOX9-AS1", "betaConv"]

# Extract its result without parallel processing
res <- results(
  dds2,
  contrast = c("tissue", "Tumor", "Normal"),
  alpha = 0.05,
  parallel = FALSE
)

print(as.data.frame(res["SOX9-AS1", ]))







# Create an output folder beside your original raw-count file
output_dir <- file.path(dirname(raw_file), "SOX9_AS1_paired_results")
dir.create(output_dir, showWarnings = FALSE)

# Save the SOX9-AS1 result
write.csv(
  as.data.frame(res["SOX9-AS1", ]),
  file.path(output_dir, "SOX9_AS1_result.csv")
)

# Save results for all tested genes
write.csv(
  as.data.frame(res),
  file.path(output_dir, "all_genes_results.csv")
)

# Save the fitted model for later plots and diagnostic checks
saveRDS(
  dds2,
  file.path(output_dir, "paired_DESeq2_model.rds")
)

# Save sample information and software versions
write.csv(
  as.data.frame(colData(dds2)),
  file.path(output_dir, "sample_metadata.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)

# Show the folder location and saved filenames
normalizePath(output_dir)
list.files(output_dir)