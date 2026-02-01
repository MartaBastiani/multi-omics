library(limma)
library(optparse)

option_list <- list(
    make_option(
        c("-i", "--id"),
        type = "character",
        default = NULL,
        help = "analysis id (required)"
    ),
    make_option(
        c("-o", "--omic"),
        type = "character",
        default = NULL,
        help = "omic type (required)"
    ),
    make_option(
        c("-t", "--p_value_thr"),
        default = 0.05,
        help = "threshold for p-value (default: 0.05)"
    ),
    make_option(
        c("-T", "--log2_fc_thr"),
        default = 1.5,
        help = "threshold for log2FC (default: 0.05)"
    )
)


opt <- parse_args(OptionParser(option_list=option_list))

p_value_thr <- opt$p_value_thr
log2_fc_thr <- opt$log2_fc_thr
analysis_id <- opt$id
omic_type <- opt$omic

metadata <- read.delim(file.path("data/", omic_type, "/", analysis_id, "/metadata.tsv"), sep="\t")
contrast_list <- read.delim(file.path("data/", omic_type, "/", analysis_id, "/contrast.tsv"), sep="\t")
samples <- metadata$Samples
conditions <- metadata$Condition
# batches <- metadata$Batches

filename <- paste0("data/", omic_type, "/", analysis_id, "/to_DEA.tsv")
data <- read.table(filename, sep="\t", header = T)
logData <- log2(data) # log transformation for limma
logData <- normalizeQuantiles(logData)

Rf <- contrast_list[1, ]$Reference
Tg <- contrast_list[1, ]$Target

groups <- factor(conditions)
groups <- relevel(groups, ref=Rf)
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(logData, design)
contrast.matrix <- makeContrasts(paste0(Tg, "-", Rf), levels=design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

res_all <- topTable(fit, adjust="BH", n = Inf)
res_all <- res_all[complete.cases(res_all), ]

colnames(res_all) <- c("log2FoldChange", "AveExpr", "t", "P.Value", "padj", "B")

suppressWarnings(dir.create(paste0("data/", omic_type, "/", analysis_id, "/DEA_results")))
write.table(res_all,
  file = file.path("data", omic_type, analysis_id, "DEA_results",
  paste0(Tg, "_vs_", Rf, "_all.tsv")),
  sep = "\t",
  quote = F,
  row.names = T
)

res <- res_all[abs(res_all$log2FoldChange) >= 1 & res_all$padj < 0.05, ]

write.table(res,
  file = file.path("data", omic_type, analysis_id, "DEA_results",
  paste0(Tg, "_vs_", Rf, "_sel.tsv")),
  sep = "\t",
  quote = F,
  row.names = T
)
