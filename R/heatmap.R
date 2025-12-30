library(getopt)

spec <- matrix(c(
    "verbose", "v", 2, "integer",
    "id", "i", 1, "character",
    "omic", "o", 1, "character",
    "reference", "R", 2, "character",
    "target", "T", 2, "character"
), byrow=TRUE, ncol=4)

opt <- getopt(spec)

analysis_id <- opt$id
omic_type <- opt$omic


metadata <- read.delim(file.path("data", omic_type, analysis_id, "metadata.tsv"), sep="\t")

if (omic_type == "proteins"){
  n_counts <- read.delim(file.path("data/proteins", analysis_id, "to_DEA.tsv"), sep="\t")
} else {
  n_counts <- readRDS(file.path("data", omic_type, analysis_id, "rds", "deseq2_data_counts_norm.rds"))
}
z_n_counts <- t(scale(t(n_counts)))

dir.create(file.path("data", omic_type, analysis_id, "plots/heatmap"), recursive = TRUE, showWarnings = FALSE)

for (rn in 1:length(rownames(contrast_list))) {
  reference <- contrast_list[rn, ]$Reference  
  target <- contrast_list[rn, ]$Target

  ann <- metadata[metadata$Condition %in% c(target, reference), c(1, 2)]
  rownames(ann) <- ann[, 1]
  ann <- ann[, 2, drop = FALSE]
  data <- z_n_counts[, rownames(ann)]
  ann_col <- list(Condition = setNames(c("#bb0c00", "#00AFBB"), c(target, reference)))

  png(file.path("data", omic_type, analysis_id, "plots/heatmap", paste0(target, "_vs_", reference, "_heatmap.png")), 1500, 1500, res = 300)
  pheatmap::pheatmap(
    mat = data,
    show_rownames = FALSE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    annotation_legend = TRUE,
    annotation_col = ann,
    show_colnames = FALSE,
    col = RColorBrewer::brewer.pal(10, "RdYlGn"),
    treeheight_row = 30,
    treeheight_col = 10,
    width = 7.5,
    height = 4,
    annotation_colors = ann_col
  )
  dev.off()
}


# Pink and green color palette
# "PiYG"
