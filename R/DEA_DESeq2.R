library(DESeq2)
library(getopt)


spec = matrix(c(
  'verbose', 'v', 2, "integer", # logical, integer, double, complex, character.
  'help'   , 'h', 0, "logical",
  'id'  , 'i', 1, "character",  
  "omic", "o", 1, "character" # 0=no argument, 1=required argument, 2=optional argument
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# p_value_thr<-opt$p_value_thr
# log2_fc_thr<-opt$log2_fc_thr
p_value_thr <-0.05
log2_fc_thr <- 1.5
analysis_id <- opt$id
omic_type <- opt$omic



metadata <- read.delim(file.path("data", omic_type, analysis_id, "metadata.tsv"), sep="\t")
contrast_list <- read.delim(file.path("data", omic_type, analysis_id, "contrast.tsv"), sep="\t")
samples <- metadata$Samples
conditions <- metadata$Condition
references <- unique(contrast_list$Reference)
batches <- metadata$Batches

# txi <- readRDS(file.path("data", omic_type, analysis_id, "rds", "txi_salmon.rds"))
count_data <- round(read.delim(file.path("data", omic_type, analysis_id, "to_DEA.tsv"), sep="\t"), 0)



condition <- factor(
	conditions,
	levels = unique(conditions)
)

#iterate over batches here



sampleTable <- data.frame(
  sample = samples,
  condition = conditions
)


previous_reference <- contrast_list[1, ]$Reference
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sampleTable, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 0, ]
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds$condition <- relevel(dds$condition, ref = contrast_list[1,]$Reference)
dds <- DESeq(dds)

dir.create(file.path("data", omic_type, analysis_id, "rds"), showWarnings = FALSE)
saveRDS(dds, file = file.path("data", omic_type, analysis_id, "rds", paste0("deseq2_data","_vs_", contrast_list[1,]$Reference, ".rds")))
saveRDS(counts(dds), file = file.path("data", omic_type, analysis_id, "rds","deseq2_data_counts.rds"))
saveRDS(counts(dds, normalized = TRUE), file = file.path("data", omic_type, analysis_id, "rds","deseq2_data_counts_norm.rds"))




for (rn in 1:length(rownames(contrast_list))) {

  row <- contrast_list[rn,]

  if (row$Reference != previous_reference) {
    dds$condition <- relevel(dds$condition, ref = row$Reference)
    dds <- DESeq(dds)
    previous_reference <- row$Reference
    saveRDS(dds, file = file.path("data", omic_type, analysis_id,"rds",paste0("deseq2_data","_vs_",row$Reference,".rds")))
  }

  filename1 <- paste(row$Target, row$Reference, sep="_vs_")
  filename2 <- paste(filename1, ".tsv", sep="")
  filename2b <- paste(filename1, "_filtered", ".tsv", sep="")
  filename2c <- paste(filename1, "_signedP", ".tsv", sep="")
  res <- results(dds,contrast=c(row$Type, row$Target, row$Reference))
  res_ordered <- res[order(res$padj), ]
  res_signed <- data.frame(
    gene_id=rownames(res_ordered),
    log2FoldChange=res_ordered$log2FoldChange,
    padj=res_ordered$padj
  )
  res_signed$signedP <- -log10(res_signed$padj) * sign(res_signed$log2FoldChange)
  res_signed <- res_signed[complete.cases(res_signed),]
  res_signed <- res_signed[order(res_signed$signedP,decreasing = TRUE),]

  res_subset <- subset(res_ordered, padj < p_value_thr & abs(log2FoldChange) > log2_fc_thr)
  dir.create(file.path("data", omic_type, analysis_id, "DEA_results"), showWarnings = FALSE)
  write.table(as.data.frame(res_ordered),file.path("data", omic_type, analysis_id,"DEA_results",filename2),sep="\t",row.names=TRUE,quote=FALSE)
  write.table(as.data.frame(res_subset),file.path("data", omic_type, analysis_id,"DEA_results",filename2b),sep="\t",row.names=TRUE,quote=FALSE)
  write.table(as.data.frame(res_signed),file.path("data", omic_type, analysis_id,"DEA_results",filename2c),sep="\t",row.names=FALSE,quote=FALSE)

  coef <- paste0(row$Type,"_",row$Target,"_vs_",row$Reference)
  resLFC <- lfcShrink(dds,coef=coef, type="apeglm")
  resLFC_ordered <- resLFC[order(resLFC$padj),]
  resLFC_subset <- subset(resLFC_ordered, padj < p_value_thr & abs(log2FoldChange) > log2_fc_thr)
  filename3<-paste0(filename1,"_shrinked",".tsv")
  filename3b<-paste0(filename1,"_shrinked_filtered",".tsv")

  write.table(as.data.frame(resLFC_ordered),file.path("data", omic_type, analysis_id,"DEA_results",filename3),sep="\t",row.names=TRUE,quote=FALSE)
  write.table(as.data.frame(resLFC_subset),file.path("data", omic_type, analysis_id,"DEA_results",filename3b),sep="\t",row.names=TRUE,quote=FALSE)

}









# for (rn in 1:length(rownames(contrast_list))) {
#   row <- contrast_list[rn,]
#   filename1<-paste(row$Target,row$Reference,sep="_vs_")
#   filename2<-paste(filename1,".tsv",sep="")
#   filename2b<-paste(filename1,"_filtered",".tsv",sep="")
#   filename2c<-paste(filename1,"_signedP",".tsv",sep="")
#   res <- results(dds,contrast=c(row$Type,row$Target,row$Reference))
#   res_ordered<-res[order(res$padj),]
#   res_signed<-data.frame(
#     gene_id=rownames(res_ordered),
#     log2FoldChange=res_ordered$log2FoldChange,
#     padj=res_ordered$padj
#   )

#   res_signed$signedP <- -log10(res_signed$padj) * sign(res_signed$log2FoldChange)
#   res_signed<-res_signed[complete.cases(res_signed),]
#   res_signed<-res_signed[order(res_signed$signedP,decreasing = TRUE),]

#   res_subset <- subset(res_ordered, padj < p_value_thr & abs(log2FoldChange) > log2_fc_thr)
#   write.table(as.data.frame(res_ordered),file.path("analysis",analysis_id,"tsv",filename2),sep="\t",row.names=TRUE,quote=FALSE)
#   write.table(as.data.frame(res_subset),file.path("analysis",analysis_id,"tsv",filename2b),sep="\t",row.names=TRUE,quote=FALSE)
#   write.table(as.data.frame(res_signed),file.path("analysis",analysis_id,"tsv",filename2c),sep="\t",row.names=FALSE,quote=FALSE)

#   if (row$Reference == conditions[1]) {
#     coef<-paste0(row$Type,"_",row$Target,"_vs_",row$Reference)
#     resLFC <- lfcShrink(dds,coef=coef, type="apeglm")
#     resLFC_ordered <- resLFC[order(resLFC$padj),]
#     resLFC_subset <- subset(resLFC_ordered, padj < p_value_thr & abs(log2FoldChange) > log2_fc_thr)
#     filename3<-paste0(filename1,"_shrinked",".tsv")
#     filename3b<-paste0(filename1,"_shrinked_filtered",".tsv")

#     write.table(as.data.frame(resLFC_ordered),file.path("analysis",analysis_id,"tsv",filename3),sep="\t",row.names=TRUE,quote=FALSE)
#     write.table(as.data.frame(resLFC_subset),file.path("analysis",analysis_id,"tsv",filename3b),sep="\t",row.names=TRUE,quote=FALSE)

#   }

# }
