library(getopt)

spec <- matrix(c(
  "id", "i", 1, "character"
  ), byrow = TRUE, ncol = 4
)

opt <- getopt(spec)
analysis_id <- opt$id
metadata <- read.delim(file.path("data/RNA", analysis_id, "metadata.tsv"), sep = "\t")
files <- paste(file.path("data/RNA", analysis_id, "mapped", paste("subset", metadata$Samples, "sorted.bam", sep = "_")))

ann <- read.delim("dbs/annotations/gene_annotations.tsv", sep ="\t")
colnames(ann) <- c("GeneID", "Chr", "Start", "End", "Strand")

print("Read counting on genes")
fc_PE <- Rsubread::featureCounts(files, annot.ext=ann, isPairedEnd=TRUE)

saveRDS(fc_PE, file.path("data/RNA", analysis_id, "rds/fc_out.rds"))

countdata <- fc_PE$counts
colnames(countdata) <- metadata$Samples
write.table(countdata, file.path("data/RNA", analysis_id, "to_DEA.tsv"), sep ="\t", quote = FALSE)