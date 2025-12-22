library(getopt)

spec <- matrix(c(
  "id", "i", 1, "character"
  ), byrow = TRUE, ncol = 4
)


files <- read.delim(pipe("ls data/RNA/", analysis_id, "mapped/ | grep _sorted.bam"), header =F)[,1]
ann <- read.delim("dbs/annotations/insert_name")
metadata <- read.delim(file.path("data/RNA", analysis_id, "metadata.tsv"), sep="\t")

fc_PE <- featureCounts(paste("mapped/", files, sep = ""), annot.ext=ann, isPairedEnd=TRUE)

saveRDS(fc_PE, file.path("data/RNA", analysis_id, "rds/count_RNA_out.rds"))

countdata <- fc_PE$counts
colnames(countdata) <- metadata$Samples
