# library(methylKit)
library(optparse)

option_list <- list(
    make_option(
        c("-i", "--id"),
        type = "character",
        default = NULL,
        help = "analysis id (required)"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
analysis_id <- opt$id

path <- file.path("data/methylations", analysis_id, "/")
metadata <- read.delim(paste0(path, "metadata.tsv"), sep = "\t")
files <- paste0(path, "mapped/", metadata$Samples, "_sorted_pe.bam")

samples <- metadata$Samples


dir.create(paste0(path, "CpG_calls"), showWarnings = FALSE)
data=list()

for (i in 1:length(files)) {
	cat("Processing", samples[i], "\n")
	data[[samples[i]]]=methylKit::processBismarkAln(files[i], samples[i], "hg19", nolap = TRUE, read.context = "CpG", save.folder = paste0(path, "CpG_calls"))
}

# save(file=paste0(path, "methyl_data.RData"), data)

