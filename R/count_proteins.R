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

analysis_id <- c(opt$id)

proteinData <- read.delim(paste0("data/proteins/", analysis_id, "/proteinGroups.txt"), sep="\t")

sel <- grep(pattern = "corrected", names(proteinData))
data <- proteinData[, sel]
data_sel <- data[apply(data, 1, sum) > 0, ]

file_out <- paste0("data/proteins/", analysis_id, "/to_DEA.tsv")
write.table(data_sel, file_out, sep="\t", quote=F, row.names=F)
