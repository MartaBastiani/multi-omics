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
data <- proteinData[proteinData$Gene.names != "", ] # removing proteins without a gene symbol
data <- data[!grepl(";", data$Gene.names), ] # important: removing proteinGroups with multiple genes identified.

sel <- grep(pattern = "corrected", names(proteinData)) 
data <- aggregate(data[, sel], by = list(Gene.names = data$Gene.names), FUN = sum) # sum of identical genes
rownames(data) <- data$Gene.names
data$Gene.names <- NULL 

countdata <- data[apply(data, 1, function(x) !any(x == 0)), ] # remove rows with 0s

colnames(countdata) <- gsub("Reporter.intensity.corrected.", "", colnames(countdata))

file_out <- paste0("data/proteins/", analysis_id, "/to_DEA.tsv")
write.table(countdata, file_out, sep="\t", quote=F, row.names=T)