library(optparse)

option_list <- list(
    make_option(
        c("-i", "--id"),
        type = "character",
        default = NULL,
        help = "analysis id (required)"
    ),
    make_option(
        c("-t", "--q_value_thr"),
        default = 0.05,
        help = "threshold for q-value (default: 0.05)"
    ),
    make_option(
        c("-T", "--meth_diff_thr"),
        default = 1.5,
        help = "threshold for methylation difference (default: 25%)"
    ),
    make_option(
        c("-R", "--ranges"),
        default = NULL,
        help = "region to perform the DEA on (1 for promoter, 2 for tss)"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$q_value_thr)) {
    q_value_thr <- 0.05
} else {
    q_value_thr <- opt$q_value_thr
}

if (is.null(opt$meth_diff_thr)) {
    meth_diff_thr <- 25
} else {
    meth_diff_thr <- opt$meth_diff_thr
}

analysis_id <- opt$id

metadata <- read.delim(file.path("data/methylations", analysis_id, "metadata.tsv"), sep = "\t")
contrast_list <- read.delim(file.path("data/methylations", analysis_id, "contrast.tsv"), sep="\t")
ann <- read.delim("dbs/annotations/gencode_annotation.txt", sep = "\t", header = F)
names(ann)=c("chr", "type", "start", "end","strand", "gene_id")

genes <- GenomicRanges::makeGRangesFromDataFrame(ann)
# adding the metadata column
GenomicRanges::values(genes) <- data.frame(gene_id = ann$gene_id)

# the promoter function considers also the strand (- or +) in calculating the ranges
promoters <- IRanges::promoters(genes, upstream = 1000, downstream = 0)
# tss <- IRanges::promoters(genes, upstream = 50, downstream = 50)


diffmet = function(indata, ann) {
	# median centering of the coverage vector
	cat(" Normalizing...\n")
	data.norm = methylKit::normalizeCoverage(indata)
    united = methylKit::unite(data.norm, destrand=FALSE)
    cat(" Regionalizing...\n")
    ready = methylKit::regionCounts(united, ann)
	# compute differential statistics
	cat(" Calculating diffs...\n")
	diff = methylKit::calculateDiffMeth(ready)
	# filters results according to given thresholds
	cat(" Filtering...\n")
    # annotate the row results with data from the GRange object
	cat(" Annotating...\n")
	ann2 = as(ann,"data.frame")
	ann2 = data.frame("range"=paste0(ann2[, 2], "-", ann2[, 3]), "name"=ann2[, 6])
	# builds a data frame with annotated results
	tmp = as(GenomicRanges::ranges(as(diff, "GRanges")),"data.frame")[, 1:2]
	# constructs the indexing labels as start-end of granges
	sel = paste0(tmp[, 1], "-", tmp[, 2])
	# associates annotations
	ann_sel = ann2[match(sel, ann2[, 1]), 2]
	# creates the final output data.frame
	out = cbind(diff, ann_sel)
	out
}

dir.create(file.path("data/methylations", analysis_id, "DEA_results"), showWarnings = FALSE)
for (i in 1:length(rownames(contrast_list))) {
  reference <- contrast_list[i, "Reference"]
  target <- contrast_list[i, "Target"]
  print(paste("Processing", target, "and", reference))
  files <- paste(file.path("data/methylations", analysis_id, "mapped",
    paste0(c(target, reference), ".CPG.txt")))
  file_list <- as.list(files)
  meths <- methylKit::methRead(file_list,
    sample.id = list(target, reference), assembly = "hg19", treatment = c(1, 0),
    context = "CpG", mincov = 10, pipeline = "bismark")
  print(paste("Calculating differential methylations between", target, "and", reference))
  diff <- diffmet(meths, promoters)
  colnames(diff)[c(6, 7)] <- c("significance", "diff_amount")
  write.table(diff, file=file.path("data/methylations", analysis_id, "DEA_results", paste0(target, "_vs_", reference, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  diff_sel <- diff[diff$significance < q_value_thr & abs(diff$diff_amount) > meth_diff_thr, ]
  write.table(diff_sel, file=file.path("data/methylations", analysis_id, "DEA_results", paste0(target, "_vs_", reference, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
