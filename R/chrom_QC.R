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


path <- paste0("data/chromatin/", analysis_id, "/mapped/")
bamfiles <- paste0(path, list.files(path=path, pattern="*.rmDup.sorted.bam$"))
fileslabels <- gsub(".rmDup.sorted.bam", "", basename(bamfiles))
dir.create(paste0("data/chromatin/", analysis_id, "/QC/fragment_size/"), showWarnings = FALSE, recursive = TRUE)

for (i in 1:length(bamfiles)) {
    print(paste("Plotting", fileslabels[i], "fragment sizes"))
    tiff(paste0("data/chromatin/", analysis_id, "/QC/fragment_size/", fileslabels[i], "_frag_size", ".tiff"), 600, 600)
    ATACseqQC::fragSizeDist(bamfiles[i], fileslabels[i])
    dev.off()
}

tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT") 
seqlvls <- c(1:22, "X") # Insert Y chromosome if necessary !

seq_obj <- GenomeInfoDb::seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
names(seq_obj) <- gsub("chr", "", names(seq_obj))
which <- as(seq_obj[seqlvls], "GRanges")
bpparam <- BiocParallel::MulticoreParam(workers=2)

gal <- list()
for (i in 1:length(bamfiles)) {
    print(paste("Processing", fileslabels[i]))
    gal[[i]] <- ATACseqQC::readBamFile(bamfiles[i], tag=tags,
                            which=which, asMates=T, bigFile=T)
    out_bam <- paste0(fileslabels[i], "_shift.bam")
    objs <- ATACseqQC::shiftGAlignmentsList(gal[[i]], outbam=paste0(path, out_bam))
}


ann <- read.delim("dbs/annotations/gencode_annotation.txt", sep = "\t", header = F)
names(ann)=c("chr", "type", "start", "end","strand", "gene_id")
ann$chr <- gsub("chr", "", ann$chr)

genes <- GenomicRanges::makeGRangesFromDataFrame(ann)


print(paste("Plotting TSS enrichment"))
tsse <- ATACseqQC::TSSEscore(objs, genes)

dir.create(paste0("data/chromatin/", analysis_id, "/QC/TSSEnrichment"), showWarnings = FALSE)
tiff(paste0("data/chromatin/", analysis_id, "/QC/TSSEnrichment/tss_enrichment.tiff"), 600, 600)
  plot(100*(-9:10-.5), tsse$values, type="b", xlab="distance to TSS", ylab="aggregate TSS score")
dev.off()