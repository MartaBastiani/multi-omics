# Proteins
Rscript R/count_proteins.R -i prova 2> count_proteins.log
Rscript R/DEA_limma.R -i prova -o proteins 2> DEA_limma_pro.log
Rscript R/volcano_plot.R -i prova -o proteins 2> volcano_plot_pro.log


# methylations
# genome preparation
# bismark_genome_preparation --verbose dbs/genomes/
Rscript bash/meth_map.sh -i prova 2> meth_map.log
Rscript R/meth_import_alignment.R -i prova 2> meth_import.log
Rscript R/meth_DEA.R -i prova 2> meth_DEA.log


# miRNA
Rscript R/count_miRNA.R -i prova 2> count_miRNA.log
Rscript R/DEA_limma.R -i prova -o miRNA 2> DEA_limma_miRNA.log
Rscript R/volcano_plot.R -i prova -o miRNA 2> volcano_plot_miRNA.log
