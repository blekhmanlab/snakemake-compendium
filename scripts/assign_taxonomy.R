library(dada2)

##########################
# Setup carried over from previous step
log <- function(message) print(paste(date(), message))
###########################

# load asv table from previous step
seqtab.nochim <- readRDS('asv.rds')
log('Assigning taxonomy...')
taxa <- assignTaxonomy(seqtab.nochim, "resources/silva_nr99_v138_train_set.fa.gz", multithread=8, tryRC=T)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# extract fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)

## output
log('Writing output files...')
write(asv_fasta, "ASVs.fa")
write.table(asv_tab, "ASVs_counts.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax, "ASVs_taxonomy.tsv",
            sep="\t", quote=F, col.names=NA)

log('DONE!!!')
