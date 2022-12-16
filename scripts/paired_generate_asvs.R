library(dada2)

##########################
# Setup carried over from previous step
log <- function(message) print(paste(date(), message))
samples <- scan("SraAccList.txt", what="character")
forward_reads <- paste0("fastq/", samples, "_1.fastq")
reverse_reads <- paste0("fastq/", samples, "_2.fastq")
filtered_forward_reads <- paste0("intermediate/", samples, ".R1.filtered.fastq.gz")
filtered_reverse_reads <- paste0("intermediate/", samples, ".R2.filtered.fastq.gz")
filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]
filtered_reverse_reads <- filtered_reverse_reads[file.exists(filtered_reverse_reads)]
samples <- gsub('intermediate/(\\w+)\\.R1.filtered.fastq.gz$', '\\1', filtered_forward_reads)
###########################

# load error models from previous step
err_forward_reads <- readRDS('err_forward_reads.rds')
err_reverse_reads <- readRDS('err_reverse_reads.rds')

#########################
# Generate count table
#########################
mergers <- vector("list", length(samples))
names(mergers) <- samples

ddF <- vector("list", length(samples))
names(ddF) <- samples
ddR <- vector("list", length(samples))
names(ddR) <- samples

# TODO: We may be able to fix the "Not all provided files exist" problem
# using this here: https://github.com/benjjneb/dada2/issues/711
for(sam in samples) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(paste("intermediate/", sam, ".R1.filtered.fastq.gz", sep=""))
  ddF[[sam]] <- dada(derepF, err=err_forward_reads, multithread=TRUE)
  derepR <- derepFastq(paste("intermediate/", sam, ".R2.filtered.fastq.gz", sep=""))
  ddR[[sam]] <- dada(derepR, err=err_reverse_reads, multithread=TRUE)
  merger <- mergePairs(ddF[[sam]], derepF, ddR[[sam]], derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

seqtab <- makeSequenceTable(mergers)

# Get rid of really short sequences that can't practically be used
# to assign taxonomy:
seqtab.noshort <- seqtab[,nchar(colnames(seqtab)) > 49]
diff <- length(colnames(seqtab)) - length(colnames(seqtab.noshort))
log(paste('Removed',diff,'ASVs for being too short.'))

# check for chimeras
log('Removing bimeras...')
seqtab.nochim <- removeBimeraDenovo(seqtab.noshort, verbose=T)

#########################
# Check reads dropped at each step
#########################
getN <- function(x) sum(getUniques(x))

print('Calculating summary stats...')
# making a little table
merged_val <- sapply(mergers, getN)
nochim_val <- rowSums(seqtab.nochim)
length_val <- rowSums(seqtab.noshort)
chim_removed_val <- round(((length_val-nochim_val)/merged_val)*100, 1)

summary_tab <- data.frame(dinput=filtered_out[,1],
                          filter=filtered_out[,2], forwd=sapply(ddF, getN),
                          revse=sapply(ddR, getN), merged=merged_val,
                          length=length_val,
                          nonchim=nochim_val,
                          chim_perc=chim_removed_val,
                          retained_perc=round((nochim_val*100)/filtered_out[,1], 1))

# OUTPUT
log('Writing summary output...')
write.table(summary_tab, "summary.tsv",
            sep="\t", quote=F, col.names=NA)
log('Writing ASV table...')
write.table(seqtab.nochim, "ASV.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(seqtab.nochim, 'asv.rds')
log('ASVs recorded.')
