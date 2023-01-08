library(dada2)

log <- function(message) print(paste(date(), message))

## first we're setting a few variables we're going to use
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("SraAccList.txt", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0("fastq/", samples, ".fastq")
# and one with the reverse
reverse_reads <- paste0("fastq/", samples, "_2.fastq")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0("intermediate/", samples, ".R1.filtered.fastq.gz")
filtered_reverse_reads <- paste0("intermediate/", samples, ".R2.filtered.fastq.gz")

# This determines whether we should use paired-end processing or single-end
paired <- sum(file.exists(reverse_reads)) == length(reverse_reads)
if(paired) {
    log('Paired-end data found!')
} else {
    log('Processing as single-end data')
}

#########################
# Quality filtering
#########################

# Filter for quality. THis also does a few DADA2-specific things
# (throwing out 'N' bases, for example) that become important later.
log('Filtering...')

tomatch = FALSE
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  if (args[1] == 'matchids') {
    tomatch = TRUE
  }
}

if(paired) {
    filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              truncQ=2, rm.phix=TRUE, multithread=8,
                              verbose=TRUE, matchIDs=tomatch)
} else {
    filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              truncQ=2, rm.phix=TRUE, multithread=8,
                              verbose=TRUE)
}

log('Filtering complete. Saving results...')

saveRDS(filtered_out, file='filtered_out.rds')

log('Filtering results saved. Done.')
