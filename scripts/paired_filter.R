library(dada2)

log <- function(message) print(paste(date(), message))

## first we're setting a few variables we're going to use
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("SraAccList.txt", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0("fastq/", samples, "_1.fastq")
# and one with the reverse
reverse_reads <- paste0("fastq/", samples, "_2.fastq")
# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0("intermediate/", samples, ".R1.filtered.fastq.gz")
filtered_reverse_reads <- paste0("intermediate/", samples, ".R2.filtered.fastq.gz")


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
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              truncQ=2, rm.phix=TRUE, multithread=8,
                              verbose=TRUE, matchIDs=tomatch)

log('Filtering complete')
