library(dada2)

##########################
# Setup carried over from previous step
log <- function(message) print(paste(date(), message))
samples <- scan("SraAccList.txt", what="character")
forward_reads <- paste0("fastq/", samples, "_1.fastq")
reverse_reads <- paste0("fastq/", samples, "_2.fastq")
filtered_forward_reads <- paste0("intermediate/", samples, ".R1.filtered.fastq.gz")
filtered_reverse_reads <- paste0("intermediate/", samples, ".R2.filtered.fastq.gz")

# This determines whether we should use paired-end processing or single-end
paired <- sum(file.exists(reverse_reads)) == length(reverse_reads)
if(paired) {
    log('Paired-end data found!')
} else {
    log('Processing as single-end data')
}
###########################

# Once filtering is done, limit the list of filtered fastq files to include
# only the ones that actually had reads pass the filter:
filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]
filtered_reverse_reads <- filtered_reverse_reads[file.exists(filtered_reverse_reads)]

#revise the list of samples to only include those
# that actually have reads now:
samples <- gsub('intermediate/(\\w+)\\.R1.filtered.fastq.gz$', '\\1', filtered_forward_reads)

#########################
# Building error models
#########################
log('Building forward error model...')
set.seed(71121)
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE, nbases = 1e8, randomize=TRUE)
pdf('forward_error_model.pdf')
plotErrors(err_forward_reads, nominalQ=TRUE)
dev.off()

if(paired) {
    log('Building reverse error model...')
    set.seed(92124)
    err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE, nbases = 1e8, randomize=TRUE)
    pdf('reverse_error_model.pdf')
    plotErrors(err_reverse_reads, nominalQ=TRUE)
    dev.off()
}

log('Saving error models...')
saveRDS(err_forward_reads, 'err_forward_reads.rds')
if(paired) {
    saveRDS(err_reverse_reads, 'err_reverse_reads.rds')
}

log('Error models saved')
