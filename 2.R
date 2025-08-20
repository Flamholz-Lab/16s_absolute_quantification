library(dada2); packageVersion("dada2")

# change path to be your folder containing trimmed reads
path <- "/Users/albertli/Desktop/lab/Flamholz/projects/absolute_16s_quantification/data/standard_sequencing_07222025/trimmed"
fnFs <- sort(list.files(path, pattern="_R1_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_trimmed.fastq", full.names = TRUE))
plotQualityProfile(fnFs[1:14]) # change based on how many samples you have
plotQualityProfile(fnRs[1:14])

samples <- scan(file.path(path, "samples"), what = "character")
filtFs <- file.path(path, "filtered", paste0(samples, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(samples, "_R_filt.fastq"))
names(filtFs) <- samples
names(filtRs) <- samples

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,210), maxN=0, maxEE=c(4,5), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
print(out)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- samples
print(track)

# change to location of your database
taxa <- assignTaxonomy(seqtab.nochim, "/Users/albertli/Desktop/lab/Flamholz/projects/absolute_16s_quantification/data/dada2_taxomony/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(Biostrings); packageVersion("Biostrings")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Subject", fill="Family") + labs(y = "Fractional Abundance")

count_table <- as.data.frame(otu_table(ps))
if (taxa_are_rows(ps)) {
  count_table <- t(count_table)
}
print(count_table)
write.csv(count_table, "../outputs/ASV_counts_per_sample.csv", row.names = TRUE)

tax_table_df <- as.data.frame(tax_table(ps))
tax_table_df$ASV <- rownames(tax_table_df)
print(tax_table_df)
write.csv(tax_table_df, "../outputs/ASV_to_taxonomy.csv", row.names = TRUE)
