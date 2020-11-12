#Script for enrichment of certain MS motifs of interest in the peaks of our second ATACseq dataset.
setwd("~/Desktop/Case_PREP/Bush_Lab_Work/ATACseq_MS_PR/peak_calls/gappedPeak/")


library(TFBSTools)
library(BSgenome)
library(tidyverse)
library(rtracklayer)
library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(MotifDb)
library(ChIPseeker)
library(TCseq)



my_dir <- "~/Desktop/Case_PREP/Bush_Lab_Work/ATACseq_MS_PR/peak_calls/gappedPeak/"
myfiles <- list.files(path = my_dir, pattern = "*.gappedPeak", full.names = T)

# This reads in and merges all of the peaks files into a single data frame
peak_merges <- plyr::ldply(myfiles, read.table)

peak_merges_subset <- peak_merges[, c("V1", "V2", "V3", "V4", "V13")]

# This removes the rows that had start coordinates with dashes in them (example: 1445-1447, instead of a position). Solving this allows me to use GRanges.
peak_merges_subset <- peak_merges_subset %>% filter(!str_detect(peak_merges_subset$V2, pattern = "-"))

peak_merges_subset <- peak_merges_subset[order(peak_merges_subset$V13, decreasing = T),]

peak_merges_subset$V2 <- as.numeric(peak_merges_subset$V2) # Cast these columns to numeric types. They were characters.
peak_merges_subset$V3 <- as.numeric(peak_merges_subset$V3)


overlapped_peak_merges_subset <- peakreference(peak_merges_subset) # Merge peaks that have genomic overlap

gr1 <- makeGRangesFromDataFrame(overlapped_peak_merges_subset, keep.extra.columns=TRUE,
                                start.field="start",
                                end.field="end",
                                seqnames.field = "chr",
                                starts.in.df.are.0based=TRUE)

seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1)

writeXStringSet(seq, filepath = "overlapped_peak_merges_sequences.fasta", format = "fasta")


# CTCF motif 
motifs = query(query(MotifDb, 'Hsapiens'), 'CTCF')

ctcf_motif  = motifs[[1]]

ctcf_pwm = PWMatrix(
  ID = 'CTCF', 
  profileMatrix = ctcf_motif
)

hits = searchSeq(ctcf_pwm, seq, min.score="80%", strand="*")
hits = as.data.frame(hits)


motif_hits_df = data.frame(
  peak_order     = 1:length(gr1)
)
motif_hits_df$contains_motif = motif_hits_df$peak_order %in% hits$seqnames
motif_hits_df = motif_hits_df[order(-motif_hits_df$peak_order),]

# calculate the percentage of peaks with motif for peaks of descending strength
motif_hits_df$perc_peaks = with(motif_hits_df, 
                                cumsum(contains_motif) / max(peak_order))
motif_hits_df$perc_peaks = round(motif_hits_df$perc_peaks, 2)

ggplot(
  motif_hits_df, 
  aes(
    x = peak_order, 
    y = perc_peaks
  )) +
  geom_line(size=2) +
  theme_minimal() +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=14),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Peak rank') +
  ylab('Percetage of peaks with motif') +
  ggtitle('Percentage of peaks with the CTCF motif')

# The rest of this script is dedicated to surveying the motifs of 6 TFs that are of relevance in MS (according to an article
# by Parnell et al.) I'll do the same as I did for the CTCF motif.

# EOMES motif
motifs = query(query(MotifDb, 'Hsapiens'), 'eomes')

eomes_motif  = motifs[[1]]

eomes_pwm = PWMatrix(
  ID = 'EOMES', 
  profileMatrix = eomes_motif
)

hits = searchSeq(eomes_pwm, seq, min.score="90%", strand="*")
hits = as.data.frame(hits)


motif_hits_df = data.frame(
  peak_order     = 1:length(gr1)
)
motif_hits_df$contains_motif = motif_hits_df$peak_order %in% hits$seqnames
motif_hits_df = motif_hits_df[order(-motif_hits_df$peak_order),]

# calculate the percentage of peaks with motif for peaks of descending strength
motif_hits_df$perc_peaks = with(motif_hits_df, 
                                cumsum(contains_motif) / max(peak_order))
motif_hits_df$perc_peaks = round(motif_hits_df$perc_peaks, 2)

ggplot(
  motif_hits_df, 
  aes(
    x = peak_order, 
    y = perc_peaks
  )) +
  geom_line(size=2) +
  theme_minimal() +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=14),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Peak rank') +
  ylab('Percetage of peaks with motif') +
  ggtitle('Percentage of peaks with the EOMES motif') # Around 2% of all peaks contain the EOMES motif


# TBX21 motif
motifs = query(query(MotifDb, 'Hsapiens'), 'tbx21')

tbx21_motif  = motifs[[1]]

tbx21_pwm = PWMatrix(
  ID = 'TBX21', 
  profileMatrix = tbx21_motif
)

hits = searchSeq(tbx21_pwm, seq, min.score="90%", strand="*")
hits = as.data.frame(hits)


motif_hits_df = data.frame(
  peak_order     = 1:length(gr1)
)
motif_hits_df$contains_motif = motif_hits_df$peak_order %in% hits$seqnames
motif_hits_df = motif_hits_df[order(-motif_hits_df$peak_order),]

# calculate the percentage of peaks with motif for peaks of descending strength
motif_hits_df$perc_peaks = with(motif_hits_df, 
                                cumsum(contains_motif) / max(peak_order))
motif_hits_df$perc_peaks = round(motif_hits_df$perc_peaks, 2)

ggplot(
  motif_hits_df, 
  aes(
    x = peak_order, 
    y = perc_peaks
  )) +
  geom_line(size=2) +
  theme_minimal() +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=14),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Peak rank') +
  ylab('Percetage of peaks with motif') +
  ggtitle('Percentage of peaks with the TBX21 motif') # Around 3% of all peaks contain the EOMES motif


