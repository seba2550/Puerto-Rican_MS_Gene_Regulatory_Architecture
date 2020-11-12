# Script for merging all of the peak files for the second set of ATAC-seq data, merging overlapping peaks, and doing some exploratory data analyses on the resulting peaks.

setwd("~/Desktop/Case PREP/Bush_Lab_Work/ATACseq_MS_PR/peak_calls/gappedPeak/")

library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("clusterProfiler")
library("tidyverse")
library("org.Hs.eg.db")
library("ggupset")
library("ggimage")
library("ReactomePA")
library("DOSE")
library("enrichplot")
library("TCseq")


my_dir <- "~/Desktop/Case PREP/Bush_Lab_Work/ATACseq_MS_PR/peak_calls/gappedPeak/"
myfiles <- list.files(path = my_dir, pattern = "*.gappedPeak", full.names = T)

# This reads in and merges all of the peaks files into a single data frame
peak_merges <- plyr::ldply(myfiles, read.table)

peak_merges_subset <- peak_merges[, c("V1", "V2", "V3", "V4", "V13")]

# This removes the rows that had start coordinates with dashes in them (example: 1445-1447, instead of a position). Solving this allows me to use GRanges.
peak_merges_subset <- peak_merges_subset %>% filter(!str_detect(peak_merges_subset$V2, pattern = "-"))

peak_merges_subset$V2 <- as.numeric(peak_merges_subset$V2) # Cast these columns to numeric types. They were characters.
peak_merges_subset$V3 <- as.numeric(peak_merges_subset$V3)


overlapped_peak_merges_subset <- peakreference(peak_merges_subset) # Merge peaks that have genomic overlap
write.table(overlapped_peak_merges_subset, file = "overlapped_merged_peaks_ATACseq_set2.txt",
            row.names = F, quote = F, col.names = F)
write_rds(overlapped_peak_merges_subset, file = "overlapped_merged_peaks_ATACseq_set2.rds")


gr1 <- makeGRangesFromDataFrame(overlapped_peak_merges_subset, keep.extra.columns=TRUE,
                                start.field="start",
                                end.field="end",
                                seqnames.field = "chr",
                                starts.in.df.are.0based=TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(gr1, windows = promoter)
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = "red", xlab = "Distance From TSS (kb)", ylab = "Peak Tags")
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5' -> 3')", ylab = "Read Count Frequency",
            conf = 0.95, resample = 500)

peakAnno <- annotatePeak(gr1, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie = T)

plotDistToTSS(peakAnno,
              title="Distribution of accessible chromatin peaks \nrelative to TSS")
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)

peakAnno_df <- as.data.frame(peakAnno)
write.table(peakAnno_df, file = "merged_peaks_annotated.txt", row.names = T, quote = F)
saveRDS(peakAnno_df, file = "peakAnno_df.rds")

gene <- seq2gene(gr1, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = txdb)
pathway2 <- enrichPathway(gene)
dotplot(pathway2)

# Enriching for disease-related genes
edo <- enrichDGN(gene)
barplot(edo, showCategory = 50)

# Now with dot plots
p1 <- dotplot(edo, showCategory = 30)
p1

# Networks
edox <- setReadable(edo, "org.Hs.eg.db")
p2 <- cnetplot(edox, showCategory = 2) # Way too dense. I'll try this out later on with less genes
p2


# Heatmap
p3 <- heatplot(edox) # Same as the networks plot
p3


# Emap
p4 <- emapplot(edo, layout = "nicely")
p4

# GO
ego <- enrichGO(gene = gene,
                universe = names(gene),
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = T)

goplot(ego)

# KEGG
kk <- enrichKEGG(gene = gene,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
head(kk)



