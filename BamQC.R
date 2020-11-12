# This is a script to run the "bamQC" function on our bam files. The function removes mitochondrial mapped reads, PCR amplification artifacts, and low-quality mappings.

BiocManager::install("ATACseqQC")

library(ATACseqQC)

bamfiles <- read.table("bam_files.txt", sep = "\t", header = F, stringsAsFactors = F)

for (i in 1:14){
        current_bamfile = bamfiles[i,]
        cat('Processing file ', i, '\n')
        print(current_bamfile)

	bamQC(current_bamfile)
	
}

