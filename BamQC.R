BiocManager::install("ATACseqQC")

library(ATACseqQC)

bamfiles <- read.table("bam_files.txt", sep = "\t", header = F, stringsAsFactors = F)

for (i in 1:14){
        current_bamfile = bamfiles[i,]
        cat('Processing file ', i, '\n')
        print(current_bamfile)

	bamQC(current_bamfile)
	
}

