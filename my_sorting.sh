#Sorting script (Files have to be sorted for peak-calling with HMMRATAC)
for filename in *_1.qc.fq.gz
do
        base=$(basename $filename .qc.fq.gz)
        echo $base

        #/storage/atac_data/original_data/wBush_RNA_ATAC_fastqs/samtools-1.9/samtools sort -@ 16 -o ${base}.sorted.bam  ${base}.bam  

	#/storage/atac_data/original_data/wBush_RNA_ATAC_fastqs/samtools-1.9/samtools index ${base}.sorted.bam -@ 16
	
	/storage/atac_data/original_data/wBush_RNA_ATAC_fastqs/samtools-1.9/samtools view -H ${base}.sorted.bam -@ 16| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > ${base}.genome.info


done
