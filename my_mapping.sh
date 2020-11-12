#Mapping script
for filename in *_1.qc.fq.gz
do
  	base=$(basename $filename .qc.fq.gz)
        echo $base

        baseR2=${base/_1/_2}
        echo $baseR2
	/home/jovyan/bwa/bwa mem -t 16 /storage/atac_data/original_data/wBush_RNA_ATAC_fastqs/ATACseq_set2/hg_38/hg38.fa.gz ${base}.qc.fq.gz ${baseR2}.qc.fq.gz | /storage/atac_data/original_data/wBush_RNA_ATAC_fastqs/samtools-1.9/samtools sort -@16 -o ${base}.bam
	
	/storage/atac_data/original_data/wBush_RNA_ATAC_fastqs/samtools-1.9/samtools index ${base}.bam	

done
