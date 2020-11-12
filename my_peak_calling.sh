#Peak calling on the cleaned (ATACseqQC) and sorted (samtools) bam files (HMMRATAC)
for filename in *.clean.sorted.bam
do
        base=$(basename $filename .clean.sorted.bam)
        echo $base

        java -jar /home/jovyan/HMMRATAC/HMMRATAC_V1.2.10_exe.jar -b ${base}.clean.sorted.bam  -i ${base}.clean.sorted.bam.bai -g ${base}.clean.sorted.genome.info --bedgraph True -o ${base}.cleaned_HMM_bedgraph_out

done
