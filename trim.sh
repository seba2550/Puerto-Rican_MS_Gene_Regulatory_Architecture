#Script for trimming adapters
for filename in *_1.txt.bz2
do
	base=$(basename $filename .txt.bz2)
	echo $base
	
	baseR2=${base/_1/_2}
	echo $baseR2


	java -jar /home/jovyan/bin/trimmomatic.jar PE ${base}.txt.bz2 ${baseR2}.txt.bz2  ${base}.qc.fq.gz s1_se ${baseR2}.qc.fq.gz s2_se ILLUMINACLIP:/home/jovyan/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10

	gzip -9c s1_se s2_se >> orphans.qc.fq.gz
	rm -f s1_se s2_se
	
done
