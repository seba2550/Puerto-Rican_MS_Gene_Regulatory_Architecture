create table ATAC_Table
AS
select 338 as Ind_ID, * from wbush_atacseqset2_201816338_01_s_6_1_peaks union
select 339, * from wbush_atacseqset2_201816339_01_s_6_1_peaks union
select 340, * from wbush_atacseqset2_201816340_01_s_6_1_peaks union
select 341, * from wbush_atacseqset2_201816341_01_s_6_1_peaks union
select 342, * from wbush_atacseqset2_201816342_01_s_6_1_peaks union
select 343, * from wbush_atacseqset2_201816343_01_s_6_1_peaks union
select 344, * from wbush_atacseqset2_201816344_01_s_6_1_peaks union
select 345, * from wbush_atacseqset2_201816345_01_s_6_1_peaks union
select 346, * from wbush_atacseqset2_201816346_01_s_6_1_peaks union
select 347, * from wbush_atacseqset2_201816347_01_s_6_1_peaks union
select 348, * from wbush_atacseqset2_201816348_01_s_6_1_peaks union
select 349, * from wbush_atacseqset2_201816349_01_s_6_1_peaks union
select 350, * from wbush_atacseqset2_201816350_01_s_6_1_peaks union
select 352, * from wbush_atacseqset2_201816352_01_s_6_1_peaks; 

select * from ATAC_Table;

select right (chrom_name, 2), * from atac_table; 

select trim(both 'chr' from chrom_name) as chrom_num from atac_table;

update ATAC_Table 
set chrom_name = trim(both 'chr' from chrom_name);



CREATE TABLE RNAseq_Table
AS
select 338 as Ind_ID, symbol, "wBush_RNASeqset2_201816427.01" from rnaseq_counts union
select 339 as Ind_ID, symbol, "wBush_RNASeqset2_201816430.01" from rnaseq_counts union
select 340 as Ind_ID, symbol, "wBush_RNASeqset2_201816433.01" from rnaseq_counts union
select 341 as Ind_ID, symbol, "wBush_RNASeqset2_201816436.01" from rnaseq_counts union
select 342 as Ind_ID, symbol, "wBush_RNASeqset2_201816439.01" from rnaseq_counts union
select 343 as Ind_ID, symbol, "wBush_RNASeqset2_201816428.01" from rnaseq_counts union
select 344 as Ind_ID, symbol, "wBush_RNASeqset2_201816431.01" from rnaseq_counts union
select 345 as Ind_ID, symbol, "wBush_RNASeqset2_201816434.01" from rnaseq_counts union
select 346 as Ind_ID, symbol, "wBush_RNASeqset2_201816437.01" from rnaseq_counts union
select 347 as Ind_ID, symbol, "wBush_RNASeqset2_201816440.01" from rnaseq_counts union
select 348 as Ind_ID, symbol, "wBush_RNASeqset2_201816429.01" from rnaseq_counts union
select 349 as Ind_ID, symbol, "wBush_RNASeqset2_201816432.01" from rnaseq_counts union
select 350 as Ind_ID, symbol, "wBush_RNASeqset2_201816435.01" from rnaseq_counts union
select 352 as Ind_ID, symbol, "wBush_RNASeqset2_201816433.01" from rnaseq_counts ;

select * from RNAseq_Table;
select * from atac_table;
select * from seq_region;
select * from gene;

create index atac_index on atac_table (ind_id, chrom_name, peak_start);
create index rna_index on RNAseq_Table (ind_id, symbol);

explain
select LEFT(symbol, 15) as modified_gene_id, * from RNAseq_Table rna inner join gene g on LEFT(rna.symbol, 15) = g.stable_id
inner join seq_region sr on g.seq_region_id = sr.seq_region_id inner join atac_table atac
on rna.ind_id = atac.ind_id and (atac.peak_start between g.seq_region_start and g.seq_region_end) and sr."name" = atac.chrom_name ;


create table atac_rna_merged_5kb_table
as
select LEFT(symbol, 15) as modified_gene_id, rna.ind_id, gene_counts, g.seq_region_id, g.seq_region_start, g.seq_region_end, g.description, g.stable_id,
chrom_name, peak_start, peak_end, peak_name, peak_score  from RNAseq_Table rna inner join gene g on LEFT(rna.symbol, 15) = g.stable_id
inner join seq_region sr on g.seq_region_id = sr.seq_region_id inner join atac_table atac
on rna.ind_id = atac.ind_id and (atac.peak_start between g.seq_region_start-500000 and g.seq_region_end+500000) and sr."name" = atac.chrom_name ;

create table atac_rna_merged_table
as
select LEFT(symbol, 15) as modified_gene_id, rna.ind_id, gene_counts, g.seq_region_id, g.seq_region_start, g.seq_region_end, g.description, g.stable_id,
chrom_name, peak_start, peak_end, peak_name, peak_score  from RNAseq_Table rna inner join gene g on LEFT(rna.symbol, 15) = g.stable_id
inner join seq_region sr on g.seq_region_id = sr.seq_region_id inner join atac_table atac
on rna.ind_id = atac.ind_id and (atac.peak_start between g.seq_region_start and g.seq_region_end) and sr."name" = atac.chrom_name ;



select * from atac_rna_merged_5kb_table;
select * from atac_rna_merged_table;

select corr(gene_counts::numeric ,peak_score), modified_gene_id, description from atac_rna_merged_5kb_table group by modified_gene_id, description;


select corr(log(gene_counts::numeric + 1) , log(peak_score + 1)), modified_gene_id, description from atac_rna_merged_5kb_table
where gene_counts::numeric > 0 and peak_score > 0 and gene_counts notnull and peak_score notnull 
group by modified_gene_id, description
order by corr(log(gene_counts::numeric + 1) , log(peak_score + 1)) desc;


select regr_r2(gene_counts::numeric, peak_score), modified_gene_id, description from atac_rna_merged_5kb_table armkt group by modified_gene_id, description;






