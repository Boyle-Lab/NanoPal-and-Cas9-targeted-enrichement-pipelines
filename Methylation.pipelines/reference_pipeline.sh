
#!/bin/bash


EXPERIMENTS=20191028_1854_MN31260_aba684_4994a63f 20191030_1921_MN31260_abb590_545cb6b8 20191106_2111_MN31260_ABB607_dff66530 20191128_2136_MN31260_ABG188_e3688847 20191226_2304_MN31260_ABN780_7e5e331f 20191228_1752_MN31260_ABO515_4ee1d8ca 20200223_2104_MN31260_ACK736_23cfc44e 20200303_2224_MN31260_FAL11389_d5790181 20200304_1521_MN31260_FAL11389_c6908db4 20200304_1532_MN31260_FAL11389_1062c44c 20191009_1753_MN31260_abb908_61bba0f8
for rep in $EXPERIMENTS

do
#Get all fast5 and trimmed fastq data
fastq=/nfs/boylelab_turbo/nanopore_data/MEI/LINE/$rep/output_methyl/pass/porechop_trimmed.fastq.gz
fast5=/nfs/boylelab_turbo/nanopore_data/MEI/LINE/$rep/output_methyl/workspace
#Align to hg38
minimap2 -ax map-ont /data/genomes/hg38/minimap_index/hg38.mmi $fastq | samtools sort -o /home/crmumm/nanopore_methylation/$rep.bam
samtools index /home/crmumm/nanopore_methylation/$rep.bam
#Index all reads to fast5s
../../software/nanopolish/nanopolish index -d $fast5 $fastq
#Call methylation
/home/crmumm/software/nanopolish/nanopolish call-methylation \
		-t 16 \
		-r reads \
		-b /home/crmumm/nanopore_methylation/$rep.bam \
		-g /data/genomes/hg38/minimap_index/hg38.mmi > /home/crmumm/nanopore_methylation/$rep.methylation_calls.tsv
./methylartist db-nanopolish -m $rep.methylation_calls.tsv -d LINE.nanopolish.db -a
done
#Merge bams using list of $rep.bam
samtools merge -b LINE_bams_to_merge.txt LINE_merged.bam
samtools sort -o LINE_merged.sorted.bam

#https://www.ncbi.nlm.nih.gov/nuccore/L19088.1
#L1 bed highlights - https://github.com/adamewing/methylartist
#Default -l = 0.95
#Default max elements = 300
./methylartist composite \
	--maxelts 300 \
	-b LINE_merged.sorted.bam \
	-m LINE.nanopolish.db \
	--sample LINE_merged.sorted_CpG \
	-s L1HS.LINE_merged_data.segmeth.tsv \
	-f L1HS \
	-r /data/genomes/hg38/seq/hg38.fa \
	-t L1.fa \
	-p 32 \
	-c '#2c5280' \
	--blocks L1.highlights.bed
