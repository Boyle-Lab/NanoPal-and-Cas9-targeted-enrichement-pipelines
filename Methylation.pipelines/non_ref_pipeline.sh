
#!/bin/bash
LIST=$(ls non_ref_contig.fasta.split/ | grep fasta$)
for item in $LIST
do
chr_item=$(echo $item | cut -d'.' -f2 | cut -d'_' -f2)
position=$(echo $item | cut -d'.' -f3)
bottom_bound=$(($position - 250))
top_bound=$(($position + 250))
coord=$chr_item:$bottom_bound-$top_bound
##
#Get reads from pooled that support element (500bp window)
samtools view -h -b -q 10 filtered.merged.sorted.bam $coord > bams/$coord.bam
samtools fastq bams/$coord.bam > reads/$coord.fastq
##
reference_contig=non_ref_contig.fasta.split/$item
reads=reads/$coord.fastq
#Map extracted reads back to non-ref contig
minimap2 -ax map-ont $reference_contig $reads | samtools sort -o /home/crmumm/nanopore_methylation/non_ref/realign/$coord.sorted.bam
samtools index /home/crmumm/nanopore_methylation/non_ref/realign/$coord.sorted.bam
#Index the extracted reads back to available fast5 directories
/home/crmumm/software/nanopolish/nanopolish index \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20191028_1854_MN31260_aba684_4994a63f/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20191030_1921_MN31260_abb590_545cb6b8/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20191106_2111_MN31260_ABB607_dff66530/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20190909_1740_MN31260_AAL970_09e46393/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20191003_1754_MN31260_aba133_adc46a9d/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20191128_2136_MN31260_ABG188_e3688847/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20200223_2104_MN31260_ACK736_23cfc44e/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20200303_2224_MN31260_FAL11389_d5790181/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20200304_1521_MN31260_FAL11389_c6908db4/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20200304_1532_MN31260_FAL11389_1062c44c/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/LINE/20191009_1753_MN31260_abb908_61bba0f8/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/Bulk/20210103_2326_MN31260_FAO08905_c9e9b039/output_methyl/workspace \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/Bulk/20201021_1942_MN34990_FAO84736_bcfa646e/output_methyl/workspace/fast5 \
	-d /nfs/boylelab_turbo/nanopore_data/MEI/Bulk/20201022_1801_MN34990_FAO84736_31aef5f8/output_methyl/workspace/fast5 \
	$reads
#echo 'Done indexing'
/home/crmumm/software/nanopolish/nanopolish call-methylation \
		-t 8 \
		-r $reads \
		-b realign/$coord.sorted.bam \
		-g $reference_contig > /home/crmumm/nanopore_methylation/non_ref/nanopolish/$coord.methylation_calls.tsv

#Aggregate nanopolish data with methylartist
./methylartist db-nanopolish \
	-m /home/crmumm/nanopore_methylation/non_ref/nanopolish/$coord.methylation_calls.tsv \
	-d /home/crmumm/nanopore_methylation/non_ref/methyl_db/non_ref.nanopolish.db -a
done

#FILES=$(ls *.sorted.bam)
for f in $FILES
do
name=$(echo $f | cut -d'.' -f1)
echo realign/$f	/home/crmumm/nanopore_methylation/non_ref/methyl_db/non_ref.nanopolish.db > $name.txt
done

###Second half of pipeline - Plotting

INTERVAL=$(ls non_ref_contig.fasta.split/ | grep fasta$)

for item in $INTERVAL:
do
chr_item=$(echo $item | cut -d'.' -f2 | cut -d'_' -f2)
position=$(echo $item | cut -d'.' -f3)
locus_bound_line=$(($position + 6000))
name_low=$(($position - 250))
name_high=$(($position + 250))
name_coord=$chr_item'_'$name_low-$name_high
strand=$(echo $item | cut -d'.' -f4)
coord=$chr_item.$position.$strand:50000-156000
echo -e $coord'\t100000 \t106000 \tL1HS \t'$strand >> beds/nonref.bed
../methylartist locus \
	-d realign/$name_coord.txt \
	-i $coord
done
#Merge bams for composite
samtools merge merge.sorted.bam realign/*.sorted.bam
#Make segmeth tsv for composite
./methylartist segmeth \
	-d nonref_merged.txt \
	-i beds/mod.nonref.bed

#Needed to manually change the strand info in segmeth output!
#Lower the minimum profiles because there is only 54 L1HS elements
./methylartist composite \
	--outelts 50 \
	-l 0.95 \
	-b realign/merge.sorted.bam \
	--sample merge.sorted_CpG \
	-m methyl_db/non_ref.nanopolish.db \
	-s beds/mod.nonref.nonref_merged.segmeth.tsv \
	-f L1HS \
	-c '#b73f75' \
	-r nonrefref.fasta \
	-t ../L1.fa \
	--blocks ../L1.highlights.bed
