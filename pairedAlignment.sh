#!/bin/bash

set -e
set -u
set -o pipefail

reference="/home/softgen/NGSTools/SeqMule/database/human_g1k_v37.fasta"
dbsnp_gz="/usr/local/share/bcbio/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz"
normal="/home/softgen/vatsFiles/NA12892_S4_files/NA12892_S4_sortedRG.bam"
normal_pileup="/home/softgen/vatsFiles/NA12892_S4_files/NA12892_latest.pileup"
dbsnp="/home/softgen/NGSTools/SeqMule/database/dbsnp_hg19_138.vcf"
bedfile="/home/vats/tumor_exon.bed"
picard="/home/softgen/vatsFiles/picard/dist/picard.jar"
cleanR1="cleanR1.fastq"
cleanR2="cleanR2.fastq"
result_file="$3_bwa.sam"
result_file2="$3_bwa.bam"
result_file3="$3_bwa_sorted.bam"
result_file4="$3_bwa_sortedRG.bam"
pileup_file="$3_bwa_sortedRG.pileup"
recabTable="$3_recabtable"
recabilatedBam="$3_recabilated.bam"
GATK="/home/softgen/vatsFiles/GenomeAnalysisTK.jar"
varscan="/home/softgen/vatsFiles/varscan/VarScan.v2.4.0.jar"
mutect="/home/softgen/vatsFiles/muTect-1.1.5.jar"
sniper="/home/softgen/vatsFiles/somatic-sniper/build/bin/bam-somaticsniper"
lofreq="/home/softgen/vatsFiles/lofreq_star-2.1.2/bin/lofreq"
strelka="/home/softgen/vatsFiles/strelka_workflow-1.0.14/src/perl/bin/configureStrelkaWorkflow.pl"
strelkaconfig="../config.ini"
sniper_vcf="$3_sniper.vcf"
varscan_vcf="$3_varscan"
lofreq_vcf="$3_lofreq"
strelka_dir="$3_strelkaAnalysis"
bfc="/home/softgen/vatsFiles/bfc/./bfc"
#bfc -t 6 $1 > $cleanR1 
#bfc -t 6 $2 > $cleanR2
#cat $cleanR2 >> $cleanR1
#bwa mem -t 6 $reference $cleanR1 > $result_file 

#samtools view -bS $result_file >$result_file2 	 

#samtools sort -f $result_file2 $result_file3		 

#samtools index $result_file3 

#java -jar $picard AddOrReplaceReadGroups INPUT=$result_file3 OUTPUT=$result_file4 RGLB=Library1 RGPL=Illumina RGPU=runbarcode RGSM=$3
#
#samtools index $result_file4
#
#java -jar $GATK -T BaseRecalibrator -nct 6 -R $reference -I $result_file4 -knownSites $dbsnp -o $recabTable
#
#java -jar $GATK -T PrintReads -nct 6 -R $reference -I $result_file4 -BQSR $recabTable -o $recabilatedBam
#
#rm *.bai $result_file $result_file2 $result_file3 $result_file4 $cleanR2 
#
#samtools index $recabilatedBam
#
#samtools mpileup -q 1 -f $reference $recabilatedBam > $pileup_file 
#
#bam-somaticsniper -q 1 -Q 20 -F vcf -f $reference $recabilatedBam $normal $sniper_vcf
#
#java -jar $varscan somatic $normal_pileup $pileup_file $varscan_vcf --min-var-freq 0.05 --output-vcf 1 
#
#mkdir lofreq
#cd lofreq
#
#lofreq somatic -n $normal -t ../$recabilatedBam -f $reference -o $lofreq_vcf --call-indels --threads 8 
#
#perl $strelka --normal=$normal --tumor=./$recabilatedBam --ref=$reference --config=$strelkaconfig --output-dir=$strelka_dir
#
#cd $strelka_dir
#make -j 6
#
java -jar $mutect -T MuTect --reference_sequence $reference --input_file:$normal --input_file:tumor $recabilatedBam -vcf MuTect.vcf -nt 4
