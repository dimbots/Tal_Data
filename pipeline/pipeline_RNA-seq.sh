#!/bin/bash

# Start in directory fastq
# Include info.tsv

#############################################################################################################################################

	tput setaf 6; tput bold; echo "			"
	tput setaf 1; tput bold; echo "INITIATE PIPELINE"
	tput setaf 2; tput bold; echo " "

	tput setaf 3; tput bold; echo "-----------------------------------"
	tput setaf 3; tput bold; echo "Start Quality Base Check & Trimming"
	tput setaf 3; tput bold; echo "-----------------------------------"
	tput setaf 2; tput bold; echo " "

#	QUALITY BASE CHECK - TRIMMING

for ((x = $1; x <= $2; x++))

do
	sample_fastq="0${x}.fastq"
	sample_trimmed="0${x}.T.fastq"
	fastq="0${x}.fastq"

	fastqc $sample_fastq

	zcat $fastq | echo $((`wc -l`/4)) >> read_count.tmp

	tput setaf 3; tput bold; echo "Trimming - $sample_fastq"
	tput setaf 2; tput bold; echo " "

	TrimmomaticSE -threads 8 $sample_fastq $sample_trimmed SLIDINGWINDOW:4:18 LEADING:28 TRAILING:28 MINLEN:36 >> log_Trimming 2>&1
	paste info.tsv read_count.tmp > fastq_info.tsv

done

	wait

#	Write number of reads in info_table

paste info.tsv read_count.tmp > fastq_info.table.tsv
rm read_count.tmp
rm info.tsv

mkdir trimmed
mv *.T.fastq trimmed/
mv log_Trimming trimmed/

mkdir base_quality
mv *.html *.zip base_quality/

#	SET DIR TRIMMED
	cd trimmed/

fastqc *.gz

#	SET DIR FASTQ
	cd ../

#############################################################################################################################################

#	MAPPING

	tput setaf 3; tput bold; echo "Start Mapping"
	tput setaf 6; tput bold; echo "Set path to reference genome. e.g (/media/dimbo/10T/data/talianidis_data/Genomes/mm10/STAR_index/)"
	tput setaf 2; tput bold; echo " "
	read genome

	mkdir mapping
	cd mapping/
	ln -s ../trimmed/* .

	tput setaf 6; tput bold; echo "Select alignment software. Type -> [HISAT2 OR STAR]"
	tput setaf 2; tput bold; echo " "
	read tool

	if [[ $tool = "HISAT2" ]]
		then

	tput setaf 2; tput bold; echo "Initiate Alignment With HISAT2"

for ((x = $1; x <= $2; x++))

do
	sample="0${x}.testing"
	summary="summary_0${x}.txt"
	fastq_input="0${x}.T.fastq"
	input_sam="0${x}.sam"
	output_bam_tmp="0${x}.bam.tmp"
	output_bam="0${x}.bam"
	output_sorted_bam="0${x}.sorted.bam"

	tput setaf 3; tput bold; echo "processing sample $fastq_input"
	tput setaf 2; tput bold; echo " "

# mapping Hisat2
hisat2 --threads 8 --summary-file $summary -x $genome -U $fastq_input -S $input_sam

# filter reads
#samtools view -@ 8 -h -F 4 $input_sam | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > $output_bam_tmp
samtools view -@ 8 -bq 30 -F 4 $input_sam > $output_bam_tmp

# sort bam
samtools sort -@ 8 $output_bam_tmp -o $output_sorted_bam

# index bam
samtools index $output_sorted_bam

# remove tmp files
rm $input_sam $output_bam.tmp $output_bam

done
	else

	tput setaf 2; tput bold; echo "Initiate Alignment With STAR"

for ((x = $1; x <= $2; x++))

do
        sample="0${x}.testing"
        fastq_input="0${x}.T.fastq"
	file_prefix="0${x}_prefix"
	prefiltered="0${x}_prefixAligned.sortedByCoord.out.bam"
	filtered="0${x}.bam"

	tput setaf 3; tput bold; echo "processing sample $fastq_input"
	tput setaf 2; tput bold; echo " "

# mapping with STAR
STAR --genomeDir $genome --runThreadN 8 --readFilesIn $fastq_input --outFileNamePrefix $file_prefix --outSAMtype BAM SortedByCoordinate

# filter bam
samtools view -@ 8 -bq 30 -F 4 $prefiltered > $filtered
rm $prefiltered

# index
samtools index $filtered

done
	fi

	wait

	mkdir STAR_OUT
	mv *.out STAR_OUT
	mv *.tab STAR_OUT

	# SET DIR TO FASTQ
	cd ../

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "Alignment To Reference Genome Complete!"
	tput setaf 2; tput bold; echo " "

#############################################################################################################################################

# DOWNSAMPLING BASED ON SAMPLE WITH THE LESS NUMBER OF READS

#############################################################################################################################################

	tput setaf 3; tput bold; echo "-------------------"
	tput setaf 3; tput bold; echo "Start Downsampling!"
	tput setaf 3; tput bold; echo "-------------------"
	tput setaf 2; tput bold; echo " "

	# CREATE DIR DOWNSAMPLING
	mkdir downsampling

	# SET DIR TO DOWNSAMPLING
	cd downsampling/

	ln -s ../mapping/*.bam .

#	cat ../mapping/*.txt > mapped_reads.tmp
#	number=$(awk 'NR%3==1' mapped_reads.tmp | grep -v "^ " | cut -d ' ' -f 1 | sort -n | head -1)
#	tput setaf 2; tput bold; echo "Sample with less number of reads is: $number"

	tput setaf 6; tput bold; echo "Type number of reads, e.g 39707152. (This number will be used to downsample all bam files!)"
	read number
	tput setaf 2; tput bold; echo " "

for ((x = $1; x <= $2; x++))

do
	input="0${x}.sorted.bam"
	header="0${x}.tmp"

	samtools view -H $input > $header

	shuffled="0${x}.sam.tmp"

	# -n (this is the sample that has the least mapped reads. (extracted from bam file))
	samtools view -@ 8 $input | shuf | head -n $number > $shuffled

	unsorted="0${x}.downsampled.tmp"

	cat $header $shuffled > $unsorted

	sorted="0${x}.downsampled.bam"

	samtools sort -@ 8 $unsorted -o $sorted

	samtools index -@ 8 $sorted

	bw="0${x}_downsampled.bw"
	bamCoverage -p 8 -b $sorted -o $bw

	rm $shuffled $unsorted
done
	rm *.bam *.bai *.tmp

	tput setaf 3; tput bold; echo "----------------------"
	tput setaf 3; tput bold; echo "Downsampling Complete!"
	tput setaf 3; tput bold; echo "----------------------"
	tput setaf 2; tput bold; echo " "

#	SET DIR FASTQ
	cd ../

