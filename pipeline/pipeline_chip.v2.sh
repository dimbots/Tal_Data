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
	sample_fastq="0${x}.fastq.gz"
	sample_trimmed="0${x}.T.fastq.gz"
	fastq="0${x}.fastq.gz"

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
mv *.T.fastq.gz trimmed/
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
	tput setaf 6; tput bold; echo "Set path to reference genome. e.g (/media/dimbo/10T/data/talianidis_data/genomes/mm10/hisat_index/mm10) "
	tput setaf 2; tput bold; echo " "
	read genome

	mkdir mapping
	ln -s trimmed/* .

	tput setaf 6; tput bold; echo "Select alignment software. Type -> [BOWTIE2 OR HISAT2]"
	tput setaf 2; tput bold; echo " "
	read tool

	if [[ $tool = "HISAT2" ]]
		then

	tput setaf 2; tput bold; echo "Initiate Alignment With HISAT2"

for ((x = $1; x <= $2; x++))

do
	sample="0${x}.testing"
	summary="summary_0${x}.txt"
	fastq_input="0${x}.T.fastq.gz"
	input_sam="0${x}.sam"
	output_bam_tmp="0${x}.bam.tmp"
	output_bam="0${x}.bam"
	output_sorted_bam="0${x}.sorted.bam"

	tput setaf 3; tput bold; echo "processing sample $fastq_input"
	tput setaf 2; tput bold; echo " "

# mapping Hisat2
hisat2 --threads 8 --summary-file $summary -x $genome -U $fastq_input -S $input_sam

# keep only unique
samtools view -@ 8 -h -F 4 $input_sam | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > $output_bam_tmp

# remove dublicates
samtools rmdup -s $output_bam_tmp $output_bam

# sort bam
samtools sort -@ 8 $output_bam -o $output_sorted_bam

# index bam
samtools index $output_sorted_bam

# remove tmp files
rm $input_sam $output_bam.tmp $output_bam

done
	else

	tput setaf 2; tput bold; echo "Initiate Alignment With BOWTIE2"

for ((x = $1; x <= $2; x++))

do
        sample="0${x}.testing"
        summary="summary_0${x}.txt"
        fastq_input="0${x}.T.fastq.gz"
        input_sam="0${x}.sam"
        unsorted_bam="0${x}.unsorted.bam.tmp"
        sorted_bam="0${x}.bam"

	tput setaf 3; tput bold; echo "processing sample $fastq_input"
	tput setaf 2; tput bold; echo " "

# mapping bowite2
bowtie2 -p 8 --sensitive-local -x $genome -U $fastq_input -S $input_sam >> $summary 2>&1

# filter reads
samtools view -h -S -b -q 30 -o $unsorted_bam $input_sam

# sort bam
samtools sort -t 8 -o $sorted_bam $unsorted_bam

# index
samtools index $sorted_bam

# remove tmp files
rm $input_sam
rm *.tmp

done

	fi

	wait

	mv *.bam *.bai *.txt mapping/

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "Alignment To Reference Genome Complete!"
	tput setaf 2; tput bold; echo " "

#############################################################################################################################################

#	CREATE BAR PLOTS (R) FOR UNIQUE-REPEAT ALIGNMENTS

#	SET DIR MAPPING
	cd mapping/

for ((x = $1; x <= $2; x++))

do
	sum_file="summary_0${x}.txt"

	uniq="uniq.tmp"
	repeat="repeat.tmp"
	uniq_repeat="uniq_repeat.tsv"
	identifier="identifier.tsv"
	id="0${x}.sorted.bam"

	awk 'NR % 4 == 0' $sum_file | awk '{print $2}' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' >> $uniq
	awk 'NR % 5 == 0' $sum_file | awk '{print $2}' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' >> $repeat
	echo $id | cut -f1 -d "." >> $identifier

done

	wait

	cat $uniq $repeat > $uniq_repeat
	ln -s ../fastq_info.table.tsv .
	rm *tmp

	num_samples=$(wc -l identifier.tsv | awk '{print $1}')

	awk '{print $14,$15}' fastq_info.table.tsv > identifier_new.tmp
	head -$num_samples identifier_new.tmp | tr -d "[:blank:]" > final_IDS.tsv

	working_directory=$(pwd)
	ncol=$(cat $uniq_repeat | echo $((`wc -l`/2)))

	paste identifier.tsv final_IDS.tsv > sample_IDs.tsv

#	Write R script for barplots and save it to file

	echo	"setwd(\"$working_directory\")
		 colors = c(\"paleturquoise3\",\"red4\")
		 names = scan(\"identifier.tsv\", character(), quote = \"\")
		 values = scan(\"uniq_repeat.tsv\")
		 values = matrix(values, nrow = 2, ncol = $ncol, byrow = TRUE)
		 pdf(file = \"$working_directory/plot.pdf\", width = 7, height = 7)
		 barplot(values, names.arg = names, ylab = \"Alignment % Rate \", col = colors, density = 300, ylim = c(0,100) ,cex.lab = 1.4, las=3, cex.names = 1)
		 dev.off()" > R_Plot.R

		chmod 755 R_Plot.R
		Rscript R_Plot.R

	rm final_IDS.tsv identifier* uniq_repeat.tsv

#	SET DIR FASTQ
	cd ../

# Merge bam files (replicate 1 & 2) Treatment and control
# samtools merge Set8KO_TCP_input.merged.bam Set8KO_TCP_A_input.bam Set8KO_TCP_B_input.bam

	tput setaf 3; tput bold; echo "---------------------------------------------------------------------------------------------------------------"
	tput setaf 3; tput bold; echo "                                 Bar Plots For Alignment Rates Created!                                        "
	tput setaf 3; tput bold; echo "---------------------------------------------------------------------------------------------------------------"
	tput setaf 2; tput bold; echo " "

#############################################################################################################################################

# PEAK CALLING

#############################################################################################################################################

	tput setaf 3; tput bold; echo "---------------------------------------------------------------------------------------------------------------"
	tput setaf 3; tput bold; echo "                               Initiate Peak Calling Analysis with MACS2!                                      "
	tput setaf 3; tput bold; echo "---------------------------------------------------------------------------------------------------------------"
	tput setaf 2; tput bold; echo " "

	mkdir peak_calling

#	SET DIR PEAK_CALLING
	cd peak_calling/

	ln -s ../mapping/*.bam .
	ln -s ../mapping/*.bam.bai .

		while [[ $treatment != "none" ]]

		do

		tput setaf 6; tput bold; echo "Type treatment bam file. else type [none]"
		read treatment

		tput setaf 6; tput bold; echo "Type input bam file. else type [none]"
		read input

		tput setaf 6; tput bold; echo "Type peak call out file. Format: [condition_rep_A]. else type [none]"
		tput setaf 2; tput bold; echo " "
		read out

			if	[[ $treatment = "none" ]]
				then
				break
			else

#			macs2 callpeak --treatment $treatment --control $input --nomodel --broad --broad-cutoff=0.01 --format BAM --gsize mm -n $out
			macs2 callpeak --treatment $treatment --control $input --nomodel --broad --format BAM --gsize mm -n $out
		fi

		done


#		MERGE PEAKS FROM REPLICATES USING BEDTOOLS INTERSECT

		mkdir merged_peaks

		while [[ $repA != "none" ]]

		do

		tput setaf 3; tput bold; echo "----------------------------------"
		tput setaf 3; tput bold; echo "Merge Peaks From Replicate A and B"
		tput setaf 3; tput bold; echo "----------------------------------"

		tput setaf 6; tput bold; echo "type peaks out file from replicate A. else type [none]"
		read repA

		tput setaf 6; tput bold; echo "type peaks out file from replicate B. else type [none]"
		read repB

		tput setaf 6; tput bold; echo "type overlapped peaks out file. Format: [overlapped_condition]. else type [none]"
		tput setaf 2; tput bold; echo " "
		read out_file

			if      [[ $repA = "none" ]]
				then
				break
		else

		bedtools intersect -a $repA -b $repB -wa > $out_file



	# Exculde black listed regions

	final_peaks="$out_file.bed"

	tput setaf 6; tput bold; echo "set path to blacklist_regions. else type [none]"
	tput setaf 2; tput bold; echo " "
	read blacklist

	bedtools intersect -v -a $out_file -b $blacklist > $final_peaks

		fi

		done

	mv *.bed merged_peaks/
	rm *.tmp


####################################################################################################################

# QC Metrics (Deeptools) + Normalization

####################################################################################################################

		tput setaf 3; tput bold; echo "--------------------------------------------------------------------"
		tput setaf 3; tput bold; echo "                   Initiate QC Metrics Analysis                     "
		tput setaf 3; tput bold; echo "--------------------------------------------------------------------"
		tput setaf 2; tput bold; echo " "

#	START IN FASTQ DIRECTORY
	cd ..
#	SET DIR QC_METRICS
	mkdir qc_metrics
	cd qc_metrics/

	ln -s ../mapping/*.bam .
	ln -s ../mapping/*.bam.bai .

	tput setaf 6; tput bold; echo "In directory (qc_metrics) rename files from (0869.sorted.bam) to (Set8KOA_Input)"
	tput setaf 6; tput bold; echo "Type [done] when renaming of files is complete."
	tput setaf 2; tput bold; echo " "
	read response

# CREATE MULTIBAM FILE

	renamed_files=$(ls --ignore=*.bai)

multiBamSummary bins --bamfiles $renamed_files -p 8 -o multiBam.npz

# PLOT SPEARMAN AND PEARSON CORRELATION
plotCorrelation -in multiBam.npz --corMethod spearman --removeOutliers --skipZeros --colorMap Blues --plotHeight 11.5 --plotWidth 13  --whatToPlot heatmap --plotNumbers -o SpearmanCor_readCounts_plot.png
plotCorrelation -in multiBam.npz --corMethod pearson --removeOutliers  --skipZeros --colorMap Blues --plotHeight 11.5 --plotWidth 13  --whatToPlot heatmap --plotNumbers -o PearsonCor_readCounts_plot.png

# PLOT READ COVERAGE
plotCoverage --bamfiles $renamed_files --skipZeros -p 8 --verbose -o Coverage_plot.png

# PLOT FINGERPRINT
plotFingerprint -b $renamed_files -p 8 -plot FingerPrint_plot.png

# PLOT PCA
plotPCA -in multiBam.npz -o PCA_readCounts.png -T "PCA of read counts"

	tput setaf 3; tput bold; echo "QC Metrics Complete!"
	tput setaf 2; tput bold; echo " "

#########################################################################################

#	NORMALIZE BAM FILES

#	CONVERT BAM TO BIGWIG (NORMALIZATION)

	tput setaf 3; tput bold; echo "--------------------------------------------------------------------------"
	tput setaf 3; tput bold; echo "			       Normalize Bam Files!				 "
	tput setaf 3; tput bold; echo "--------------------------------------------------------------------------"
	tput setaf 2; tput bold; echo " "

	tput setaf 3; tput bold; echo "Type Method of Normalization. [RPKM or BPM or RPGC]"
	tput setaf 2; tput bold; echo " "
	read method

	if
	[[ $method = "BPM" ]]
			then

			for i in $(ls -I "*.bai" -I "*.npz" -I "*.png")

			do
				bw="$i.bw"
				bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing BPM -p 8 --extendReads 200
			done
	fi



		if

		[[ $method = "RPKM" ]]
		then
			for i in $(ls -I "*.bai" -I "*.npz" -I "*.png")

			do

			bw="$i.bw"
			bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPKM -p 8 --extendReads 200

			done
		fi



	if
	[[ $method = "RPGC" ]]
	then

		for i in $(ls -I "*.bai" -I "*.npz" -I "*.png")

		do
		bw="$i.bw"
		bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPGC -p 8 --effectiveGenomeSize 2407883318 --extendReads 200

		done
	fi

	mkdir normalization
	mv *.bw normalization/

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "Normalization Comlete!"
	tput setaf 2; tput bold; echo " "


#########################################################################################

# 			COMPUTEMATRIX

#	SET DIR NORMALIZATION
	cd normalization/

	tput setaf 3; tput bold; echo "          Start ComputeMatrix            "
	tput setaf 3; tput bold; echo " "
	tput setaf 3; tput bold; echo "            Create Heatmaps              "
	tput setaf 3; tput bold; echo " "

	while [[ $num_samples != "none" ]]

		do

	tput setaf 6; tput bold; echo "Type how many samples to run [2 or 4]. When finished type [none]"
	read num_samples
	tput setaf 2; tput bold; echo " "

	if      [[ $num_samples = "none" ]]
		then
		break
	else
		if
		[[ $num_samples = "4" ]]
		then

		echo "Run compute matrix with 4 samples!"

	# RUN WITH 4 SAMPLES
	tput setaf 6; tput bold; echo "Type treatment bw file for replicate A"
	read bw1
	tput setaf 6; tput bold; echo "Type input bw file for replicate A"
	read bw2
	tput setaf 6; tput bold; echo "Type treatment bw file for replicate B"
	read bw3
	tput setaf 6; tput bold; echo "Type input bw file for replicate B"
	read bw4
	tput setaf 6; tput bold; echo "Type matrix out file. e.g (matrix_Set8KO)"
	read out_matrix
	tput setaf 6; tput bold; echo "Type regions out file. e.g (out_regions_Set8KO)"
	read out_regions
	tput setaf 6; tput bold; echo "Set path to genes.bed file. e.g (/media/dimbo/10T/data/talianidis_data/Genomes/mm10/genes_info/mm10_GENES_GENCODE.VM23.bed)"
	read genes
	tput setaf 6; tput bold; echo "Type number of base pairs before Transcription Start Site"
	read before_tss
	tput setaf 6; tput bold; echo "Type number of base pairs after Transcription End Site"
	read after_tss
	tput setaf 6; tput bold; echo "Type number of base pairs for Region Body Length"
	read region_body
	tput setaf 2; tput bold; echo " "

	# Scale regions
	computeMatrix scale-regions --startLabel "TSS" --endLabel "TES" -b $before_tss -a $after_tss --regionBodyLength $region_body -R $genes -S $bw1 $bw2 $bw3 $bw4 --smartLabels --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 8

	# Reference point
#	computeMatrix reference-point -b $before_tss -a $after_tss -R $genes -S $bw1 $bw2 $bw3 $bw4 --skipZeros --smartLabels -o $out_matrix --outFileSortedRegions $out_regions -p 8

		out_plot="Heatmap_$out_matrix"
		plotHeatmap -m $out_matrix -out $out_plot --colorMap Blues --boxAroundHeatmaps no --missingDataColor 1

		out_profile_plot="ProfilePlot_$out_matrix"
		plotProfile -m $out_matrix -out $out_profile_plot --colors blue blue

		else

		echo "Run compute matrix with 2 samples!"

	# RUN WITH 2 SAMPLES
	tput setaf 6; tput bold; echo "Type treatment bw file"
        read bw1
        tput setaf 6; tput bold; echo "Type input bw file"
	read bw2
	tput setaf 6; tput bold; echo "Type matrix out file. e.g (matrix_Set8KOA)"
	read out_matrix
	tput setaf 6; tput bold; echo "Type regions out file. e.g (out_regions_Set8KOA)"
	read out_regions
	tput setaf 6; tput bold; echo "Set path to genes.bed file. e.g (/media/dimbo/10T/data/talianidis_data/Genomes/mm10/genes_info/mm10_GENES_GENCODE.VM23.bed)"
	read genes
	tput setaf 6; tput bold; echo "Type number of base pairs before Transcription Start Site"
	read before_tss
	tput setaf 6; tput bold; echo "Type number of base pairs after Transcription End Site"
	read after_tss
	tput setaf 6; tput bold; echo "Type number of base pairs for Region Body Length"
	read region_body
	tput setaf 2; tput bold; echo " "

	# Scale regions
	computeMatrix scale-regions --startLabel "TSS" --endLabel "TES" -b $before_tss -a $after_tss --regionBodyLength $region_body -R $genes -S $bw1 $bw2 --smartLabels --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 8

	# Reference point
#	computeMatrix reference-point -b $before_tss -a $after_tss -R $genes -S $bw1 $bw2 --skipZeros --smartLabels -o $out_matrix --outFileSortedRegions $out_regions -p 8

		out_plot="Heatmap_$out_matrix"
		plotHeatmap -m $out_matrix -out $out_plot --colorMap Blues --boxAroundHeatmaps no --missingDataColor 1

		out_profile_plot="ProfilePlot_$out_matrix"
		plotProfile -m $out_matrix -out $out_profile_plot --colors blue blue

		fi

	fi

	done

mv $out_plot ../



#########################################################################################

#	ANNOTATING REGION IN THE GENOME

#       HOMER

#	SET DIR FASTQ
	cd ../../
#	MAKE DIR ANNOTATIONS
	mkdir annotation
	cd annotation/

	tput setaf 6; tput bold; echo "SET PATH TO REFERENCE GENOME"
	read ref_genome
	tput setaf 6; tput bold; echo "SET PATH TO GENES.GTF FILE"
	read genes_gtf
	tput setaf 2; tput bold; echo " "

	ln -s ../peak_calling/merged_peaks/overlapped_* .

	while [[ $PEAKS != "none" ]]

		do

	tput setaf 6; tput bold; echo "TYPE PEAKS.BED FILES ELSE TYPE [none]"
	read PEAKS
	tput setaf 6; tput bold; echo "TYPE ANNOATED FILE E.G(annotations_Set8Ko.tsv) ELSE TYPE [none]"
	read out_annotation
	tput setaf 2; tput bold; echo " "

		if      [[ $PEAKS = "none" ]]
			then
			break
		else

		annotatePeaks.pl $PEAKS $ref_genome -gtf $genes_gtf > $out_annotation
	fi
	done


	tput setaf 3; tput bold; echo "-------------------"
	tput setaf 3; tput bold; echo "ANNOTATION COMPLETE"
	tput setaf 3; tput bold; echo "-------------------"
	tput setaf 2; tput bold; echo " "

#########################################################################################

# DOWNSAMPLING BASED ON SAMPLE WITH THE LESS NUMBER OF READS

	tput setaf 3; tput bold; echo "-------------------"
	tput setaf 3; tput bold; echo "START DOWNSAMPLING "
	tput setaf 3; tput bold; echo "-------------------"
	tput setaf 2; tput bold; echo " "

# SET DIR TO FASTQ
cd ..
# CREATE DIR DOWNSAMPLING
mkdir downsampling
# SET DIR TO DOWNSAMPLING
cd downsampling/

ln -s ../mapping/*.bam .

	cat ../mapping/*.txt > mapped_reads.tmp
	number=$(awk 'NR%3==1' mapped_reads.tmp | grep -v "^ " | cut -d ' ' -f 1 | sort -n | head -1)

	tput setaf 2; tput bold; echo "MAPPED SAMPLE WITH LESS NUMBER OF READS IS: $number"

for ((x = $1; x <= $2; x++))

do

	input="${x}.sorted.bam"
	header="${x}.tmp"

	samtools view -H $input > $header

	shuffled="${x}.sam.tmp"

	# -n (this is the sample that has the least mapped reads. (extracted from bam file))
	samtools view -@ 7 $input | shuf | head -n $number > $shuffled

	unsorted="${x}.downsampled.tmp"

	cat $header $shuffled > $unsorted

	sorted="${x}.downsampled.bam"

	samtools sort -@ 7 $unsorted -o $sorted

	samtools index -@ 7 $sorted

	bw="${x}_downsampled.bw"
	bamCoverage -p 7 -b $sorted -o $bw 

	rm $shuffled $unsorted

done

	rm *.bam *.bai *.tmp

	tput setaf 6; tput bold; echo "---------------------"
	tput setaf 2; tput bold; echo "DOWNSAMPLING COMPLETE"
	tput setaf 6; tput bold; echo "---------------------"
	tput setaf 2; tput bold; echo " "

#############################################################################################################################################

# IDENTIFY SUPER ENHANCERS (ROSE)

#############################################################################################################################################


#!/bin/bash

tput setaf 2; tput bold; echo "                     "
tput setaf 2; tput bold; echo "                     "
tput setaf 2; tput bold; echo "                     "

tput setaf 6; tput bold; echo "INITIATE SUPER ENHANCER ANALYSIS? TYPE Y or N"
tput setaf 2; tput bold; echo "                     "

	read SE

	if
	[[ $SE = "Y" ]]

	then

#	SET DIR FASTQ
	cd ../
#	CREATE SE DIR
	mkdir SE
#	SET SE DIR
	cd SE/

	ln -s ../peak_calling/merged_peaks/*.bed .
	ln -s ../mapping/*.bam .

	tput setaf 6; tput bold; echo "RENAME BAM FILES WHEN FINISH TYPE OK"
	read resp

# MERGE BAM FILES

	while [[ $merged != "none" ]]

		do

		tput setaf 6; tput bold; echo "TYPE MERGED BAM FILE. E.G Set8Ko_merged.bam.  WHEN FINISHED TYPE none"
		read merged
		tput setaf 6; tput bold; echo "TYPE REPLICATE A BAM FILE. WHEN FINISHED TYPE none"
		read rep_A
		tput setaf 6; tput bold; echo "TYPE REPLICATE B BAM FILE. WHEN FINISHED TYPE none"
		read rep_B
		tput setaf 2; tput bold; echo " "

		if	[[ $merged = "none" ]]
			then
			break
		else

		samtools merge -@ 8 $merged $rep_A $rep_B
		samtools index $merged
		fi
		done

# MAKE GFF FILES

	while [[ $peaks_bed != "none" ]]

		do

		tput setaf 6; tput bold; echo "TYPE PEAKS BED FILE. ELSE TYPE none"
		tput setaf 2; tput bold; echo " "
		read peaks_bed

		if     [[ $peaks_bed = "none" ]]
			then
			break

		else

# Remove unique IDS (in order to run ROSE)
		sort -u -k4 $peaks_bed > "$peaks_bed_uniq.tmp"

# convert to ROSE gff Format
		gff="$peaks_bed.gff"
		awk '{print $1,$4,$10,$2,$3,$10,$6,$10,$4}' OFS='\t' "$peaks_bed_uniq.tmp" > $gff

		fi
		done
# Run ROSE - IDENTIFY SUPER ENHACERS

	mkdir run_rose
	cd run_rose
	cp -r /home/dimbo/rose/* .
	ln -s ../* .

	while [[ $merged_gff != "none" ]]

	do

	tput setaf 6; tput bold; echo "TYPE GFF FILE. ELSE TYPE none"
	read merged_gff

	tput setaf 6; tput bold; echo "TYPE MERGED TREATMENT BAM FILE. ELSE TYPE none"
	read treatment_merged_bam

	tput setaf 6; tput bold; echo "TYPE MERGED INPUT BAM FILE. ELSE TYPE none"
	read input_merged_bam

	tput setaf 6; tput bold; echo "SE OUT FILE. ELSE TYPE none"
	tput setaf 2; tput bold; echo " "
	read out

	if      [[ $merged_gff = "none" ]]
		then
		break

		else

		mkdir $out

# Run ROSE. Identify SUPER ENHANCERS
python ROSE_main.py --genome MM10 --i $merged_gff --rankby $treatment_merged_bam --control $input_merged_bam --out $out
		fi
		done


else
	echo "SUPER ENHANCER ANALYSIS ABORTED"
fi








































#############################################################################################################################################

# IDENTIFY SUPER ENHANCERS (ROSE)

#############################################################################################################################################



#	SET DIR FASTQ
#	cd ../
#	CREATE SE DIR
#	mkdir SE
#	SET SE DIR
#	cd SE/

#	ln -s ../peak_calling/*.bed .

#	Remove unique IDS (in order to run ROSE)
#	sort -u -k4 rep_Filtered_H3K27ac_WTuntr_overlap_peaks.bed > rep_Filtered_Unique_H3K27ac_WTuntr_overlap_peaks.bed (H3K27ac_Set8KO_FU.mergedPeaks.bed)

# convert to ROSE gff Format
#awk '{print $1,$4,$10,$2,$3,$10,$6,$10,$4}' OFS='\t' rep_H3K27ac_Set8KO.overlaps.broadPeak > rep_H3K27ac_Set8KO.overlaps.ROSE.gff (H3K27ac_Set8KO_FU.mergedPeaks.gff)

# Run ROSE. Identify SUPER ENHANCERS
#python ROSE_main.py --genome MM10 --i H3K27ac_Set8KO_FU.mergedPeaks.gff --rankby Set8KO_TCP_H3K27ac.merged.bam --control Set8KO_TCP_input.merged.bam --out SE_Set8KO/
#python ROSE_main.py --genome MM10 --i H3K27ac_WTuntr_FU.mergedPeaks.gff --rankby WTuntr_H3K27ac.merged.bam --control WTuntr_input.merged.bam --out SE_WTuntr/


######## Super enhancer plots

#awk '{ total += $6 } END { print total/NR }' SE_sorted.tsv


#computeMatrix scale-regions -b 5000 -a 5000 --regionBodyLength 2780 -R TE.bed -S Set8KO_TCP_A_H3K27ac.bw --skipZeros -o output_matrix --outFileSortedRegions output_regions -p 6



#plotProfile --averageType mean --regionsLabel Set8KO --plotType lines --colors red --startLabel start --endLabel end --samplesLabel "Genome Wide Average" --yAxisLabel "reads per genome coverage, RPGC" --yMax 60 --plotWidth 6 -m output_matrix -out example
