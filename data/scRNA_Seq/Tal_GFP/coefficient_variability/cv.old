#!/bin/bash

# move manually all rows one to the right

# remove first column
cut -f2- KO_Gene_Count_per_Cell.1.tsv > KO_Gene_Count_per_Cell.2.tsv
cut -f2- WT_Gene_Count_per_Cell.1.tsv > WT_Gene_Count_per_Cell.2.tsv
rm KO_Gene_Count_per_Cell.1.tsv WT_Gene_Count_per_Cell.1.tsv

# sum all number of reads
awk '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' KO_Gene_Count_per_Cell.2.tsv > KO_label.tsv
awk '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' WT_Gene_Count_per_Cell.2.tsv > WT_label.tsv
paste KO_Gene_Count_per_Cell.2.tsv KO_label.tsv > KO_Gene_Count_per_Cell.3.tsv
paste WT_Gene_Count_per_Cell.2.tsv WT_label.tsv > WT_Gene_Count_per_Cell.3.tsv
rm KO_Gene_Count_per_Cell.2.tsv WT_Gene_Count_per_Cell.2.tsv KO_label.tsv WT_label.tsv

# read first line and save it
awk 'NR==1 {print; exit}' KO_Gene_Count_per_Cell.3.tsv > ko_label.tmp
awk 'NR==1 {print; exit}' WT_Gene_Count_per_Cell.3.tsv > wt_label.tmp

# exclude genes where the sum of reads of all cells is greater than 50
awk '$NF >10 {print $0}' KO_Gene_Count_per_Cell.3.tsv > KO_Gene_Count_per_Cell.4.tsv   # 8556 genes remain out of 17094
awk '$NF >10 {print $0}' WT_Gene_Count_per_Cell.3.tsv > WT_Gene_Count_per_Cell.4.tsv   # 9351 genes remain out of 17094
rm KO_Gene_Count_per_Cell.3.tsv WT_Gene_Count_per_Cell.3.tsv

# paste label to files
cat wt_label.tmp WT_Gene_Count_per_Cell.4.tsv > WT_matrix_Gene_Count.tsv
cat ko_label.tmp KO_Gene_Count_per_Cell.4.tsv > KO_matrix_Gene_Count.tsv
rm wt_label.tmp ko_label.tmp WT_Gene_Count_per_Cell.4.tsv KO_Gene_Count_per_Cell.4.tsv
# remove intermediate
#rm *.tmp

# remove last column
grep -Po '.*(?=\s+[^\s]+$)' KO_matrix_Gene_Count.tsv > KO_matrix.tsv
grep -Po '.*(?=\s+[^\s]+$)' WT_matrix_Gene_Count.tsv > WT_matrix.tsv
rm KO_matrix_Gene_Count.tsv WT_matrix_Gene_Count.tsv

