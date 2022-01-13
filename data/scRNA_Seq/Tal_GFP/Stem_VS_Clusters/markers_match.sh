#!/bin/bash

# $1: top markers from clusters
# $2: average expression from Setdb1KO
# $3: average expression from WT

while IFS= read -r line

	do

	gene=$(echo $line | awk '{print $1}')

	cat $2 | grep -w $gene >> KO.tmp
	cat $3 | grep -w $gene >> WT.tmp

done < $1

	paste KO.tmp WT.tmp | awk '{print $1 "\t" $2 "\t" $4}' > Markers_KO_WT.tsv

	rm *.tmp
