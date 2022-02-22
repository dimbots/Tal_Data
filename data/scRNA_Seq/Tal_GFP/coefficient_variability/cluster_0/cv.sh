#!/bin/bash

awk '{print $NF}' KO.tmp > KO_cv.tmp
awk '{print $NF}' WT.tmp > WT_cv.tmp

paste KO_cv.tmp WT_cv.tmp > cov_values.tsv.tmp

awk '$1 != "NA" && $2 != "NA" {print $0}' cov_values.tsv.tmp > coefficient_variation.tsv


rm *.tmp
