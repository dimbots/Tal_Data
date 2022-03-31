# read peaks
setwd("/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/peak_calling/")
ko = read.csv("width_Setdb1KO.bed")
wt = read.csv("width_Lgr5Cre.bed")

ko = ko$X999
wt = wt$X395

# make list
c = list(ko, wt)

names(c) <- c(paste("Setdb1KO\n peaks=" , length(ko) , sep=""), paste("Lgr5Cre\npeaks=" , length(wt) , sep=""))
# Change the mgp argument: avoid text overlaps axis
par(mgp=c(3,2,0))

boxplot(c  , ylab="Assecibility Width (Base Pairs)", outline = FALSE, col = "#69b3a2")


