library(ggplot2)
library(dplyr)
library(hrbrthemes)

data1 = read.csv("/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/Peaks_Compare/final_peaks.tsv", sep = "\t")

# Represent it
p1 <- data1 %>%
  ggplot( aes(x=X452, fill=Setdb1KO)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("dark blue", "dark orange")) +
  theme_classic() +
  labs(fill="") +
  xlim(c(400,2000)) +
  xlab("Accesibility Width (base pairs)") +
  ylab("Peaks") +
pdf("Peak_width.pdf", width = 12,
    height = 8)
p1
dev.off()

