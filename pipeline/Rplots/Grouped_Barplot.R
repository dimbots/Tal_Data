# Grouped BARPLOTS

# Create data
colors <- c("goldenrod1","dodgerblue3")
names <- c("Stem I", "Stem II", "TA", "Progenitor I", "S phace Cells", "Progenitor II", "Progenitor III", "Goblet/Paneth")
legends <- c("Lgr5Cre","Setdb1KO")

data <- data.frame(values = c(20.3, 21.6, 19.7, 11.9, 15, 14.3, 15, 12.2, 10, 16.4, 9.3, 10, 8.5, 11, 2, 2),  # Create example data
                   group = rep(c("Stem I",
                                 "Stem II",
                                 "TA",
                                 "Progenitor I",
                                 "S phase Cells",
                                 "Progenitor II",
                                 "Progenitor III",
                                 "Goblet/Paneth"),
                               each = 2),
                   subgroup = LETTERS[1:2])



data_base <- reshape(data,                        # Modify data for Base R barplot
                     idvar = "subgroup",
                     timevar = "group",
                     direction = "wide")
row.names(data_base) <- data_base$subgroup
data_base <- data_base[ , 2:ncol(data_base)]
colnames(data_base) <- c("group 1", "group 2", "group 3", "group 4")
data_base <- as.matrix(data_base)
data_base   



pdf("Cells_per_Cluster.pdf", width=11, height=8)
barplot(height = data_base, col = colors, names.arg = names, las =1,   ylim=c(0,25),ylab = "% Cells Normalized", xlab = "Cluster ID",  beside = TRUE)
legend("topright", legends, cex = 1.1, fill = colors)
dev.off()
