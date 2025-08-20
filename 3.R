library(ggplot2)
library(dplyr)
library(tidyr)
# install.packages("vegan")
#library(vegan)
library(patchwork)

taxonomy <- read.csv("../outputs/for_standard/ASV_to_taxonomy.csv", header = TRUE, row.names = 1)
sequencing_data <- read.csv("../outputs/for_standard/ASV_counts_per_sample.csv", header = TRUE, row.names = 1)

spike_in_Imtechella <- subset(taxonomy,Family == "Flavobacteriaceae" & Genus == "Imtechella" & Species == "halotolerans") 
spike_in_Allobacillus <- subset(taxonomy,Family == "Bacillaceae" & Genus == "Allobacillus" & Species == "halotolerans") 

spike_in_asv_Imtechella <- rownames(spike_in_Imtechella)
spike_in_asv_Allobacillus <- rownames(spike_in_Allobacillus)
spike_in_asv_both <- c(spike_in_asv_Imtechella, spike_in_asv_Allobacillus)
print(spike_in_asv_both)

sequencing_data$spike_in_sum <- rowSums(sequencing_data[, spike_in_asv_both, drop = FALSE])
sequencing_data <- sequencing_data[, !colnames(sequencing_data) %in% spike_in_asv_both]
print(sequencing_data)
sequencing_data$spike_in_sum

spike_in_cell_count = 4*(10^7) # changing based on spike-in volume, this is for 20 ul
sequencing_data_normalized <- sequencing_data
other_cols <- colnames(sequencing_data)[colnames(sequencing_data) != "spike_in_sum"]
sequencing_data_normalized[, other_cols] <- spike_in_cell_count / sequencing_data$spike_in_sum * sequencing_data[, other_cols]
sequencing_data_normalized <- sequencing_data_normalized[, colnames(sequencing_data_normalized) != "spike_in_sum"]
print(sequencing_data_normalized)
write.csv(sequencing_data_normalized, "../outputs/for_soil/ASV_spikein_removed.csv", row.names = TRUE)

# normalize based on sample mass (optional)
mass_data <- read.csv("../RU_forest_sample_mass.csv")
asv_data_with_mass <- merge(sequencing_data_normalized, mass_data,
                           by.x = "row.names", by.y = "sample", all.x = TRUE)
print(asv_data_with_mass)
rownames(asv_data_with_mass) <- asv_data_with_mass$Row.names
asv_data_with_mass$Row.names <- NULL
print(asv_data_with_mass)
other_cols <- setdiff(colnames(asv_data_with_mass), c("mass.g."))
asv_data_normalized_by_mass <- asv_data_with_mass[, other_cols] / asv_data_with_mass$mass.g.
asv_data_normalized_by_mass$Sample <- rownames(asv_data_normalized_by_mass)
print(asv_data_normalized_by_mass)
write.csv(asv_data_normalized_by_mass, "normalized_sequencing_data.csv", row.names = TRUE)

other_cols <- setdiff(colnames(asv_data_normalized_by_mass), c("Sample"))
asv_data <- asv_data_normalized_by_mass[, other_cols]
asv_totals <- colSums(asv_data)
top_20_asvs <- names(sort(asv_totals, decreasing = TRUE)[1:20])
plot_data <- asv_data
plot_data$Sample <- rownames(plot_data)
plot_data$Other <- rowSums(plot_data[, !colnames(plot_data) %in% c(top_20_asvs, "Sample")])
plot_data <- plot_data[, c("Sample", top_20_asvs, "Other")]
plot_data_long <- pivot_longer(plot_data, 
                              cols = -Sample, 
                              names_to = "ASV", 
                              values_to = "Abundance")                           
taxonomy$ASV <- rownames(taxonomy)
tax_labels <- setNames(taxonomy$Genus, taxonomy$ASV) 
plot_data_long$Taxonomy <- ifelse(plot_data_long$ASV == "Other", 
                                  "Other", 
                                  tax_labels[plot_data_long$ASV])
plot_data_long$Taxonomy[is.na(plot_data_long$Taxonomy)] <- plot_data_long$ASV[is.na(plot_data_long$Taxonomy)]
top_20_taxonomy <- tax_labels[top_20_asvs]
taxonomy_counts <- table(top_20_taxonomy)
duplicated_taxa <- names(taxonomy_counts[taxonomy_counts > 1])
for (taxa in duplicated_taxa) {
  matching_asvs <- names(top_20_taxonomy[top_20_taxonomy == taxa])
  for (i in seq_along(matching_asvs)) {
    asv_name <- matching_asvs[i]
    plot_data_long$Taxonomy[plot_data_long$ASV == asv_name] <- 
      paste(taxa, paste0("(", asv_name, ")"))
  }
}
unique_taxonomy_levels <- unique(c("Other", plot_data_long$Taxonomy[plot_data_long$ASV %in% top_20_asvs]))
plot_data_long$Taxonomy <- factor(plot_data_long$Taxonomy, 
                                  levels = unique_taxonomy_levels)
colors <- c("lightgray", 
                  colorRampPalette(c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
  "#1FF8FF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#4B6A53", "#B249D5", "#7EDC45",
  "#5C47B8", "#CFD251", "#FF69B4", "#69C86C", "#CD3E50", "#83D5AF",
  "#DA6130", "#5E79B2", "#C29545", "#532A5A", "#5F7B35", "#C497CF",
  "#773A27", "#7CB9CB", "#594E50", "#D3C4A8", "#C17E7F"
))(20))
ggplot(plot_data_long, aes(x = Sample, y = Abundance, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Cell Counts by Sample", 
       x = "Sample", 
       y = "Absolute Abundance") +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),  # keep only this one
    legend.position = "right"
  ) +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(title = "Genus", ncol = 1))

