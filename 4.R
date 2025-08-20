# accessing spike-in bias in my soil samples

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

taxonomy <- read.csv("../outputs/for_soil/ASV_to_taxonomy.csv", header = TRUE, row.names = 1)
sequencing_data <- read.csv("../outputs/for_soil/ASV_counts_per_sample.csv", header = TRUE, row.names = 1)
spike_in_Imtechella <- subset(taxonomy,Family == "Flavobacteriaceae" & Genus == "Imtechella" & Species == "halotolerans") 
spike_in_Allobacillus <- subset(taxonomy,Family == "Bacillaceae" & Genus == "Allobacillus" & Species == "halotolerans") 

spike_in_asv_Imtechella <- rownames(spike_in_Imtechella)
spike_in_asv_Allobacillus <- rownames(spike_in_Allobacillus)
sequencing_data$Allobacillus<- rowSums(sequencing_data[, spike_in_asv_Allobacillus, drop = FALSE])
sequencing_data <- sequencing_data[, !colnames(sequencing_data) %in% spike_in_asv_Allobacillus]
print(sequencing_data$Allobacillus)

sequencing_data$Imtechella<- rowSums(sequencing_data[, spike_in_asv_Imtechella, drop = FALSE])
sequencing_data <- sequencing_data[, !colnames(sequencing_data) %in% spike_in_asv_Imtechella]
print(sequencing_data$Imtechella)

sample_cols <- setdiff(colnames(sequencing_data), c("Imtechella","Allobacillus"))
sequencing_data$sample_asv<- rowSums(sequencing_data[, sample_cols, drop = FALSE])
sequencing_data <- sequencing_data[, !colnames(sequencing_data) %in% sample_cols]
sequencing_data$sum<- rowSums(sequencing_data)
print(sequencing_data)

Allobacillus_absolute = 1.4*(10^8)
Imtechella_absolute = 6*(10^7)
sequencing_data$Allobacillus_absolute <- Allobacillus_absolute
sequencing_data$Imtechella_absolute <- Imtechella_absolute

sequencing_data$Allobacillus_fraction <- sequencing_data$Allobacillus / sequencing_data$sum 
sequencing_data$Imtechella_fraction <- sequencing_data$Imtechella / sequencing_data$sum 
sequencing_data <- sequencing_data[!rownames(sequencing_data) %in% c("a2_B_S10", "a2_B_SI_S11"), ]
sequencing_data$fractions <- sequencing_data$Imtechella / sequencing_data$Allobacillus
AL = sum(sequencing_data$Allobacillus)
IM = sum(sequencing_data$Imtechella)
x = IM/AL
print(sequencing_data)
write.csv(sequencing_data, "../outputs/spike_in_comparison.csv", row.names = TRUE)

sample_metadata <- data.frame(
  sample = c("a2_H1_S1", "a2_H2_S2", "a2_H3_S3",
             "a2_M1_S4", "a2_M2_S5", "a2_M3_S6",
             "a2_L1_S7", "a2_L2_S8", "a2_L3_S9"),
  treatment = c("Site H", "Site H", "Site H",
                "Site M", "Site M", "Site M",
                "Site L", "Site L", "Site L")
)
plot_data <- sequencing_data %>%
  mutate(sample = rownames(sequencing_data)) %>%
  left_join(sample_metadata, by = "sample") %>%
  pivot_longer(
    cols = c(Allobacillus_absolute, Imtechella_absolute, 
             Allobacillus_fraction, Imtechella_fraction),
    names_to = "measure",
    values_to = "value"
  ) %>%
  separate(measure, into = c("species", "type"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  filter(!is.na(absolute), !is.na(fraction))

ggplot(plot_data, aes(x = absolute, y = fraction, color = treatment, shape = species)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(aes(group = sample), alpha = 0.5, color = "black") +
  # Add the theoretical red line
  geom_segment(x = 6e7, y = 0.15, xend = 1.4e8, yend = 0.348, 
               color = "red", size = 1.2, linetype = "dashed") +
  # Add text annotation
  annotate("text", x = 1.2e8, y = 0.28, label = "Theoretical Slope", 
           color = "red", size = 5, angle = 38, hjust = 0.5) +
  theme_bw() +
  labs(
    x = "Theoretical Spike-in 16S Reads Added",
    y = "Actual Reads Fraction",
    color = "Site",
    shape = "Species"
  ) +
  scale_x_continuous(labels = scales::scientific) + theme(
    axis.text = element_text(size = 12),        
    axis.title = element_text(size = 16),       
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 14, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold"),    
    legend.position = "right"
  )






# using correlation to compare community composition similarities

df <- read.csv("../outputs/for_soil/ASV_counts_per_sample.csv", header = TRUE, row.names = 1)
df <- df[, !colnames(df) %in% spike_in_asv_both]
df <- df[, colSums(df) >= 10]
df <- sweep(df, 1, rowSums(df), FUN = "/")
df1 <- df[rownames(df) %in% c("a2_H1_S1", "a2_H2_S2","a2_H3_S3"), ]
print(df1)

samples <- rownames(df1)
plots <- list()
make_plot <- function(x_raw, y_raw, x_label, y_label, n_bootstrap = 1000) {
  mask <- (x_raw > 0) | (y_raw > 0)
  x <- log10(x_raw + 1e-6)
  y <- log10(y_raw + 1e-6)
  cor_result <- cor.test(x, y, method = "pearson")
  r <- cor_result$estimate
  p <- cor_result$p.value
  n <- length(x)
  r_boot <- replicate(n_bootstrap, {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    cor(x[idx], y[idx])
  })
  ci <- quantile(r_boot, probs = c(0.025, 0.975))
  
  ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      x = paste0(x_label),
      y = paste0(y_label),
      title = paste0(
        "r = ", round(r, 3),
        ", p = ", signif(p, 3),
        "\n95% CI: [", round(ci[1], 3), ", ", round(ci[2], 3), "]"
      )
    ) +
    theme_minimal() +
    coord_fixed()
}
plots[[1]] <- make_plot(
  x_raw = as.numeric(df[samples[1], ]),
  y_raw = as.numeric(df[samples[2], ]),
  x_label = "H Site Rep1",
  y_label = "H Site Rep2"
)
plots[[2]] <- make_plot(
  x_raw = as.numeric(df[samples[1], ]),
  y_raw = as.numeric(df[samples[3], ]),
  x_label = "H Site Rep1",
  y_label = "H Site Rep3"
)
plots[[3]] <- make_plot(
  x_raw = as.numeric(df[samples[2], ]),
  y_raw = as.numeric(df[samples[3], ]),
  x_label = "H Site Rep2",
  y_label = "H Site Rep3"
)

grid.arrange(grobs = plots, ncol = 3)\df4 <- df[rownames(df) %in% c("a2_H3_S3","a2_M3_S6","a2_L3_S9"),]
df4 <- df4[, colSums(df4) >= 10]
print(df4)

samples <- rownames(df4)
plots <- list()
make_plot <- function(x_raw, y_raw, x_label, y_label, n_bootstrap = 1000) {
  mask <- (x_raw > 0) | (y_raw > 0)
  x <- log10(x_raw + 1e-6)
  y <- log10(y_raw + 1e-6)
  cor_result <- cor.test(x, y, method = "pearson")
  r <- cor_result$estimate
  p <- cor_result$p.value
  n <- length(x)
  r_boot <- replicate(n_bootstrap, {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    cor(x[idx], y[idx])
  })
  ci <- quantile(r_boot, probs = c(0.025, 0.975))
  
  ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      x = paste0(x_label),
      y = paste0(y_label),
      title = paste0(
        "r = ", round(r, 3),
        ", p = ", signif(p, 3),
        "\n95% CI: [", round(ci[1], 3), ", ", round(ci[2], 3), "]"
      )
    ) +
    theme_minimal() +
    coord_fixed()
}
plots[[1]] <- make_plot(
  x_raw = as.numeric(df[samples[1], ]),
  y_raw = as.numeric(df[samples[2], ]),
  x_label = "H Site Rep1",
  y_label = "M Site Rep1"
)
plots[[2]] <- make_plot(
  x_raw = as.numeric(df[samples[1], ]),
  y_raw = as.numeric(df[samples[3], ]),
  x_label = "H Site Rep1",
  y_label = "L Site Rep1"
)
plots[[3]] <- make_plot(
  x_raw = as.numeric(df[samples[2], ]),
  y_raw = as.numeric(df[samples[3], ]),
  x_label = "M Site Rep1",
  y_label = "L Site Rep1"
)
grid.arrange(grobs = plots, ncol = 3)






# accessing biases in community standard species

taxa_comp <- read.csv("../outputs/for_standard/ASV_counts_per_sample.csv", header = TRUE, row.names = 1)
spike_in_Imtechella <- subset(taxonomy,Family == "Flavobacteriaceae" & Genus == "Imtechella" & Species == "halotolerans") 
spike_in_Allobacillus <- subset(taxonomy,Family == "Bacillaceae" & Genus == "Allobacillus" & Species == "halotolerans")
spike_in_asv_both <- c(spike_in_asv_Imtechella, spike_in_asv_Allobacillus)
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% spike_in_asv_both]

taxa_comp$tot_reads <- rowSums(taxa_comp)

Pseudomonas <- subset(taxonomy,Genus == "Pseudomonas") 
Pseudomonas <- rownames(Pseudomonas)
taxa_comp$Pseudomonas_sum <- rowSums(sequencing_data[, Pseudomonas, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Pseudomonas]

Escherichia <- subset(taxonomy,Genus == "Escherichia-Shigella") 
Escherichia <- rownames(Escherichia)
taxa_comp$Escherichia_sum <- rowSums(sequencing_data[, Escherichia, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Escherichia]
print(taxa_comp$Escherichia_sum)

Salmonella <- subset(taxonomy,Genus == "Salmonella") 
Salmonella <- rownames(Salmonella)
taxa_comp$Salmonella_sum <- rowSums(sequencing_data[, Salmonella, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Salmonella]

Limosilactobacillus <- subset(taxonomy,Genus == "Limosilactobacillus") 
Limosilactobacillus <- rownames(Limosilactobacillus)
taxa_comp$Limosilactobacillus_sum <- rowSums(sequencing_data[, Limosilactobacillus, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Limosilactobacillus]

Staphylococcus <- subset(taxonomy,Genus == "Staphylococcus") 
Staphylococcus <- rownames(Staphylococcus)
taxa_comp$Staphylococcus_sum <- rowSums(sequencing_data[, Staphylococcus, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Staphylococcus]

Listeria <- subset(taxonomy,Genus == "Listeria") 
Listeria <- rownames(Listeria)
taxa_comp$Listeria_sum <- rowSums(sequencing_data[, Listeria, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Listeria]

Bacillus <- subset(taxonomy,Genus == "Bacillus") 
Bacillus <- rownames(Bacillus)
taxa_comp$Bacillus_sum <- rowSums(sequencing_data[, Bacillus, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Bacillus]

Enterococcus <- subset(taxonomy,Genus == "Enterococcus") 
Enterococcus <- rownames(Enterococcus)
taxa_comp$Enterococcus_sum <- rowSums(sequencing_data[, Enterococcus, drop = FALSE])
taxa_comp <- taxa_comp[, !colnames(taxa_comp) %in% Enterococcus]

taxa_comp <- taxa_comp %>% select(-contains("ASV"))
taxa_comp <- taxa_comp[!rownames(taxa_comp) %in% c("Q_B_S7", "Z_B_S14"), ]                       
print(taxa_comp)

df <- read.csv("../outputs/taxa_comp.csv", header = TRUE, row.names = 1)

# Add genus names as a column
df$Genus <- rownames(df)

# Reshape data from wide to long format
df_long <- df %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Value")

# Create sample groups for better color coding
df_long$Sample_Group <- case_when(
  grepl("Q_2_", df_long$Sample) ~ "Qiagen 2-Fold",
  grepl("Q_20_", df_long$Sample) ~ "Qiagen 20-Fold",
  grepl("Q_200_", df_long$Sample) ~ "Qiagen 200-Fold",
  grepl("Z_2_", df_long$Sample) ~ "Zymo 2-Fold",
  grepl("Z_20_", df_long$Sample) ~ "Zymo 20-Fold",
  grepl("Z_200_", df_long$Sample) ~ "Zymo 200-Fold"
)

# Create short labels for points
df_long$Short_Label <- case_when(
  df_long$Sample == "Q_2_1_S1" ~ "Q2-1",
  df_long$Sample == "Q_2_2_S2" ~ "Q2-2",
  df_long$Sample == "Q_20_1_S3" ~ "Q20-1",
  df_long$Sample == "Q_20_2_S4" ~ "Q20-2",
  df_long$Sample == "Q_200_1_S5" ~ "Q200-1",
  df_long$Sample == "Q_200_2_S6" ~ "Q200-2",
  df_long$Sample == "Z_2_1_S8" ~ "Z2-1",
  df_long$Sample == "Z_2_2_S9" ~ "Z2-2",
  df_long$Sample == "Z_20_1_S10" ~ "Z20-1",
  df_long$Sample == "Z_20_2_S11" ~ "Z20-2",
  df_long$Sample == "Z_200_1_S12" ~ "Z200-1",
  df_long$Sample == "Z_200_2_S13" ~ "Z200-2"
)

df_long$Genus_Clean <- gsub("_sum", "", df_long$Genus)


df_long$Gram_Status <- case_when(
  df_long$Genus_Clean %in% c("Pseudomonas", "Escherichia", "Salmonella") ~ "Gram-Negative",
  TRUE ~ "Gram-Positive"
)

df_long$Genus_Clean <- factor(
  df_long$Genus_Clean,
  levels = c(
    unique(df_long$Genus_Clean[df_long$Gram_Status == "Gram-Negative"]),
    unique(df_long$Genus_Clean[df_long$Gram_Status == "Gram-Positive"])
  )
)

p1 <- ggplot(df_long, aes(x = Genus_Clean, y = Value)) +
  geom_boxplot(aes(fill = Gram_Status), alpha = 0.6, outlier.shape = NA) +
  geom_point(aes(color = Sample_Group), size = 3, alpha = 0.8, 
             position = position_jitter(width = 0.2, seed = 123)) +
  labs(
    title = "Distribution of Values by Genus",
    x = "Genus",
    y = "Fold Diviation From Expected",
    color = "Sample Group",
    fill = "Gram Stain"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.title = element_text(size = 11),
    legend.position = "right"
  ) +
  scale_color_brewer(type = "qual", palette = "Set2") +
  scale_fill_manual(values = c("Gram-Negative" = "#4F81BD", "Gram-Positive" = "#C0504D")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7)
print(p1)



# comparing theoretical and actual spike-in ratios in soil and standard samples

df <- read.csv("../outputs/spike_in_comparison.csv", header = TRUE)
df$Group <- case_when(
  grepl("Qiagen.*2-Fold", df$Condition) ~ "Qiagen 2-Fold Dilution",
  grepl("Qiagen.*20-Fold", df$Condition) ~ "Qiagen 20-Fold Dilution", 
  grepl("Qiagen.*200-Fold", df$Condition) ~ "Qiagen 200-Fold Dilution",
  grepl("Zymo.*2-Fold", df$Condition) ~ "Zymo 2-Fold Dilution",
  grepl("Zymo.*20-Fold", df$Condition) ~ "Zymo 20-Fold Dilution",
  grepl("Zymo.*200-Fold", df$Condition) ~ "Zymo 200-Fold Dilution",
  grepl("Soil", df$Condition) ~ "Soil"
)
df$Group <- factor(df$Group, levels = c(
  "Qiagen 2-Fold Dilution", "Qiagen 20-Fold Dilution", "Qiagen 200-Fold Dilution",
  "Zymo 2-Fold Dilution", "Zymo 20-Fold Dilution", "Zymo 200-Fold Dilution",
  "Soil"
))

ggplot(df, aes(x = Actual, y = Theoretical, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_hline(yintercept = 2.33, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_point(x = 2.33, y = 2.33, color = "red", size = 5, shape = 16) +
  geom_text(x = 2.33, y = 2.33, label = "Theoretical", 
            vjust = -0.8, hjust = 0.5, color = "red", size = 4, fontface = "bold") +
  labs(
    title = "Allobacillus to Imtechella Ratio",
    x = "Actual Ratio",
    y = "Theoretical Ratio",
    color = "Experimental Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 9)
  ) +
  scale_color_brewer(type = "qual", palette = "Set3") +
  guides(color = guide_legend(override.aes = list(size = 3)))


