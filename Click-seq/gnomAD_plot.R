# Load libraries
library(ggplot2)
library(scales)
library(dplyr)
library(stringr)

# Define your working directory
path = "/Users/gabrielleferra/Desktop/CYP-project/20250828_clickseq_dataset_MINQUAL67_final_datasets_and_analysis/"

# Load the final click-seq DMS activity scores
CYP2D6_data <- read.table(file = file.path(path, "Final_Scores/Pacybara_CYP2D6_activity_scores_singles_fs_and_doubles_12_1e-5_5681e-5_2expts_softfilter_MINQUAL67.csv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load the gnomAD data
gnomAD_data <- read.table(file = "/Users/gabrielleferra/Desktop/CYP-project/gnomAD_analysis/gnomAD_singlesitesub_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Add variant column to gnomAD_data to match CYP2D6_data input
gnomAD_data <- gnomAD_data %>%
  mutate(
    variant = str_replace(Protein.Consequence, "^p\\.", ""),         # remove leading 'p.'
  ) %>%
  mutate(
    variant = str_replace_all(variant, 
                              c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", 
                                "Cys"="C", "Gln"="Q", "Glu"="E", "Gly"="G", 
                                "His"="H", "Ile"="I", "Leu"="L", "Lys"="K", 
                                "Met"="M", "Phe"="F", "Pro"="P", "Ser"="S", 
                                "Thr"="T", "Trp"="W", "Tyr"="Y", "Val"="V"))
  )

df_merge <- merge(CYP2D6_data, gnomAD_data, by = "variant") 

df_merge$activity_class <- factor(df_merge$activity_class, 
                                  levels = c('possibly_increased', 'increased', 'possibly_wt-like', 'wt-like', 'possibly_decreased', 'decreased', 'possibly_nonsense-like', 'nonsense-like'))

top_3_rows <- df_merge %>%
  arrange(desc(Allele.Frequency)) %>%  # Arrange the dataframe by descending Allele.Frequency
  head(3)  # Select the top 3 rows

# Combine AFs and keep positive finite values (needed for log scale)
af_all <- c(df_merge$Allele.Frequency, top_3_rows$Allele.Frequency)
af_all <- af_all[is.finite(af_all) & !is.na(af_all) & af_all > 0]

min_af <- min(af_all)
max_af <- max(af_all)

# Rightmost must cover ~0.55 -> use 1
rightmost_tick <- 1

# Decade ticks only
left_decade_exp <- ceiling(log10(min_af))          # first decade >= min
decade_breaks   <- 10^(seq.int(from = left_decade_exp, to = 0, by = 1))

gnomAD_plot <- ggplot(df_merge, aes(x = score, y = Allele.Frequency, color = activity_class)) +
  geom_point(size = 3.5) +
  geom_text(
    data = subset(top_3_rows, is.finite(Allele.Frequency) & Allele.Frequency > 0),
    aes(x = score, y = Allele.Frequency, label = variant),
    fontface = "bold",
    vjust = -1,
    hjust = 0.5,
    size = 5,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_y_log10(   # <-- move log scale to y because y is allele frequency now
    limits = c(min_af, rightmost_tick),
    breaks = decade_breaks,
    labels = function(x) {
      out <- label_scientific(digits = 1)(x)
      out[x == 1] <- "1"
      out
    },
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(
    values = c('#28C765','#B1D7FA','#3888D1','#FFE600','#FFAB00','#FF6E6E','#D53E4F'),
    labels = c("Possibly Increased","Possibly WT-like","WT-like",
               "Possibly Decreased","Decreased",
               "Possibly Nonsense-like","Nonsense-like"),
    name   = "Click-seq Interpretation"
  ) +
  labs(
    title = "Click-seq Activity Scores vs. gnomAD Allele Frequency",
    x = "Click-seq Activity Score",
    y = "gnomAD allele frequency"
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text    = element_text(size = 15),
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 16),
    legend.text  = element_text(size = 14)
  )

# Print the plot
print(gnomAD_plot)

# I saved the plots manually but feel free to change the code here
