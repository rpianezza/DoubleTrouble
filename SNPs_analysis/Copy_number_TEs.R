library(tidyverse)
theme_set(theme_bw())

df_mel <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dmel.csv", header = TRUE)
df_sim <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dsim.csv", header = TRUE)
df_sec <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dsec.csv", header = TRUE)
df_mau <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dmau.csv", header = TRUE)
df_tei <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dtei.csv", header = TRUE)
df_yak <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dyak.csv", header = TRUE)
df_san <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dsan.csv", header = TRUE)
df_ere <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dere.csv", header = TRUE)
df_ore <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dore.csv", header = TRUE)

df_mel <- df_mel %>%filter(TE != "TE")
df_sim <- df_sim %>%filter(TE != "TE")
df_sec <- df_sec %>%filter(TE != "TE")
df_mau <- df_mau %>%filter(TE != "TE")
df_tei <- df_tei %>%filter(TE != "TE")
df_yak <- df_yak %>%filter(TE != "TE")
df_san <- df_san %>%filter(TE != "TE")
df_ere <- df_ere %>%filter(TE != "TE")
df_ore <- df_ore %>%filter(TE != "TE")

df_overview <- data.frame(
  Species = c("D. mel", "D. sim", "D. sec", "D. mau", "D. tei", "D. yak", "D. san", "D. ere", "D. ore"),
  Samples = numeric(9),
  Shellder_positive_samples = numeric(9),
  spoink_positive_samples = numeric(9),
  PPI251_positive_samples = numeric(9),
  Shellder_copies_haploid_genome = numeric(9),
  spoink_copies_haploid_genome = numeric(9),
  PPI251_copies_haploid_genome = numeric(9)
)

mean_HQ_reads <- function(data_frame, TE_name) {
  mean_HQ_reads <- data_frame %>%
    filter(TE == TE_name) %>%
    filter(HQ_reads > 1) %>%
    summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
    pull(mean_HQ_reads)
  
  # If mean_HQ_reads is NaN, set it to 0
  if (is.nan(mean_HQ_reads)) {
    mean_HQ_reads <- 0
  } else {
    mean_HQ_reads <- round(mean_HQ_reads, digits = 0)
  }
  
  return(mean_HQ_reads)
}

df_overview$Samples[1] <- df_mel %>% filter(TE == "PPI251") %>% count()
df_overview$Samples[2] <- df_sim %>% filter(TE == "PPI251") %>% count()
df_overview$Samples[3] <- df_sec %>% filter(TE == "Shellder") %>% count()
df_overview$Samples[4] <- df_mau %>% filter(TE == "Shellder") %>% count()
df_overview$Samples[5] <- df_tei %>% filter(TE == "Shellder") %>% count()
df_overview$Samples[6] <- df_yak %>% filter(TE == "Shellder") %>% count()
df_overview$Samples[7] <- df_san %>% filter(TE == "Shellder") %>% count()-6
df_overview$Samples[8] <- df_ere %>% filter(TE == "Shellder") %>% count()/2
df_overview$Samples[9] <- df_ore %>% filter(TE == "Shellder") %>% count()


TE_n = "Shellder"

df_overview$Shellder_copies_haploid_genome[1] <- mean_HQ_reads(df_mel, TE_n)
df_overview$Shellder_copies_haploid_genome[2] <- mean_HQ_reads(df_sim, TE_n)
df_overview$Shellder_copies_haploid_genome[3] <- mean_HQ_reads(df_sec, TE_n)
df_overview$Shellder_copies_haploid_genome[4] <- mean_HQ_reads(df_mau, TE_n)
df_overview$Shellder_copies_haploid_genome[5] <- mean_HQ_reads(df_tei, TE_n)
df_overview$Shellder_copies_haploid_genome[6] <- mean_HQ_reads(df_yak, TE_n)
df_overview$Shellder_copies_haploid_genome[7] <- mean_HQ_reads(df_san, TE_n)
df_overview$Shellder_copies_haploid_genome[8] <- mean_HQ_reads(df_ere, TE_n)
df_overview$Shellder_copies_haploid_genome[9] <- mean_HQ_reads(df_ore, TE_n)

df_overview$Shellder_positive_samples[1] <- df_mel %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[2] <- df_sim %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[3] <- df_sec %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[4] <- df_mau %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[5] <- df_tei %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[6] <- df_yak %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[7] <- df_san %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[8] <- df_ere %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$Shellder_positive_samples[9] <- df_ore %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()


TE_n = "spoink"

df_overview$spoink_copies_haploid_genome[1] <- mean_HQ_reads(df_mel, TE_n)
df_overview$spoink_copies_haploid_genome[2] <- mean_HQ_reads(df_sim, TE_n)
df_overview$spoink_copies_haploid_genome[3] <- mean_HQ_reads(df_sec, TE_n)
df_overview$spoink_copies_haploid_genome[4] <- mean_HQ_reads(df_mau, TE_n)
df_overview$spoink_copies_haploid_genome[5] <- mean_HQ_reads(df_tei, TE_n)
df_overview$spoink_copies_haploid_genome[6] <- mean_HQ_reads(df_yak, TE_n)
df_overview$spoink_copies_haploid_genome[7] <- mean_HQ_reads(df_san, TE_n)
df_overview$spoink_copies_haploid_genome[8] <- mean_HQ_reads(df_ere, TE_n)
df_overview$spoink_copies_haploid_genome[9] <- mean_HQ_reads(df_ore, TE_n)

df_overview$spoink_positive_samples[1] <- df_mel %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[2] <- df_sim %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[3] <- df_sec %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[4] <- df_mau %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[5] <- df_tei %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[6] <- df_yak %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[7] <- df_san %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[8] <- df_ere %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$spoink_positive_samples[9] <- df_ore %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()


TE_n = "PPI251"

df_overview$PPI251_copies_haploid_genome[1] <- mean_HQ_reads(df_mel, TE_n)
df_overview$PPI251_copies_haploid_genome[2] <- mean_HQ_reads(df_sim, TE_n)
df_overview$PPI251_copies_haploid_genome[3] <- mean_HQ_reads(df_sec, TE_n)
df_overview$PPI251_copies_haploid_genome[4] <- mean_HQ_reads(df_mau, TE_n)
df_overview$PPI251_copies_haploid_genome[5] <- mean_HQ_reads(df_tei, TE_n)
df_overview$PPI251_copies_haploid_genome[6] <- mean_HQ_reads(df_yak, TE_n)
df_overview$PPI251_copies_haploid_genome[7] <- mean_HQ_reads(df_san, TE_n)
df_overview$PPI251_copies_haploid_genome[8] <- mean_HQ_reads(df_ere, TE_n)
df_overview$PPI251_copies_haploid_genome[9] <- mean_HQ_reads(df_ore, TE_n)

df_overview$PPI251_positive_samples[1] <- df_mel %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[2] <- df_sim %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[3] <- df_sec %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[4] <- df_mau %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[5] <- df_tei %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[6] <- df_yak %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[7] <- df_san %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[8] <- df_ere %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
df_overview$PPI251_positive_samples[9] <- df_ore %>% filter(TE == TE_n) %>% filter(HQ_reads > 1) %>% count()
