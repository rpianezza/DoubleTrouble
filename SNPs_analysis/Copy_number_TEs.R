library(tidyverse)
theme_set(theme_bw())

col_names <- c("run_accession","TE", "All_reads", "HQ_reads")

df_overview <- data.frame(
  Species = c("D. mel", "D. sim", "D. sec", "D. mau", "D. tei", "D. yak", "D. san", "D. ere", "D. ore"),
  Samples = numeric(9),
  Spoink_copies_haploid_genome = numeric(9),
  Shellder_copies_haploid_genome = numeric(9),
  P_element_copies_haploid_genome = numeric(9)
)

# Drosophila melanogaster

df_mel <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dmel.csv", header = FALSE)
names(df_mel) <- col_names

df_mel <- df_mel %>%
  filter(TE != "TE")

df_overview$Samples[1] <- df_mel %>%
  filter(TE == "spoink") %>%
  count()

df_overview$Spoink_copies_haploid_genome[1] <- df_mel %>%
  filter(TE == "spoink") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[1] <- df_mel %>%
  filter(TE == "PPI251") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)


# Drosophila simulans

df_sim <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dsim.csv", header = FALSE)
names(df_sim) <- col_names

df_sim <- df_sim %>%
  filter(TE != "TE")

df_overview$Samples[2] <- df_sim %>%
  filter(TE == "Shellder") %>%
  count()

df_overview$Spoink_copies_haploid_genome[2] <- df_sim %>%
  filter(TE == "spoink") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[2] <- df_sim %>%
  filter(TE == "Shellder") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[2] <- df_sim %>%
  filter(TE == "PPI251") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

# Drosophila sechelia

df_sec <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dsec.csv", header = FALSE)
names(df_sec) <- col_names

df_sec <- df_sec %>%
  filter(TE != "TE")


df_overview$Samples[3] <- df_sec %>%
  filter(TE == "Shellder") %>%
  count()

df_overview$Spoink_copies_haploid_genome[3] <- df_sec %>%
  filter(TE == "spoink") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[3] <- df_sec %>%
  filter(TE == "Shellder") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[3] <- df_sec %>%
  filter(TE == "PPI251") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

# Drosophila mauritiana

df_mau <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dmau.csv", header = FALSE)
names(df_mau) <- col_names

df_mau <- df_mau %>%
  filter(TE != "TE")

df_overview$Samples[4] <- df_mau %>%
  filter(TE == "Shellder") %>%
  count()

df_overview$Spoink_copies_haploid_genome[4] <- df_mau %>%
  filter(TE == "spoink") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[4] <- df_mau %>%
  filter(TE == "Shellder") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[4] <- df_mau %>%
  filter(TE == "PPI251") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)


# Drosophila teissieri
df_tei <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dtei.csv", header = FALSE)
names(df_tei) <- col_names

df_tei <- df_tei %>%
  mutate(TE = ifelse(TE == "gypsy-29-dsim", "Shellder", ifelse(TE == "gypsy-7-sim1", "Spoink", TE)))


df_tei <- df_tei %>%
  filter(TE != "TE")

df_overview$Samples[5] <- df_tei %>%
  filter(TE == "Shellder") %>%
  count()

df_overview$Spoink_copies_haploid_genome[5] <- df_tei %>%
  filter(TE == "Spoink") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[5] <- df_tei %>%
  filter(TE == "Shellder") %>%
  filter(HQ_reads > 1 ) %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)


# Drosophila yakuba

df_yak_Santome <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/Dyak_Santome.csv", header = FALSE)
df_yak_Africa <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/Dyak_Africa.csv", header = FALSE)

names(df_yak_Santome) <- col_names
names(df_yak_Africa) <- col_names

df_yak <- rbind(df_yak_Santome, df_yak_Africa)

df_yak <- df_yak %>%
  filter(TE != "TE")

df_overview$Samples[6] <- df_yak %>%
  filter(TE == "Shellder") %>%
  count()

df_overview$Spoink_copies_haploid_genome[6] <- df_yak %>%
  filter(TE == "Spoink") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[6] <- df_yak %>%
  filter(TE == "Shellder") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[6] <- df_yak %>%
  filter(TE == "PPI251") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)


# Drosophila santomea
df_san <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dsan.csv", header = FALSE)
names(df_san) <- col_names


df_san <- df_san %>%
  filter(TE != "TE")

df_overview$Samples[7] <- 17

df_overview$Spoink_copies_haploid_genome[7] <- df_san %>%
  filter(TE == "Spoink") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[7] <- df_san %>%
  filter(TE == "Shellder") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[7] <- df_san %>%
  filter(TE == "PPI251") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)


# Drosophila erecta
df_ere <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dere.csv", header = FALSE)
names(df_ere) <- col_names


df_ere <- df_ere %>%
  filter(TE != "TE")

df_overview$Samples[8] <- df_ere %>%
  filter(TE == "Shellder") %>%
  count()/2

df_overview$Spoink_copies_haploid_genome[8] <- df_ere %>%
  filter(TE == "Spoink") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$Shellder_copies_haploid_genome[8] <- df_ere %>%
  filter(TE == "Shellder") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)

df_overview$P_element_copies_haploid_genome[8] <- df_ere %>%
  filter(TE == "PPI251") %>%
  summarise(mean_HQ_reads = mean(as.numeric(HQ_reads), na.rm = TRUE)) %>%
  pull(mean_HQ_reads) %>%
  round(digits = 0)


#D ore 1 sample, no Spoink, no Shellder, no P-element
df_ore <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/csv/dore.csv", header = FALSE)
names(df_ore) <- col_names


df_overview$Samples[9] <- df_ore %>%
  filter(TE == "Shellder") %>%
  count()
