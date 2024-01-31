library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

df_Shellder_dsim <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/Shellder-dsim.tsv", header = TRUE, sep = "\t")
colnames(df_Shellder_dsim) <- sub("^Shellder_", "", colnames(df_Shellder_dsim))
v_Shellder_dsim <- names(df_Shellder_dsim)[-1]
v_Shellder_dsim <- as.numeric(v_Shellder_dsim)

df_Shellder_dsec <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/Shellder-dsec.tsv", header = TRUE, sep = "\t")
colnames(df_Shellder_dsec) <- sub("^Shellder_", "", colnames(df_Shellder_dsec))
v_Shellder_dsec <- names(df_Shellder_dsec)[-1]
v_Shellder_dsec <- as.numeric(v_Shellder_dsec)

df_Shellder_dmau <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/Shellder-dmau.tsv", header = TRUE, sep = "\t")
colnames(df_Shellder_dmau) <- sub("^Shellder_", "", colnames(df_Shellder_dmau))
v_Shellder_dmau <- names(df_Shellder_dmau)[-1]
v_Shellder_dmau <- as.numeric(v_Shellder_dmau)

df_Shellder_dtei <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/Shellder-dtei.tsv", header = TRUE, sep = "\t")
colnames(df_Shellder_dtei) <- sub("^Shellder_", "", colnames(df_Shellder_dtei))
v_Shellder_dtei <- names(df_Shellder_dtei)[-1]
v_Shellder_dtei <- as.numeric(v_Shellder_dtei)

unique_Shellder <- sort(unique(c(v_Shellder_dsim, v_Shellder_dsec, v_Shellder_dmau, v_Shellder_dtei)))
common_Shellder <- intersect(intersect(intersect(v_Shellder_dsim, v_Shellder_dsec), v_Shellder_dmau), v_Shellder_dtei)


# Create a list of vectors
vectors <- list(v_Shellder_dsim = v_Shellder_dsim, v_Shellder_dsec = v_Shellder_dsec, v_Shellder_dmau = v_Shellder_dmau, v_Shellder_dtei = v_Shellder_dtei)

# Initialize empty lists to store data for the data frame
pos_list <- c()
species_list <- character()

# Iterate over the list items and populate the lists
for (species in names(vectors)) {
  pos_list <- c(pos_list, vectors[[species]])
  species_list <- c(species_list, rep(species, length(vectors[[species]])))
}

df_Shellder <- data.frame(pos = pos_list, species = species_list)


g_Shellder_old<-ggplot(data=df_Shellder,aes(x=pos, fill=species))+geom_histogram(binwidth=1)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  xlab("position")+ylab("SNPs")+
  ggtitle("Shellder's SNPs in different species")+
  scale_fill_manual(values = c("v_Shellder_dsim" = "#009900", "v_Shellder_dmau" = "#0033cc", "v_Shellder_dsec" = "#FF9900", "v_Shellder_dtei" = "#ff6666"),
                    labels = c("v_Shellder_dsim" = "D. simulans",
                               "v_Shellder_dsec" = "D. sechelia",
                               "v_Shellder_dmau" = "D. mauritiana",
                               "v_Shellder_dtei" = "D. teissieri"))

plot(g_Shellder_old)


g_Shellder<-ggplot(data=df_Shellder,aes(x=pos,y=species, color=species))+geom_point(alpha = 0.4)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  xlab("position")+ylab("SNPs")

plot(g_Shellder)



df_spoink_dsim <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/spoink-dsim.tsv", header = TRUE, sep = "\t")
colnames(df_spoink_dsim) <- sub("^spoink_", "", colnames(df_spoink_dsim))
v_spoink_dsim <- names(df_spoink_dsim)[-1]
v_spoink_dsim <- as.numeric(v_spoink_dsim)

df_spoink_dsec <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/spoink-dsec.tsv", header = TRUE, sep = "\t")
colnames(df_spoink_dsec) <- sub("^spoink_", "", colnames(df_spoink_dsec))
v_spoink_dsec <- names(df_spoink_dsec)[-1]
v_spoink_dsec <- as.numeric(v_spoink_dsec)

df_spoink_dmau <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/spoink-dmau.tsv", header = TRUE, sep = "\t")
colnames(df_spoink_dmau) <- sub("^spoink_", "", colnames(df_spoink_dmau))
v_spoink_dmau <- names(df_spoink_dmau)[-1]
v_spoink_dmau <- as.numeric(v_spoink_dmau)

df_spoink_dmel <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/SNPs/spoink-dmel.tsv", header = TRUE, sep = "\t")
colnames(df_spoink_dmel) <- sub("^spoink_", "", colnames(df_spoink_dmel))
v_spoink_dmel <- names(df_spoink_dmel)[-1]
v_spoink_dmel <- as.numeric(v_spoink_dmel)

unique_spoink <- sort(unique(c(v_spoink_dsim, v_spoink_dsec, v_spoink_dmau, v_spoink_dmel)))
common_spoink <- intersect(intersect(intersect(v_spoink_dsim, v_spoink_dsec), v_spoink_dmau), v_spoink_dmel)

# Create a list of vectors
vectors <- list(v_spoink_dsim = v_spoink_dsim, v_spoink_dsec = v_spoink_dsec, v_spoink_dmau = v_spoink_dmau, v_spoink_dmel = v_spoink_dmel)

# Initialize empty lists to store data for the data frame
pos_list <- c()
species_list <- character()

# Iterate over the list items and populate the lists
for (species in names(vectors)) {
  pos_list <- c(pos_list, vectors[[species]])
  species_list <- c(species_list, rep(species, length(vectors[[species]])))
}

df_spoink <- data.frame(pos = pos_list, species = species_list)


g_spoink_old<-ggplot(data=df_spoink,aes(x=pos, fill=species))+geom_histogram(binwidth=1)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  xlab("position")+ylab("SNPs")+
  ggtitle("Spoink's SNPs in different species")+
  scale_fill_manual(values = c("v_spoink_dsim" = "#009900", "v_spoink_dmau" = "#0033cc", "v_spoink_dsec" = "#FF9900", "v_spoink_dmel" = "#CC0033"),
                    labels = c("v_spoink_dsim" = "D. simulans",
                               "v_spoink_dsec" = "D. sechelia",
                               "v_spoink_dmau" = "D. mauritiana",
                               "v_spoink_dmel" = "D. melanogaster"))


plot(g_spoink_old)


g_spoink<-ggplot(data=df_spoink,aes(x=pos,y=species, color=species))+geom_point(alpha = 0.4)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  xlab("position")+ylab("SNPs")

plot(g_spoink)


gg_SNPs <- ggarrange(g_Shellder_old,g_spoink_old, ncol = 1, nrow = 2,align = ("v"),
                     labels = c("A", "B"), heights = c(2,2), widths = c(2,2))

plot(gg_SNPs)

#Proble the previous analysis include the ancient shared SNPs

Shellder_SNPs_exclude <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/shellder-ancient-coverage-0.txt", header = TRUE, sep = "\t")
spoink_SNPs_exclude <- read.csv("/Users/ascarpa/Downloads/Double_trouble_local/spoink-ancient-coverage-03.txt", header = TRUE, sep = "\t")

df_Shellder_clean <- df_Shellder %>%
  anti_join(Shellder_SNPs_exclude, by = c("pos" = "base"))

df_spoink_clean <- df_spoink %>%
  anti_join(spoink_SNPs_exclude, by = c("pos" = "base"))


g_Shellder_clean<-ggplot(data=df_Shellder_clean,aes(x=pos, fill=species))+geom_histogram(binwidth=1)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  xlab("position")+ylab("SNPs")+
  ggtitle("Shellder's SNPs in different species")+
  scale_fill_manual(values = c("v_Shellder_dsim" = "#009900", "v_Shellder_dmau" = "#0033cc", "v_Shellder_dsec" = "#FF9900", "v_Shellder_dtei" = "#ff6666"),
                    labels = c("v_Shellder_dsim" = "D. simulans",
                               "v_Shellder_dsec" = "D. sechelia",
                               "v_Shellder_dmau" = "D. mauritiana",
                               "v_Shellder_dtei" = "D. teissieri"))

plot(g_Shellder_clean)


g_spoink_clean<-ggplot(data=df_spoink_clean,aes(x=pos, fill=species))+geom_histogram(binwidth=1)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  xlab("position")+ylab("SNPs")+
  ggtitle("Spoink's SNPs in different species")+
  scale_fill_manual(values = c("v_spoink_dsim" = "#009900", "v_spoink_dmau" = "#0033cc", "v_spoink_dsec" = "#FF9900", "v_spoink_dmel" = "#CC0033"),
                    labels = c("v_spoink_dsim" = "D. simulans",
                               "v_spoink_dsec" = "D. sechelia",
                               "v_spoink_dmau" = "D. mauritiana",
                               "v_spoink_dmel" = "D. melanogaster"))


plot(g_spoink_clean)


gg_SNPs_clean <- ggarrange(g_Shellder_clean, df_spoink_clean, ncol = 1, nrow = 2,align = ("v"),
                     labels = c("A", "B"), heights = c(2,2), widths = c(2,2))

plot(gg_SNPs_clean)


#SNPs

df_Shell_sim_clean <- df_Shellder_clean %>%
  filter(species == "v_Shellder_dsim")
df_Shell_sec_clean <- df_Shellder_clean %>%
  filter(species == "v_Shellder_dsec")
df_Shell_mau_clean <- df_Shellder_clean %>%
  filter(species == "v_Shellder_dmau")
df_Shell_tei_clean <- df_Shellder_clean %>%
  filter(species == "v_Shellder_dtei")

m_sim_sec <- inner_join(df_Shell_sim_clean, df_Shell_sec_clean, by = "pos")
m_sim_mau <- inner_join(df_Shell_sim_clean, df_Shell_mau_clean, by = "pos")
m_sim_tei <- inner_join(df_Shell_sim_clean, df_Shell_tei_clean, by = "pos")

m_sec_mau <- inner_join(df_Shell_sec_clean, df_Shell_mau_clean, by = "pos")
m_sec_tei <- inner_join(df_Shell_sec_clean, df_Shell_tei_clean, by = "pos")

m_mau_tei <- inner_join(df_Shell_mau_clean, df_Shell_tei_clean, by = "pos")



df_overview_Shellder <- data.frame(
  D_sec = numeric(3),
  D_mau = numeric(3),
  D_tei = numeric(3)
)
rownames(df_overview_Shellder) <- c("D_sim", "D_sec", "D_mau")


df_overview_Shellder["D_sim", "D_sec"] <- nrow(m_sim_sec)
df_overview_Shellder["D_sim", "D_mau"] <- nrow(m_sim_mau)
df_overview_Shellder["D_sim", "D_tei"] <- nrow(m_sim_tei)

df_overview_Shellder["D_sec", "D_mau"] <- nrow(m_sec_mau)
df_overview_Shellder["D_sec", "D_tei"] <- nrow(m_sec_tei)

df_overview_Shellder["D_mau", "D_tei"] <- nrow(m_mau_tei)

print(df_overview_Shellder)

tot <- rbind(m_sim_sec, m_sim_mau, m_sim_tei, m_sec_mau, m_sec_tei, m_mau_tei)
write.csv(tot, "/Users/ascarpa/Downloads/Double_trouble_local/tot_SNPs_Shellder.csv", row.names=FALSE)



df_spoink_sim_clean <- df_spoink_clean %>%
  filter(species == "v_spoink_dsim")
df_spoink_sec_clean <- df_spoink_clean %>%
  filter(species == "v_spoink_dsec")
df_spoink_mau_clean <- df_spoink_clean %>%
  filter(species == "v_spoink_dmau")
df_spoink_mel_clean <- df_spoink_clean %>%
  filter(species == "v_spoink_dmel")

m_sim_sec_spoink <- inner_join(df_spoink_sim_clean, df_spoink_sec_clean, by = "pos")
m_sim_mau_spoink <- inner_join(df_spoink_sim_clean, df_spoink_mau_clean, by = "pos")
m_sim_mel_spoink <- inner_join(df_spoink_sim_clean, df_spoink_mel_clean, by = "pos")

m_sec_mau_spoink <- inner_join(df_spoink_sec_clean, df_spoink_mau_clean, by = "pos")
m_sec_mel_spoink <- inner_join(df_spoink_sec_clean, df_spoink_mel_clean, by = "pos")

m_mau_mel_spoink <- inner_join(df_spoink_mau_clean, df_spoink_mel_clean, by = "pos")


df_overview_spoink <- data.frame(
  D_sec = numeric(3),
  D_mau = numeric(3),
  D_mel = numeric(3)
)
rownames(df_overview_spoink) <- c("D_sim", "D_sec", "D_mau")


df_overview_spoink["D_sim", "D_sec"] <- nrow(m_sim_sec_spoink)
df_overview_spoink["D_sim", "D_mau"] <- nrow(m_sim_mau_spoink)
df_overview_spoink["D_sim", "D_mel"] <- nrow(m_sim_mel_spoink)

df_overview_spoink["D_sec", "D_mau"] <- nrow(m_sec_mau_spoink)
df_overview_spoink["D_sec", "D_mel"] <- nrow(m_sec_mel_spoink)

df_overview_spoink["D_mau", "D_mel"] <- nrow(m_mau_mel_spoink)

print(df_overview_spoink)

tot_spoink <- rbind(m_sim_sec_spoink, m_sim_mau_spoink, m_sim_mel_spoink, m_sec_mau_spoink, m_sec_mel_spoink, m_mau_mel_spoink)
write.csv(tot_spoink, "/Users/ascarpa/Downloads/Double_trouble_local/tot_SNPs_spoink.csv", row.names=FALSE)

