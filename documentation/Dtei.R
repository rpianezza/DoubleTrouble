library(tidyverse)
theme_set(theme_bw())

df_matute_supp <- read.csv("/Users/ascarpa/DoubleTrouble/data/dtei/Matute_clean.csv", header = TRUE, sep = ",")

df_TE <- read.csv("/Users/ascarpa/DoubleTrouble/data/dtei/Dtei_new_TEs.csv", header = TRUE, sep = ",")
df_TE <- df_TE%>%
  filter(Sample != "Sample")%>%
  filter(TE == "gypsy-29-dsim")

df_tsv <- read.csv("/Users/ascarpa/DoubleTrouble/data/dtei/filereport_read_run_PRJNA395473_tsv.txt", header = TRUE, sep = "\t")


result <- merge(df_TE, df_tsv, by.x = "Sample", by.y = "run_accession", all = FALSE)
result2 <- merge(result, df_matute_supp, by.x = "sample_alias", by.y = "Line", all = FALSE)

ggplot(result2, aes(x = Sample, y = as.numeric(All_reads), fill = Population)) +
  geom_bar(stat = "identity")

countries <- c("Bioko", "Equatorial Guinea", "Gabon", "Zimbabwe")
lat <- c(3.558445,1.547533, -0.989408, -19.150480)
long <- c(8.731123, 10.416840,11.561101,30.339835)
coordinates <- tibble(Population=countries, lat=lat, long=long)


result3 <-inner_join(result2, coordinates, by = "Population")


plot_map <- function(dataset) {
  dataset$TE <- factor(dataset$TE, levels = c("gypsy-29-dsim"))
  world_map <- map_data("world")
  africa_map <- subset(world_map, region %in% c(
    "Algeria", "Angola", "Benin", "Botswana", "Burundi", "Cameroon", "Cape Verde", 
    "Central African Republic", "Chad", "Comoros", "Congo", "Democratic Republic of the Congo",
    "Djibouti", "Egypt", "Equatorial Guinea", "Eritrea", "Swaziland", "Ethiopia", "Gabon",
    "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Ivory Coast", "Kenya", "Lesotho", 
    "Liberia", "Libya", "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius", "Morocco", 
    "Mozambique", "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal", 
    "Seychelles", "Sierra Leone", "Somalia", "South Africa", "South Sudan", "Sudan", 
    "Tanzania", "Togo", "Tunisia", "Uganda", "Zambia", "Zimbabwe", "Western Sahara", 
    "Mauritania", "Mali", "Niger", "Chad", "Sudan", "Libya", "Egypt", 
    "Burkina Faso", "Lesotho", "Republic of Congo"
  ))
  
  ggplot() +
    geom_map(data = africa_map, map = africa_map,
             aes(long, lat, map_id = region),
             color = "white", fill = "burlywood3", size = 0) +
    geom_point(data = dataset, aes(x = long, y = lat, color = as.numeric(HQ_reads)), size = 2, position = position_jitter(width = 0.1, height = 0.1), alpha = 0.7) +
    scale_colour_gradient(low = "darkgreen", high = "red") + 
    theme(plot.title = element_text(hjust = 0.5), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none")
}

plot_map(result3)



pca_dtei_eigenval <- read.csv("/Users/ascarpa/DoubleTrouble/data/dtei/dtei.pca.eigenval", header = FALSE)
names(pca_dtei_eigenval) <- c("col_val")

pca_dtei_eigenvec <- read.csv("/Users/ascarpa/DoubleTrouble/data/dtei/dtei.pca.eigenvec", sep = "", header = FALSE)
PCs <-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13")
names(pca_dtei_eigenvec) <- c("ID","ID2","PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13")

pca_dtei_eigenvec <- pca_dtei_eigenvec %>%
  mutate(ID = str_extract(ID, "(?<=/)[^/]+(?=.fastq)"))

miniinfo <- subset(result3, select = c("Sample", "Population", "HQ_reads"))
pca_dtei_eigenvec <- merge(pca_dtei_eigenvec, miniinfo, by.x = "ID", by.y = "Sample", all = FALSE)

pca_dtei_eigenvec <- pca_dtei_eigenvec %>%
  mutate(Presence = ifelse(HQ_reads > 1, "present", "absent"))


eigen2pca <- function(vec, val, title) {
  
  eigenval <-  val %>% mutate(PC = PCs, variability_explained = paste0(round((col_val/sum(col_val)*100),2), "%"))
  
  var_explained <- eigenval %>% select(variability_explained) %>% pull()
  
  pca_plot_islands <- ggplot(vec, aes(x=PC1, y=PC3, color=Population)) + geom_point(alpha=0.5, size=4) +
    xlab(paste0("PC1: ", var_explained[1])) + ylab(paste0("PC3: ", var_explained[3])) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  
  pca_plot_presence <- ggplot(vec, aes(x=PC1, y=PC3, color=Presence)) + geom_point(alpha=0.5, size=4) +
    xlab(paste0("PC1: ", var_explained[1])) + ylab(paste0("PC3: ", var_explained[3])) +
    scale_color_manual(values = c("darkgreen", "red")) + labs(color = "Shellder") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  
  pca_plot_merged <- ggplot(vec, aes(x=PC1, y=PC3, color=Presence, fill=Population)) + geom_point(shape=21, alpha=0.5, size=4, stroke=2) +
    xlab(paste0("PC1: ", var_explained[1])) + ylab(paste0("PC3: ", var_explained[3])) +
    scale_color_manual(values = c("darkgreen", "red", "orange")) + labs(color = "TE carried") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 2), legend.position = "bottom")
  
  list(islands = pca_plot_islands, presence = pca_plot_presence, merged = pca_plot_merged)
  
}

eigen2pca(pca_dtei_eigenvec, pca_dtei_eigenval, "PCA 7412533 SNPs Drosophila teissieri")

