rm(list = ls()) #clear data

#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "vegan", "ggpubr","janitor", "ape")

#import metadata
metadata2 <- readr::read_csv("data/sample_metadata.csv") %>%
  dplyr::mutate(short_identifier = tolower(stringr::str_c(location, site, paste0("r",replicate), sep = "_"))) 
metadata2$short_identifier[metadata2$short_identifier == "alt_up2_r3"] <- "alt_up2_r2" # naming inconsistency fixed
metadata2$short_identifier[metadata2$short_identifier == "alt_dw1_r3"] <- "alt_dw1_r2" # naming inconsistency fixed
metadata3 <- readr::read_csv("data/94_metadata_with_biosamples_final_samples_annotated.csv")  %>%
  dplyr::mutate(short_identifier = tolower(paste0(word(sample_name, 2, sep = "_"), "_", word(sample_name, 3, sep = "_"), "_", word(sample_name, 4, sep = "_")))) %>%


# Biofilm data
bf <- read_csv("data/BF_data_heatmap.csv") 
bf$short_identifier <- bf$loc_batch
bfwithmet <- bf %>%
  left_join(., metadata3, by = "short_identifier") %>%
  dplyr::select(short_identifier, substance, k) %>%
  dplyr::filter(grepl("Sulf", substance)) %>%
  dplyr::filter(!grepl("thia", substance))
bfwithmet
bfwide <- reshape2::dcast(bfwithmet, short_identifier ~ substance, fill = "k") %>%
  column_to_rownames(var = "short_identifier") 
bfnume <- as.data.frame(sapply(bfwide, as.numeric))
rownames(bfnume) <- rownames(bfwide)

# KO data
ko_data <- readr::read_csv("data/20250128_biofilm_samples_sulfonamide_ko.csv")
colnames(ko_data)[1] <- "accession"
metadata <- metadata2 %>%
  left_join(metadata3, by = c("short_identifier")) 
kom <- ko_data %>%
  left_join(., metadata3, by = "accession")
kowithmet <- kom %>%
  dplyr::select(short_identifier, contains("ko")) %>%
  column_to_rownames(var = "short_identifier")

# Check overlap
metadata2$short_identifier %in% bf$loc_batch
unique(bf$loc_batch)[!unique(bf$loc_batch) %in% metadata2$short_identifier]


#square root transformation
ko_data <- kowithmet %>%
  dplyr::mutate(across(.cols = everything(), sqrt))
bfwide
bf_data <- bfnume %>%
  dplyr::mutate(across(.cols = everything(), sqrt))

#distance matrix
dist_method <- "bray"
ko_dist_mat <- ko_data %>%
  vegan::vegdist(method = dist_method)
bf_dist_mat <- bf_data %>%
  vegan::vegdist(method = dist_method)
#pcoa
ko_pcoa <- ape::pcoa(ko_dist_mat)
bf_pcoa <- ape::pcoa(bf_dist_mat)
#procrustes
crustypro <- vegan::procrustes(X = ko_pcoa$vectors, Y = bf_pcoa$vectors, symmetric = TRUE)
ctest <- tibble(sample_code = rownames(crustypro$X),
                yrda1 = crustypro$Yrot[,1],
                yrda2 = crustypro$Yrot[,2],
                xrda1 = crustypro$X[,1],
                xrda2 = crustypro$X[,2]) %>%
  dplyr::left_join(dplyr::select(metadata, sample_code, site, location), by = "sample_code")
plot <- ggplot(ctest) +
  geom_point(aes(x=xrda1, y=xrda2), size = 5) +
  geom_segment(aes(x=xrda1,y=xrda2,xend=yrda1,yend=yrda2, color = site), arrow=arrow(length=unit(0.5,"cm")), linewidth = 1.5) +
  theme_pubr() +
  ggplot2::theme(legend.position = "right",
                 legend.text = element_text(size = 28),
                 legend.title = element_text(size = 32),
                 strip.text.x = element_text(size = 32, face = "italic"),
                 axis.title = element_text(size = 32),
                 axis.text = element_text(size = 28)) +
  ggplot2::scale_color_manual(values = c("brown4", "orange", "purple", "green", "blue1", "blue3", "blue4", "gray", "cyan", "darkgreen", "red1", "red3")) +
  ggplot2::scale_shape_manual(values = c(15, 16, 17, 0, 1, 2)) +
  ggplot2::labs(x = "Dimension 1", y = "Dimension 2")
plot
#save plot
ggsave("plots/procrustes.jpg",
       dpi = 300,
       device = "jpeg",
       units = "cm",
       width = 75/2,
       height = (75/2)/1.5)
#statistical test
stat_test <- protest(X = bf_pcoa$vectors, Y = ko_pcoa$vectors, scores = "sites", permutations = 999)
stat_test
