##

library(tidyverse)
library(dplyr)
library(here)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(circlize)
library(viridis)
library(gridExtra)
library(magrittr)
library(Cairo)
library(stringr)
library(VennDiagram)
library(grid)
library(plotly)
library(cowplot)
library(purrr)
library(fs)
# library(ComplexHeatmap)
library(nsprcomp)
library(ggrepel)

##

panscreen_long <- readRDS("/Users/triandafillou/RajLab Dropbox/Catherine Triandafillou/Shared_CatT/misc-data/2024-10_figure-workshop/panscreen_long.rds")
lab_validation_data <- readRDS("/Users/triandafillou/RajLab Dropbox/Catherine Triandafillou/Shared_CatT/misc-data/2024-10_figure-workshop/validationdata_20240403.rds")

# file_path <- dir_ls(here(), recurse = TRUE) %>%
#   keep(~ str_detect(., fixed("validationdata_20240403.rds")))

# Check if at least one file was found and read it
# if (length(file_path) > 0) {
#   lab_validation_data <- readRDS(file_path[1]) # Reading the first match
#   head(lab_validation_data)
# } else {
#   cat("File not found.")
# }

#Rename the Drug column in the lab validatation dataset to match the taxonomy used for everything else
lab_validation_data = lab_validation_data %>% rename(Name = Drug)

##

custom_colors <- c("Parental" = "#5E5E5E", "D13" = "#FF7E79", "E2" = "#4294F8", "C12" = "#6539F8", "A15" = "#3E8D27", "E2 + D13" = "#8E1356", "C14" = "#953877", "D5" = "#6ADCC1", "D4" = "#9E6D15",
                   "A2" = "#EF8BF9", "A11" = "#CC5C76", "C3" = "#F57946", "E8" = "#F9AD2A", "F9" = "#1D457F", "I11" = "#625A94", "G12" = "#0013F0")

custom_shapes <- c("DMSO" = 16, "Vemurafenib" = 17, "DabTram" = 15)

custom_theme <- function() {
  theme(
    text = element_text(size = 14, family = "Arial", color = "black"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "white", linewidth = 0.8),
    panel.grid.minor = element_line(color = "white", linewidth = 0.4),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 12, face = "bold", hjust = 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
}

## Task 7 ----

#create a data frame of the HTS data with only the drugs that have been used for validation
hts_validation_data = panscreen_long %>% filter(Name %in% c("(+)-JQ1", "Dasatinib", "Clofarabine", "Saracatinib (AZD0530)", "OTX015", "Epirubicin HCl", "S3I-201")) %>% select(-Target) %>% select(-Pathway)

#add a viability column to the HTS data by mutating the Toxicity column and then removing it
hts_validation_data = hts_validation_data %>% mutate(Viability_HTS = 100 - Toxicity) %>% select(-Toxicity)

#Mutate the Viability column in the lab data so it's in percent format rather than ratio format to match the HTS data
validation_data = lab_validation_data %>% mutate(Viability_LAB = 100 * Viability) %>% select(-Viability)

#set up a dataframe with all of the naming inconsistencies between the HTS data and the lab data
name_mapping <- data.frame(NonStandard = c("JQ1", "Saracatinib"), Standard = c("(+)-JQ1", "Saracatinib (AZD0530)"))

#loop through the mapping data frame and update names in the lab data frame

for (i in 1:nrow(name_mapping)) {
  validation_data$Name[validation_data$Name == name_mapping$NonStandard[i]] <- name_mapping$Standard[i]
}

# Note this can be acomplished with an ifelse so can be done in the tidyverse
validation_data2 <- validation_data %>%
  mutate(Name2 = case_when(Name == "JQ1" ~ "(+)-JQ1",
                           Name == "Saracatinib" ~ "Saracatinib (AZD0530)",
                           TRUE ~ Name))

#filter only the cell lines that have been used in HTS from the lab validation data
validation_data = validation_data %>% filter(CellLine %in% c("Parental", "D13", "E2", "C12", "D4"))

#Calculate the average for each cell line/dose/drug/treatment combination for each data frame
avg_hts_validation_data = hts_validation_data %>% group_by(CellLine, Name, Dose, Treatment) %>% summarize(AVG_Viability_HTS = mean(Viability_HTS, na.rm = TRUE), .groups = 'drop')

avg_lab_validation_data = validation_data %>% group_by(CellLine, Name, Dose, Treatment) %>% summarize(AVG_Viability_LAB = mean(Viability_LAB, na.rm = TRUE), .groups = 'drop')

#Use full_join to combine the two data frames while retaining all rows
all_validation_data = full_join(avg_lab_validation_data, avg_hts_validation_data, by = c("Name", "Dose", "Treatment", "CellLine"))

ggplot(all_validation_data, aes(x = AVG_Viability_HTS, y = AVG_Viability_LAB, color = CellLine)) +
  geom_point(aes(shape = Treatment), position = position_dodge(width = 0.1), size = 2) +
  scale_color_manual(values = custom_colors) +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Adds the y=x line 
  labs(title = "Comparison of Toxicity Between HTS and In-Lab Screens",
       x = "High-Throughput Screen % Viability",
       y = "In-Lab Screen % Viability",
       color = "Cell Line") +
  custom_theme() +
  coord_fixed()

ggplot(all_validation_data, aes(x = AVG_Viability_HTS, y = AVG_Viability_LAB, color = CellLine)) +
  geom_point(aes(shape = Treatment), position = position_dodge(width = 0.1), size = 2) +
  scale_color_manual(values = custom_colors) +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Adds the y=x line 
  facet_grid( ~Dose) +
  labs(title = "Comparison of Toxicity Between HTS and In-Lab Screens",
       x = "High-Throughput Screen % Viability",
       y = "In-Lab Screen % Viability",
       color = "Cell Line") +
  custom_theme() +
  coord_fixed()  
