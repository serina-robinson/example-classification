# Load package
pacman::p_load("tidyverse")

# Read in the dataset and sort by fuzzy matching into four groups
# Can also try 'agrepl' for fuzzy pattern matching, but I've had limited success with it
dat <- read_csv("data/toy_dataset.csv")  %>%
  dplyr::mutate(description = tolower(description)) %>% # makes everything lower case
  dplyr::mutate(label = case_when(grepl(paste0(c("vehicle", "car", "truck"), collapse = "|"), description) ~ "vehicle",
                                  grepl(paste0(c("person", "man"), collapse = "|"), description) ~ "human",
                                  grepl(paste0(c("rubber", "ducky"), collapse = "|"), description) ~ "rubber duck",
                                  is.na(description) ~ NA_character_,
                                  TRUE ~ "other"))
