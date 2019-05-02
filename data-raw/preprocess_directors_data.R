## Code used to create the 'directors' data
## Data is proprietary and cannot be shared

rm(list=ls())
# Subset
library(dplyr)
directors <- read.csv2("testing/FinalData.csv") %>%
  # Select only independent directors
  filter(isED == 0) %>%
  # Select only these variables
  select(Recent_Comp_Thous, Sector, isM, Age, Company, Industry) %>%
  # Remove NA
  filter(!is.na(Recent_Comp_Thous),
         !is.na(Age)) %>%
  # Filter small sectors
  filter(!Sector %in% c("Money Center Banks", "Property & Casualty Insurance")) %>%
  # Filter sectors and add re-factor
  #filter(Sector %in% c("Financial", "Services", "Basic Materials")) %>%
  mutate(Sector = factor(as.character(Sector)),
         Industry = factor(as.character(Industry)),
                         #levels = c("Financial", "Services", "Basic Materials")),
         Male = factor(isM, levels=c(0,1)),
         Company = as.character(Company)) %>%
  select(-isM) %>%
  rename(Compensation = Recent_Comp_Thous)

# Save data
devtools::use_data(directors, overwrite = TRUE)
