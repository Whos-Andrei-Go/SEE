# Load the AdhereR package
library(AdhereR)
library(this.path)
setwd(this.path::this.dir())  # Sets to the script's directory

set.seed(1234)  # For reproducibility

# Generate 1080 rows of sample data
med.events <- data.frame(
  PATIENT_ID = sample(1:100, 1080, replace = TRUE),  # 100 unique patients
  DATE = sample(seq(as.Date("2020-01-01"), as.Date("2023-01-01"), by = "day"), 1080, replace = TRUE),
  PERDAY = sample(1:10, 1080, replace = TRUE),
  CATEGORY = sample(c("medA", "medB", "medC"), 1080, replace = TRUE),
  DURATION = sample(5:50, 1080, replace = TRUE)
)

# Get the current working directory (where your R script is located)
current_dir <- file.path(getwd(), "csv-files")

# Define the file name
file_name <- "med_events.csv"

# Create the full file path
file_path <- file.path(current_dir, file_name)

# Write the data to a CSV file in the current working directory
write.csv(med.events, file = file_path, row.names = FALSE)

cat("CSV file 'med_events.csv' created successfully in", current_dir, "\n")
