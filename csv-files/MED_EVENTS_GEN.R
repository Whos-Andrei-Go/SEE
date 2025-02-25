# Load the necessary libraries
library(AdhereR)
library(this.path)

# Set working directory to the script's directory
setwd(this.path::this.dir())

# Define the file path (saved in the script's directory inside "csv-files")
current_dir <- getwd()

# Define the file name
file_name <- "med_events.csv"

# Create the full file path
file_path <- file.path(current_dir, file_name)

# Write the sorted data to a CSV file
write.csv(med.events, file = file_path, row.names = FALSE)

cat("CSV file 'med_events.csv' created successfully in", current_dir, "\n")
