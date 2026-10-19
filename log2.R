# ---------------------------------------------
# 🔹 LOG2 TRANSFORMATION FOR PROTEIN ABUNDANCE DATA
# ---------------------------------------------

library(dplyr)

# STEP 1: Choose your input CSV
file_path <- file.choose()
data <- read.csv(file_path, stringsAsFactors = FALSE)

# STEP 1.5 (optional but recommended):
# Convert numeric-looking columns to numeric automatically
data <- data %>%
  mutate(across(everything(), ~ type.convert(.x, as.is = TRUE)))

# STEP 2: Apply log2(+1) to selected numeric columns
data_log2 <- data %>%
  mutate(across(
    .cols = 2:ncol(.),       
    .fns = ~ if (is.numeric(.x)) log2(.x + 1) else .x
  ))

# STEP 3: Save output in the same folder with "_log2" added
save_path <- sub("\\.csv$", "_log2.csv", file_path)
write.csv(data_log2, save_path, row.names = FALSE)

cat("✅ Log2 transformation complete!\n",
    "File saved to:", save_path, "\n")