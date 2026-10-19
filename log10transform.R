# ---------------------------------------------
# 🔹 LOG10 TRANSFORMATION FOR PROTEIN ABUNDANCE DATA
# ---------------------------------------------

library(dplyr)

# STEP 1: Choose your input CSV
file_path <- file.choose()
data <- read.csv(file_path, stringsAsFactors = FALSE)

# STEP 1.5 (optional but recommended):
# Convert numeric-looking columns to numeric automatically
data <- data %>%
  mutate(across(everything(), ~ type.convert(.x, as.is = TRUE)))

# STEP 2: Apply log10(+1) to ALL numeric columns except the first
data_log10 <- data %>%
  mutate(across(
    .cols = 2:ncol(.),       # apply to columns 2 through last
    .fns = ~ if (is.numeric(.x)) log10(.x + 1) else .x
  ))

# STEP 3: Save output in the same folder with "_log10" added
save_path <- sub("\\.csv$", "_log10.csv", file_path)
write.csv(data_log10, save_path, row.names = FALSE)

cat("✅ Log10 transformation complete!\n",
    "File saved to:", save_path, "\n")
