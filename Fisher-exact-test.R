## Perform Fisher's exact test for a specific treatment and write a csv file with the p-values ##

# File to perform Fisher's exact test
df <- read.csv2("data_example_simulated_R.csv")

# Name of the row in column 1
rownames(df) <- df[, 1]
df <- df[, -1]

# Initialize a vector to store the p-values
pv <- numeric(nrow(df))

# Fisher test for each line
for (i in 1:nrow(df)) {
  row_data <- df[i, ]
  sum_row_data <- c(84, 253, 50, 1152, 33463, 1340, 748, 358, 97)
  
  # Create cotingency table 2x9
  contingency_table <- rbind(row_data, sum_row_data)
  
  # Compute Fisher's exact test
  result <- fisher.test(contingency_table)
  
  # Store p-values
  p_value <- result$p.value
  
  # Store p-values in pv vector
  pv[i] <- p_value
}

# Write the new csv file with p-values
write.csv2(pv, "data_example_simulated_p_values.csv")
