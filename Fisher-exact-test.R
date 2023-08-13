## Perform Fisher's exact test for a specific treatment and write a csv file with the p-values ##

############# Test pour une ligne 
# Données pour la 3e ligne et la somme des colonnes
# row_data <- c(0, 0, 0, 1, 4, 1, 0, 0, 0)
# sum_row_data <- c(10, 22, 99, 941, 36560, 499, 1, 0, 0)

# Créer la matrice de contingence 2x9
# contingency_table <- rbind(row_data, sum_row_data)

# Effectuer le test exact de Fisher
# result <- fisher.test(contingency_table)

# Récupérer la p-value
# p_value <- result$p.value

# Afficher la p-value
# print(p_value)

############# Test avec le .csv de l'Atrazin 1
# File to perform Fisher's exact test
df <- read.csv2("test-melanie-atrazin-1.csv")

# Je passe col 1 en nom de ligne
rownames(df) <- df[, 1]
df <- df[, -1]

# On initialise un vecteur pour stocker les p-values
pv <- numeric(nrow(df))

# On effectue le test exact de Fisher pour chaque paire de lignes
for (i in 1:nrow(df)) {
  row_data <- df[i, ]
  sum_row_data <- c(10, 22, 99, 941, 36560, 499, 1, 0, 0)
  
  # Créer la matrice de contingence 2x9
  contingency_table <- rbind(row_data, sum_row_data)
  
  # Effectuer le test exact de Fisher
  result <- fisher.test(contingency_table)
  
  # Récupérer la p-value
  p_value <- result$p.value
  
  # Stocker la p-value dans le vecteur pv
  pv[i] <- p_value
}

# Write the new csv file with p-values
write.csv2(pv, "Diuron_1_simulated_p_values.csv")


############# Code Juliette Douillet
# File to perform Fisher's exact test
# df <- read.csv2("/Users/MelaniePietri/Desktop/test-melanie-atrazin-1.csv")

# Je passe col 1 en nom de ligne
# rownames(df) <- df[,1]
# df <- df[,-1]

# On applique la fonction fisher.test (apply) à toutes les lignes (1)
# (On a dû augmenter l'espace de travail à cause des dimensions du workspace)
# On récupère ensuite les p.value
# pv <- apply(df, 1, function(x) 
#  fisher.test(matrix(x, ncol=2), workspace=1e9)$p.value)

# Write the new csv file with p-values
# write.csv2(pv, "/Users/MelaniePietri/Desktop/test-melanie-atrazin-1-R.csv")


