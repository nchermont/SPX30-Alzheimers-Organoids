library(readxl)
library(dplyr)

# Caminho para o arquivo Excel
info_file <- read_excel("")


# Número de replicatas por condição
n_replicates <- 6

# Separar metadados (as colunas que nao sao abundancias)
meta_info <- info_file[, 1:2]

# Dados de abundância (trocar 4 pelo numero de colunas que nao sao abundancias)
abund_data <- info_file[, -(1:2)]

# Número total de colunas de abundância
n_cols <- ncol(abund_data)

# Criar grupos de colunas por condição
group_indices <- split(seq_len(n_cols), ceiling(seq_len(n_cols) / n_replicates))

# Função que verifica se cada linha tem pelo menos 2 quantificações em cada grupo
keep_protein <- apply(abund_data, 1, function(row) {
  all(sapply(group_indices, function(idxs) {
    sum(!is.na(row[idxs]) & row[idxs] != "0") > 3
  }))
})

# Filtrar o data frame
filtered_info <- cbind(meta_info[keep_protein, ], abund_data[keep_protein, ])

# Visualizar resultado
print(filtered_info)

library(openxlsx)

# Se quiser salvar:
write.xlsx(filtered_info, "/Users/namoreira/Desktop/MGL/data_NAfiltered3.xlsx")
