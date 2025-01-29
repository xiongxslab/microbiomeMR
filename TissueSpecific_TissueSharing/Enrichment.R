######## Enrichment of microbiome in tissue-specific and -sharing groups
shared_three_1 <- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/shared_three_new_1.txt", sep = "\t", header = FALSE)

t_and_s_1 <- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/t_and_s_new_1.txt", sep = "\t", header = FALSE)
s_and_i_1 <- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/s_and_i_new_1.txt", sep = "\t", header = FALSE)
i_and_t_1<- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/i_and_t_new_1.txt", sep = "\t", header = FALSE)
t_specific_1<- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/t_specific_new_1.txt", sep = "\t", header = FALSE)
s_specific_1<- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/s_specific_new_1.txt", sep = "\t", header = FALSE)
i_specific_1<- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/i_specific_new_1.txt", sep = "\t", header = FALSE)



combined_data <- rbind(shared_three_1,t_and_s_1,s_and_i_1,i_and_t_1,t_specific_1,s_specific_1,i_specific_1)  
colnames(combined_data) <- c("part1", "part2", "new_column")
combined_data <- combined_data%>% dplyr::select(microbiome = part1, gene= part2, character = new_column)



my_dataframe <- combined_data %>%
  mutate(microbiome = as.character(microbiome),
         character = as.character(character))

# Create the matrix using group_by and summarize
heatmap_matrix <- my_dataframe %>%
  group_by(microbiome, character) %>%
  summarize(n = n()) %>%  # Count the occurrences
  spread(character, n, fill = 0)  # Spread into matrix format
heatmap_matrix <- as.data.frame(heatmap_matrix)


# Set the rownames to be the microbiome types
rownames(heatmap_matrix) <- heatmap_matrix$microbiome
heatmap_matrix <- heatmap_matrix[, -1]  
desired_order <- c( "t_specific", "s_specific", "i_specific","t_and_s", "s_and_i", "i_and_t", "shared_three")

# Reorder the columns
heatmap_matrix <- heatmap_matrix %>%
  dplyr::select(all_of(desired_order))


odds_ratio_matrix <- matrix(0, nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix))
rownames(odds_ratio_matrix) <- rownames(heatmap_matrix)
colnames(odds_ratio_matrix) <- colnames(heatmap_matrix)

p_value_matrix <- matrix(0, nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix))
rownames(p_value_matrix) <- rownames(heatmap_matrix)
colnames(p_value_matrix) <- colnames(heatmap_matrix)

# Loop through each cell in the heatmap matrix
for (i in 1:nrow(heatmap_matrix)) {
  for (j in 1:ncol(heatmap_matrix)) {
    
    # The value of the current cell
    cell_value <- heatmap_matrix[i, j]
    
    # The total number in the column
    column_total <- sum(heatmap_matrix[, j])
    
    # The total number in the row
    row_total <- sum(heatmap_matrix[i, ])
    
    # The total number of all cells
    grand_total <- sum(heatmap_matrix)
    
    # Construct the contingency table
    contingency_table <- matrix(c(cell_value, 
                                  column_total - cell_value, 
                                  row_total - cell_value, 
                                  grand_total - row_total - (column_total - cell_value)),
                                nrow = 2, byrow = TRUE)
    
    # Perform Fisher's exact test
    fisher_test_result <- fisher.test(contingency_table)
    
    # Store the odds ratio in the corresponding cell in the odds_ratio_matrix
    odds_ratio_matrix[i, j] <- fisher_test_result$estimate
    
    p_value_matrix[i, j] <- fisher_test_result$p.value
  }
}
