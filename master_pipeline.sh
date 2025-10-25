#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
output_heatmap <- args[2]
output_table <- args[3]
gtf_file <- args[4]  # Add GTF file as parameter

# Create directories for outputs
dir.create(dirname(output_heatmap), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)

# Function to safely load packages
safe_library <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat(paste("Could not install package:", package_name, "\n"))
      return(FALSE)
    }
  }
  return(TRUE)
}

# Try to load required libraries
dplyr_loaded <- safe_library("dplyr")
pheatmap_loaded <- safe_library("pheatmap")

# Exit if libraries aren't available
if (!dplyr_loaded || !pheatmap_loaded) {
  cat("Required R packages could not be loaded. Exiting.\n")
  quit(status = 1)
}

# Function to extract gene names from GTF
extract_gene_names <- function(gtf_file) {
  if (!file.exists(gtf_file)) {
    cat("GTF file not found. Using gene IDs as names.\n")
    return(NULL)
  }
  
gtf_lines <- readLines(gtf_file)
gene_map <- data.frame(gene_id = character(), gene_name = character(), stringsAsFactors = FALSE)
  
  for (line in gtf_lines) {
    if (grepl("gene_id", line)) {
      gene_id_match <- regmatches(line, regexpr('gene_id "[^"]*"', line))
      gene_name_match <- regmatches(line, regexpr('gene_name "[^"]*"', line))
      
      if (length(gene_id_match) > 0) {
        gene_id <- gsub('gene_id "|"', '', gene_id_match)
        gene_name <- if (length(gene_name_match) > 0) {
          gsub('gene_name "|"', '', gene_name_match)
        } else {
          gene_id
        }
        
        if (!gene_id %in% gene_map$gene_id) {
          gene_map <- rbind(gene_map, data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE))
        }
      }
    }
  }
  
  return(gene_map)
}

# Try to load raw counts
counts_loaded <- tryCatch({
  counts <- read.delim(counts_file, comment.char = "#")
  TRUE
}, error = function(e) {
  cat(paste("Error loading counts file:", e$message, "\n"))
  FALSE
})

if (!counts_loaded) {
  cat("Could not load counts file. Exiting.\n")
  quit(status = 1)
}

# Safety check for empty or invalid data
if (nrow(counts) == 0 || ncol(counts) <= 1) {
  cat("Counts file appears to be empty or invalid. Exiting.\n")
  quit(status = 1)
}

tryCatch({
  # Load gene name mapping from GTF
  gene_map <- extract_gene_names(gtf_file)
  
  # Get count columns (all numeric columns except the first)
  count_cols <- which(sapply(counts, is.numeric))
  if (length(count_cols) == 0) {
    cat("Could not identify count columns. Using all columns except first.\n")
    count_cols <- 2:ncol(counts)
  }
  
  # Exclude first column (Geneid) from count_cols if it's included
  if (1 %in% count_cols) {
    count_cols <- count_cols[count_cols != 1]
  }
  
  count_data <- counts[, count_cols, drop = FALSE]
  
  # Extract sample names from BAM file paths (SRR9613403, SRR9613404, SRR9613405)
  samples <- gsub(".*\/(SRR[0-9]+)\.bam", "\1", colnames(count_data))
  colnames(count_data) <- samples
  
  # Keep gene IDs separately
  gene_ids <- counts$Geneid
  
  # Total counts per sample
  total_counts <- colSums(count_data)
  
  # Safety check for zero total counts
  if (any(total_counts == 0)) {
    cat("Warning: Some samples have zero total counts. Using pseudocount for normalization.\n")
    total_counts[total_counts == 0] <- 1
  }
  
  # CPM normalization
  cpm_data <- sweep(count_data, 2, total_counts, FUN=function(x, y) (x / y) * 1e6)
  cpm_data <- data.frame(Geneid = gene_ids, cpm_data, stringsAsFactors = FALSE)
  
  # Map gene IDs to gene names
  if (!is.null(gene_map)) {
    cpm_data <- merge(gene_map, cpm_data, by.x = "gene_id", by.y = "Geneid", all.y = TRUE)
    cpm_data$gene_name[is.na(cpm_data$gene_name)] <- cpm_data$gene_id[is.na(cpm_data$gene_name)]
  } else {
    cpm_data$gene_name <- cpm_data$Geneid
    cpm_data$gene_id <- cpm_data$Geneid
  }
  
  # Identify top 10 genes by mean CPM
  top_genes <- cpm_data %>%
    mutate(mean_cpm = rowMeans(select(., starts_with("SRR")))) %>%
    arrange(desc(mean_cpm)) %>%
    head(10)
  
  # Create table with gene names for output
  top_genes_output <- top_genes %>%
    select(gene_id, gene_name, everything())
  
  # Extract just the numerical columns for heatmap with gene names as row names
  heatmap_data <- top_genes %>% 
    select(starts_with("SRR"))
  rownames(heatmap_data) <- top_genes$gene_name
  
  # Save top genes table
  write.csv(top_genes_output, file = output_table, row.names = FALSE)
  
  # Create heatmap with proper labels
  pdf(output_heatmap, width = 10, height = 8)
  pheatmap(log2(heatmap_data + 1),
           main = "Top 10 Expressed Genes (CPM, log2 scale)",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 10,
           fontsize_col = 12,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           border_color = "grey60")
  dev.off()
  
  cat("Analysis complete!\n")
  cat("Top 10 genes table saved to:", output_table, "\n")
  cat("Heatmap saved to:", output_heatmap, "\n")
  cat("\nTop 10 genes:\n")
  print(top_genes_output %>% select(gene_id, gene_name, mean_cpm))
}, error = function(e) {
  cat(paste("Error in R analysis:", e$message, "\n"))
  quit(status = 1)
}
