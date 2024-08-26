select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
reduce <- purrr::reduce
##Create a ensemble ID to gene symbol table
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

options(future.globals.maxSize = 1e9)

make_gene_id_table <- function(all_ensembl_ids) {
  # Fetch the gene annotation for all unique Ensembl IDs
  gene_id_table <- getBM(filters = "ensembl_gene_id",
                         attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                         values = all_ensembl_ids,
                         mart = mart) %>%
    rename(ENSEMBL_ID = ensembl_gene_id,
           HGNC_SYMBOL = hgnc_symbol,
           GENE_BIOTYPE = gene_biotype) %>%
    distinct()  # Ensure uniqueness of rows
  
  return(gene_id_table)
}


#Standardize rownames (because some were converted to symbols by Azimuth)
standardize_row_names <- function(seu_obj, gene_id_table) {
  gene_map <- setNames(gene_id_table$HGNC_SYMBOL, gene_id_table$ENSEMBL_ID)
  current_rownames <- rownames(seu_obj)
  standardized_rownames <- ifelse(current_rownames %in% names(gene_map),
                                  gene_map[current_rownames],
                                  current_rownames)
  rownames(seu_obj) <- standardized_rownames
  return(seu_obj)
}
# Step 1-4: Function to process markers for each sample
filter_markers <- function(markers_df, sample_name, fc_threshold = 1, pct_diff = 0.07, p_val_threshold = 0.05) {
  markers_df %>%
    filter(avg_log2FC > fc_threshold, 
           pct.1 - pct.2 > pct_diff, 
           p_val_adj < p_val_threshold) %>%
    arrange(desc(avg_log2FC)) %>%
    mutate(sample = sample_name)
}

process_seurat <- function(seu_obj) {
  tryCatch({
    # Basic cell-level quality filtering
    seu_obj <- seu_obj[rowSums(seu_obj@assays$RNA@counts) > 0, ]  # drop truly all-zero genes
    seu_obj <- seu_obj[, seu_obj$nCount_RNA > 200]  # at least 200 counts detected in every cell 
    
    # Perform SCTransform, PCA, UMAP, and clustering
    seu_obj <- SCTransform(seu_obj, vars.to.regress = "spliced_percent_MT", verbose = FALSE)
    seu_obj <- RunPCA(seu_obj, verbose = FALSE)
    seu_obj <- RunUMAP(seu_obj, dims = 1:30, verbose = FALSE)
    seu_obj <- RunTSNE(seu_obj)
    seu_obj <- FindNeighbors(seu_obj, dims = 1:30, verbose = FALSE)
    seu_obj <- FindClusters(seu_obj, resolution = 1, verbose = FALSE)
    
    return(seu_obj)
  }, error = function(e) {
    message("Processing failed: ", conditionMessage(e))
    return(NULL)
  })
  
  find_cluster_markers <- function(seu_obj) {
  tryCatch({
    # Find cluster markers
    obj_markers <- FindAllMarkers(seu_obj, only.pos = TRUE)
    obj_markers <- obj_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
    
    return(obj_markers)
  }, error = function(e) {
    message("Marker identification failed: ", conditionMessage(e))
    return(NULL)
  })
}

}

