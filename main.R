library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
    # Read in the counts file and metadata.
    meta <- read_csv(metafile_csv, show_col_types = FALSE)
  
    # Keep desired columns and subset to the desired timepoints.
    meta_sub <- meta |> 
      dplyr::select(samplename, timepoint) |> 
      dplyr::filter(timepoint %in% selected_times)
    
    # Set 'vP0' as the reference factor level for the timepoint column.
    meta_sub$timepoint <- factor(meta_sub$timepoint, levels = c(("vP0"), setdiff(selected_times, "vP0")))
  
    # Read count data.
    counts <- read_tsv(counts_csv, show_col_types = FALSE)
    counts_df <- as.data.frame(counts)
    rownames(counts_df) <- counts_df[[1]]                       # Make gene IDs (first column) the row names.
    counts_df <- counts_df[, -1, drop = FALSE]                  # Drop the gene IDs column afterwards.
    counts_df <- counts_df[, meta_sub$samplename, drop = FALSE] # Subset counts to selected samples.
    counts_mtx <- as.matrix(counts_df)
    
    # Build SummarizedExperiment object.
    se <- SummarizedExperiment(
      assays = list(counts = counts_mtx),    # Subsetted counts matrix
      colData = meta_sub                     # Matching sample metadata colData
    )
        
    return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  
    # Convert SummarizedExperiment to DESeqDataSet.
    dsq <- DESeqDataSet(se, design = design)
  
    # Run DESeq2.
    dsq <- DESeq(dsq)
    
    # Extract results and convert to dataframe.
    res <- results(dsq)
    res_df <- as.data.frame(res)
    
    # Return named list with both objects.
    return(list(
      dsq = dsq,
      results_df = res_df
    ))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
    # deseq2 results df = rows(genes), cols(stats(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj))
    
    # Convert DESeq2 results to a tibble.
    res_tibble <- deseq2_res |> 
      
      # Ensure col for rownames is called genes.
      rownames_to_column(var = "genes") |> 
      as_tibble() |> 
      
      # Add a column to denote plotting status (volc_plot_status).
      mutate(
        volc_plot_status = case_when(
          !is.na(padj) & padj < padj_threshold & log2FoldChange > 0 ~ "UP",
          !is.na(padj) & padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
          TRUE ~ "NS"
        )
      )

    return(res_tibble)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  
    plot <- labeled_results |> 
      ggplot(aes(x = pvalue)) +
      geom_histogram(bins = 50) +
      labs(
        title = "Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)",
        x = "pvalue",
        y = "count"
      ) +
      theme_minimal()
  
    return(plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
    
    plot <- labeled_results |> 
      filter(padj < padj_threshold) |> 
      ggplot(aes(x = log2FoldChange)) +
      geom_histogram(bins = 100) +
      labs(
        title = "Histogram of Log2FoldChanges for DE Genes (vP0 vs. vAd)",
        x = "log2FoldChange",
        y = "count"
      ) +
      theme_light()
  
    return(plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
    # Get top ten genes ranked by ascending padj.
    top_genes <- labeled_results |> 
      filter(!is.na(padj)) |> 
      arrange(padj) |> 
      slice_head(n = num_genes) |> 
      pull(genes)
    
    # Get normalized counts for each sample for each of the top genes.
    norm_counts <- counts(dds_obj, normalized = TRUE)
    top_counts <- norm_counts[top_genes, , drop = FALSE]
    
    # Convert to data frame and reshape to long format.
    plot_df <- as.data.frame(top_counts) |> 
      rownames_to_column(var = "genes") |> 
      pivot_longer(
        cols = -genes,
        names_to = "samplenames",
        values_to = "norm_counts"
      ) |> 
      mutate(log10_norm_counts = log10(norm_counts))
    
    # Generate the scatterplot.
    final_plot <- ggplot(plot_df, aes(x = genes, y = log10_norm_counts, color = samplenames)) +
      geom_point(position = position_dodge(width = 0.2)) +
      labs(
        title = paste("Plot of Log10(normalized counts) for top", num_genes, "DE genes"),
        x = NULL,
        y = "log10(norm_counts)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
    
    return(final_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
    
    plot <- labeled_results |> 
      filter(!is.na(padj), !is.na(log2FoldChange)) |> 
      ggplot(aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
      geom_point(size = 1.5) +
      # Mark the significance threshold (padj = 0.10) on the -log10 scale.
      geom_hline(yintercept = -log10(0.10), linetype = "dashed") +
      scale_color_manual(values = c("DOWN" = "salmon", "NS" = "green", "UP" = "cornflowerblue")) +
      labs(
        title = "Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)",
        x = "log2FoldChange",
        y = "-log10(padj)",
        color = "volc_plot_status"
      ) +
      theme_minimal()
    
    return(plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
    # Read in the mapping file.
    id2gene <- read_tsv(id2gene_path, col_names = c("ensembl_gene_id", "gene_symbol"))
  
    # Join with the DESeq2 results and remove missing values.
    merged <- labeled_results |> 
      inner_join(id2gene, by = c("genes" = "ensembl_gene_id")) |> 
      filter(!is.na(log2FoldChange), !is.na(gene_symbol))
    
    # Convert to named vector as specified, with log2FoldChange values in descending order.
    ranked_vec <- merged$log2FoldChange
    names(ranked_vec) <- merged$gene_symbol
    ranked_vec <- sort(ranked_vec, decreasing = TRUE)
    
    return(ranked_vec)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
    
    # Read the GMT file.
    pathways <- gmtPathways(gmt_file_path)
    
    # Remove duplicate gene symbols from the rank list.
    rnk_list <- rnk_list[!duplicated(names(rnk_list))]
  
    # Run FGSEA.
    fgsea_results <- fgsea(
      pathways = pathways,
      stats = rnk_list,
      minSize = min_size,
      maxSize = max_size
    )

    # Convert to a tibble.
    fgsea_results <- as_tibble(fgsea_results)
    
    return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
    # Clean the FGSEA results.
    fgsea_results <- fgsea_results |> 
      filter(!is.na(NES))
    
    # Get the top most positive NES pathways.
    top_pos <- fgsea_results |> 
      arrange(desc(NES)) |> 
      slice_head(n = num_paths)
    
    # Get the top most negative NES pathways.
    top_neg <- fgsea_results |> 
      arrange(NES) |> 
      slice_head(n = num_paths)
    
    # Combine the top positive and negative results into a dataframe.
    plot_df <- bind_rows(top_pos, top_neg) |> 
      arrange(NES) |> 
      mutate(
        pathway = gsub("_", " ", pathway),
        pathway = str_wrap(pathway, width = 40),
        pathway = factor(pathway, levels = pathway),
        direction = ifelse(NES > 0, "Positive", "Negative")
      )
    
    final_plot <- ggplot(plot_df, aes(x = NES, y = pathway, fill = direction)) +
      geom_col() +
      scale_fill_manual(
        values = c("Positive" = "red", "Negative" = "blue")
      ) +
      labs(
        title = "fgsea results for Hallmark MSigDB gene set",
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      ) +
      theme_minimal() + 
      theme(
        legend.position = "none",
        axis.text.y = element_text(size = 8)
      )
    
    return(final_plot)
}

