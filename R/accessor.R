#' Get cluster differential expression results
#'
#' @param circObj circObj with completed cluster DE analysis
#' @param cluster_type Type of cluster to retrieve ("A5BS", "A3BS", or "both")
#'
#' @return Data frame with cluster DE results
#' @export
getClusterDE <- function(circObj, cluster_type = "both") {

  # Check if cluster DE results exist
  if (!("circCluster.DE" %in% slotNames(circObj)) || is.null(circObj@circCluster.DE)) {
    stop("No cluster DE results found. Please run circClusterDE() first.")
  }

  results <- list()

  # Get A5BS results
  if (cluster_type %in% c("both", "A5BS") && "A5BS.cluster" %in% names(circObj@circCluster.DE)) {
    a5bs_results <- circObj@circCluster.DE$A5BS.cluster
    a5bs_results$cluster_type <- "A5BS"
    results$A5BS <- a5bs_results
  }

  # Get A3BS results
  if (cluster_type %in% c("both", "A3BS") && "A3BS.cluster" %in% names(circObj@circCluster.DE)) {
    a3bs_results <- circObj@circCluster.DE$A3BS.cluster
    a3bs_results$cluster_type <- "A3BS"
    results$A3BS <- a3bs_results
  }

  if (length(results) == 0) {
    stop(paste("No", cluster_type, "cluster DE results found."))
  }

  # Combine results
  if (length(results) == 1) {
    final_results <- results[[1]]
  } else {
    final_results <- do.call(rbind, results)
    rownames(final_results) <- NULL
  }

  # Reorder columns with cluster_type as last column
  col_order <- c("juncid", "numcircs", "m0", "m1", "log2fc", "pvalue", "fdr")
  existing_cols <- intersect(col_order, colnames(final_results))
  other_cols <- setdiff(colnames(final_results), c(existing_cols, "cluster_type"))
  final_results <- final_results[, c(existing_cols, other_cols, "cluster_type")]

  return(final_results)
}


#' Get cluster composition and annotation details
#'
#' @param circObj circObj with cluster information
#' @param cluster_type Type of cluster ("A5BS" or "A3BS")
#' @param juncid Vector of junction IDs to get annotations for
#'
#' @return Data frame with circRNA-level annotation for specified clusters
#' @export
#' @importFrom dplyr group_by summarise first
#' @importFrom magrittr %>%
getClusterAnnotation <- function(circObj, cluster_type, juncid) {

  # Validate cluster_type
  if (!cluster_type %in% c("A5BS", "A3BS")) {
    stop("cluster_type must be either 'A5BS' or 'A3BS'")
  }

  # Get the appropriate cluster slot
  slot_name <- paste0(cluster_type, ".cluster")
  if (!(slot_name %in% slotNames(circObj))) {
    stop(paste("Slot", slot_name, "not found. Please run getCircCluster() first."))
  }

  cluster_data <- slot(circObj, slot_name)

  # Filter by specified juncid(s)
  cluster_subset <- cluster_data[cluster_data$juncid %in% juncid, ]

  if (nrow(cluster_subset) == 0) {
    stop(paste("No data found for juncid(s):", paste(juncid, collapse = ", ")))
  }

  # Extract chromosome from 'from' field
  cluster_subset$chr <- sub(":.*", "", cluster_subset$from)

  # Create circRNA-level summary (deduplicate by circid)
  annotation_summary <- cluster_subset %>%
    group_by(circid, juncid) %>%
    summarise(
      chr = first(chr),
      start = first(from.pos),
      end = first(end.pos),
      host_gene = first(genename),
      .groups = 'drop'
    )

  # Add cluster information
  annotation_summary$cluster_type <- cluster_type

  # Calculate circRNA width
  annotation_summary$width <- annotation_summary$end - annotation_summary$start + 1

  # Reorder columns with cluster_type at the end
  col_order <- c("juncid", "circid", "chr", "start", "end", "width", "host_gene", "cluster_type")
  annotation_summary <- annotation_summary[, col_order]

  # Sort by juncid then circid
  annotation_summary <- annotation_summary[order(annotation_summary$juncid, annotation_summary$circid), ]

  return(as.data.frame(annotation_summary))
}
