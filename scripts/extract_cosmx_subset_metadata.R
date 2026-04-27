#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(jsonlite)))
suppressWarnings(suppressPackageStartupMessages(library(Seurat)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: extract_cosmx_subset_metadata.R <seurat_rds> <fov_spec>")
}

rds_path <- args[[1]]
fov_spec <- args[[2]]

trim_or_empty <- function(x) {
  if (length(x) == 0 || is.null(x)) {
    return("")
  }
  text <- trimws(as.character(x[[1]]))
  if (is.na(text)) {
    return("")
  }
  text
}

unique_nonempty <- function(x) {
  values <- trimws(as.character(x))
  unique(values[!is.na(values) & nzchar(values)])
}

fmt_num <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) {
    return("")
  }
  format(x, scientific = FALSE, trim = TRUE, digits = 15)
}

parse_fovs <- function(text) {
  text <- trimws(text)
  if (!nzchar(text)) {
    return(integer())
  }

  pieces <- unlist(strsplit(text, ",", fixed = TRUE))
  values <- integer()
  for (piece in pieces) {
    piece <- trimws(piece)
    if (!nzchar(piece)) {
      next
    }
    if (grepl("-", piece, fixed = TRUE)) {
      bounds <- trimws(unlist(strsplit(piece, "-", fixed = TRUE)))
      if (length(bounds) == 2) {
        start <- suppressWarnings(as.integer(bounds[[1]]))
        stop_value <- suppressWarnings(as.integer(bounds[[2]]))
        if (!is.na(start) && !is.na(stop_value) && stop_value >= start) {
          values <- c(values, seq.int(start, stop_value))
        }
      }
    } else {
      single_value <- suppressWarnings(as.integer(piece))
      if (!is.na(single_value)) {
        values <- c(values, single_value)
      }
    }
  }

  unique(values)
}

pick_dimensionality_reduction <- function(reduction_names) {
  if (length(reduction_names) == 0) {
    return("")
  }
  lowered <- tolower(reduction_names)
  if (any(grepl("pca", lowered))) {
    return("PCA")
  }
  if (any(grepl("umap", lowered))) {
    return("UMAP")
  }
  if (any(grepl("tsne", lowered))) {
    return("t-SNE")
  }
  trim_or_empty(reduction_names)
}

obj <- readRDS(rds_path)
meta <- obj[[]]
warnings <- character()

keep <- rep(TRUE, nrow(meta))
requested_fovs <- parse_fovs(fov_spec)
if (length(requested_fovs) > 0) {
  if (!"fov" %in% colnames(meta)) {
    stop("The Seurat object does not contain an fov column for CosMx subsetting")
  }
  fov_values <- trimws(as.character(meta$fov))
  requested_labels <- as.character(requested_fovs)
  missing_fovs <- setdiff(requested_labels, unique(fov_values))
  if (length(missing_fovs) > 0) {
    warnings <- c(warnings, sprintf("Requested FOVs were not present in the object: %s", paste(missing_fovs, collapse = ",")))
  }
  keep <- fov_values %in% requested_labels
}

if (!any(keep)) {
  stop("No cells matched the requested FOV subset")
}

meta_subset <- meta[keep, , drop = FALSE]
assay_name <- DefaultAssay(obj)
assay_names <- Assays(obj)
feature_count <- nrow(obj[[assay_name]])
panel_values <- if ("Panel" %in% colnames(meta_subset)) unique_nonempty(meta_subset$Panel) else character()
chemistry_values <- if ("version" %in% colnames(meta_subset)) unique_nonempty(meta_subset$version) else character()
slide_values <- if ("slide_ID_numeric" %in% colnames(meta_subset)) unique_nonempty(meta_subset$slide_ID_numeric) else character()
segmentation_values <- if ("cellSegmentationSetName" %in% colnames(meta_subset)) unique_nonempty(meta_subset$cellSegmentationSetName) else character()
mean_columns <- colnames(meta_subset)[startsWith(colnames(meta_subset), "Mean.")]
channels <- gsub("^Mean\\.", "", mean_columns)
reduction_names <- Reductions(obj)

read_counts <- if ("nCount_RNA" %in% colnames(meta_subset)) suppressWarnings(as.numeric(meta_subset$nCount_RNA)) else numeric()

result <- list(
  assay_chemistry_version = trim_or_empty(chemistry_values),
  cell_segmentation_method = trim_or_empty(segmentation_values),
  cell_segmented_object_type = "Whole cell",
  dimensionality_reduction_method = pick_dimensionality_reduction(reduction_names),
  has_cell_segmentation = TRUE,
  has_dimensionality_reduction = length(reduction_names) > 0,
  number_of_segmented_cells = as.character(nrow(meta_subset)),
  panel_name = paste(panel_values, collapse = "; "),
  panel_size_total_targets = as.character(feature_count),
  protein_measured = any(grepl("protein|adt", assay_names, ignore.case = TRUE)),
  qc_feature_number = as.character(nrow(meta_subset)),
  qc_mean_reads_per_feature = if (length(read_counts) > 0) fmt_num(mean(read_counts, na.rm = TRUE)) else "",
  qc_spatial_unit = "cell",
  qc_total_genes_detected = as.character(feature_count),
  qc_total_number_of_reads = if (length(read_counts) > 0) fmt_num(sum(read_counts, na.rm = TRUE)) else "",
  rna_measured = any(grepl("rna", assay_names, ignore.case = TRUE)) || identical(tolower(assay_name), "rna"),
  same_section_imaging_channels = paste(channels, collapse = ","),
  slide_serial_number = trim_or_empty(slide_values),
  software_and_version = "",
  spatial_assay_type = "In situ",
  transcriptome_type = "Targeted",
  warnings = unname(warnings)
)

cat(toJSON(result, auto_unbox = TRUE))
