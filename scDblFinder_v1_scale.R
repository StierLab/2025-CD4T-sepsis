library(Matrix)
library(Seurat)
library(scCustomize)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(DropletUtils)
library(scDblFinder)
library(reticulate)
library(sceasy)

# Set SeuratObject v5.0.0 to the historical v3 assay
options(Seurat.object.assay.version = "v3")

# Ensure reticulate uses the correct Conda environment
use_condaenv("myenv", required = TRUE)

# Define arguments
args <- commandArgs(trailingOnly = TRUE)

h5_file <- args[1]
metadata_csv <- args[2]
sample_name <- args[3]
output_seurat <- args[4]
output_adata <- args[5]
percent_mito_threshold <- as.numeric(args[6])
min_features <- as.numeric(args[7])
max_features <- as.numeric(args[8])

# Ensure that there are enough arguments
if (length(args) < 8) {
  stop("Not enough arguments provided to the script")
}

# Trim white spaces from paths
h5_file <- trimws(h5_file)
metadata_csv <- trimws(metadata_csv)
output_seurat <- trimws(output_seurat)
output_adata <- trimws(output_adata)

# Read the CellBender output .h5 file
count_matrix <- Read_CellBender_h5_Mat(
  file_name = h5_file,
  use.names = TRUE,
  unique.features = TRUE
)

# Create a Seurat Object with Orig.Ident as Sample Name
sobj <- CreateSeuratObject(counts = count_matrix)
sobj$orig.ident <- sample_name

# Define ribosomal and heme-associated gene patterns
ribosomal_pattern <- "^RPS|^RPL"
heme_pattern <- "^HBA|^HBB"

# Calculate the percentage of ribosomal gene content
sobj <- PercentageFeatureSet(sobj, pattern = ribosomal_pattern, col.name = "percent.ribosomal")

# Calculate the percentage of heme-associated gene content
sobj <- PercentageFeatureSet(sobj, pattern = heme_pattern, col.name = "percent.heme")

# Calculate the percentage of mitochondrial gene content
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")

# Filter cells based on mitochondrial content and min/max features
initial_cell_count <- ncol(sobj)
sobj <- subset(sobj, subset = percent.mt < percent_mito_threshold)
filtered_mito_cell_count <- initial_cell_count - ncol(sobj)
sobj <- subset(sobj, subset = nFeature_RNA > min_features & nFeature_RNA < max_features)
feature_filtered_count <- initial_cell_count - filtered_mito_cell_count - ncol(sobj)
total_cells_after_filtering <- ncol(sobj)

# Normalization and finding variable features
sobj <- SCTransform(sobj, method="glmGamPoi", vst.flavor="v2", verbose = FALSE)
sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
sobj <- ScaleData(sobj, verbose = FALSE)
sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)

# Remove unnecessary assays and sample information
DefaultAssay(sobj) <- "RNA"
sobj[["SCT"]] <- NULL
sobj[["pca"]] <- NULL

# Create SCE object and assign cluster values
sce_obj <- as.SingleCellExperiment(sobj)
clusters <- sce_obj$seurat_clusters

# Doublet detection, counting, and subsetting to retain singlets
sce_obj <- scDblFinder(sce_obj, clusters = clusters)
doublet_filtered_count <- sum(sce_obj$scDblFinder.class == "doublet")
sce_obj <- sce_obj[, sce_obj$scDblFinder.class == "singlet"]

# Add filtering information to metadata
sce_obj$mito_filtered_cells <- filtered_mito_cell_count
sce_obj$feature_filtered_cells <- feature_filtered_count
sce_obj$doublet_filtered_cells <- doublet_filtered_count
sce_obj$percent_mito_filtered <- (filtered_mito_cell_count / initial_cell_count) * 100
sce_obj$percent_feature_filtered <- (feature_filtered_count / (initial_cell_count - filtered_mito_cell_count)) * 100
sce_obj$percent_doublets <- (doublet_filtered_count / total_cells_after_filtering) * 100

# Define a vector of colData names to remove
columns_to_remove <- c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters",
                       "scDblFinder.cluster", "scDblFinder.class", "scDblFinder.score",
                       "scDblFinder.weighted", "scDblFinder.difficulty", "scDblFinder.cxds_score",
                       "scDblFinder.mostLikelyOrigin", "scDblFinder.originAmbiguous", "ident", "filtered_mito_cell_count", "doublet_filtered_count", "feature_filtered_count", "mito_filtered_cells", "feature_filtered_cells", "doublet_filtered_cells")

# Remove the specified columns from colData
colData(sce_obj) <- colData(sce_obj)[, !colnames(colData(sce_obj)) %in% columns_to_remove]

# Add metadata
metadata <- read.csv(metadata_csv, row.names = 1, header = TRUE)
sce_sample_ids <- sce_obj$orig.ident
metadata_sample_ids <- rownames(metadata)
colData(sce_obj) <- cbind(colData(sce_obj), metadata[sce_obj$orig.ident, ])

# Dynamically set the output file names based on sample_name
output_seurat <- paste0(sample_name, ".rds")
output_adata <- paste0(sample_name, ".h5ad")

# Extract counts matrix
counts <- counts(sce_obj)

# Extract cell metadata
cell_metadata <- colData(sce_obj)

# Extract gene metadata
gene_metadata <- rowData(sce_obj)

# Save Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, 
                                 meta.data = as.data.frame(cell_metadata),
                                 row.names = rownames(gene_metadata))
saveRDS(seurat_obj, file = output_seurat)

# Function convert SCE to data using sceasy raw code
regularise_df <- function(df, drop_single_values = FALSE) {
  if (ncol(df) == 0) {
    df[["name"]] <- rownames(df)
  }
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0) {
      warning(
        "Dropping single category variables:",
        paste(colnames(df)[k_singular], collapse = ", ")
      )
      df <- df[, !k_singular, drop = FALSE]
    }
    if (ncol(df) == 0) {
      df[["name"]] <- rownames(df)
    }
  }
  return(df)
}

sce2anndata <- function(obj, outFile = NULL, main_layer = "counts", transfer_layers = NULL, drop_single_values = FALSE) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("This function requires the 'SummarizedExperiment' package.")
  }
  if (!requireNamespace("SingleCellExperiment")) {
    stop("This function requires the 'SingleCellExperiment' package.")
  }
  assay_names <- SummarizedExperiment::assayNames(obj)
  main_layer <- match.arg(main_layer, assay_names)
  transfer_layers <- transfer_layers[transfer_layers %in% assay_names]
  transfer_layers <- transfer_layers[transfer_layers != main_layer]
  
  X <- SummarizedExperiment::assay(obj, main_layer)
  
  obs <- regularise_df(as.data.frame(SummarizedExperiment::colData(obj)), drop_single_values = drop_single_values)
  
  var <- regularise_df(as.data.frame(SummarizedExperiment::rowData(obj)), drop_single_values = drop_single_values)
  
  obsm <- NULL
  reductions <- SingleCellExperiment::reducedDimNames(obj)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) {
        as.matrix(
          SingleCellExperiment::reducedDim(obj, type = name)
        )
      },
      simplify = FALSE
    )
    names(obsm) <- paste0(
      "X_", tolower(SingleCellExperiment::reducedDimNames(obj))
    )
  }
  
  layers <- list()
  for (layer in transfer_layers) {
    mat <- SummarizedExperiment::assay(obj, layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }
  
  anndata <- reticulate::import("anndata", convert = FALSE)
  
  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )
  
  if (!is.null(outFile)) {
    adata$write(outFile, compression = "gzip")
  }
  
  adata
}

# Convert SCE to AnnData file using raw code from sceasy as above
adata_obj <- sce2anndata(sce_obj, outFile = output_adata)

# List all objects in the environment
all_objects <- ls()

# Define objects to keep (file paths, not R objects)
objects_to_keep <- c("output_seurat", "output_adata")

# Remove all objects except the specified ones
rm(list=setdiff(all_objects, objects_to_keep))