import pandas as pd
import scanpy as sc
import anndata
import scvi
import argparse

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Post-merging scVI analysis script")
    parser.add_argument('--input_adata', type=str, required=True)
    parser.add_argument('--HVG', type=int)
    parser.add_argument('--batch_key', type=str)
    parser.add_argument('--continuous_covariate_keys')
    parser.add_argument('--categorical_covariate_keys')
    parser.add_argument('--gene_likelihood', type=str)
    parser.add_argument('--n_latent', type=int)
    parser.add_argument('--n_layers', type=int)
    parser.add_argument('--n_hidden', type=int)
    parser.add_argument('--dispersion', type=str)
    parser.add_argument('--check_val_every_n_epoch', type=int)
    parser.add_argument('--max_epochs', type=int)
    parser.add_argument('--early_stopping', type=bool)
    parser.add_argument('--early_stopping_patience', type=int)
    parser.add_argument('--early_stopping_monitor', type=str)
    parser.add_argument('--lr', type=float)
    parser.add_argument('--optimizer', type=str)
    parser.add_argument('--dropout_rate', type=float)
    parser.add_argument('--output_model_path')
    parser.add_argument('--output_adata_path')
    return parser.parse_args()

# Parse arguments
args = parse_args()

# Convert comma-separated strings to lists
continuous_covariates = args.continuous_covariate_keys.split(',') if args.continuous_covariate_keys else []
categorical_covariates = args.categorical_covariate_keys.split(',') if args.categorical_covariate_keys else []

# Load the merged AnnData object
all_data = sc.read_h5ad(args.input_adata)

# Preprocessing steps
# Ensure counts layer exists for scVI model setup
all_data.layers["counts"] = all_data.X.copy()
# Normalize and log transform
sc.pp.normalize_total(all_data, target_sum=1e4)
sc.pp.log1p(all_data)
# Store the raw data for later use
all_data.raw = all_data

# Identify highly variable genes
sc.pp.highly_variable_genes(
    all_data, 
    n_top_genes=args.HVG, 
    subset=True, 
    layer="counts", 
    flavor="seurat_v3", 
    batch_key=args.batch_key
)

# Setup the AnnData object for scVI
scvi.model.SCVI.setup_anndata(
    all_data, 
    layer="counts", 
    continuous_covariate_keys=continuous_covariates, 
    categorical_covariate_keys=categorical_covariates, 
    batch_key=args.batch_key
)

# Create and train the SCVI model
model = scvi.model.SCVI(
    all_data,
    gene_likelihood=args.gene_likelihood,
    n_latent=args.n_latent,
    n_layers=args.n_layers,
    n_hidden=args.n_hidden,
    dropout_rate=args.dropout_rate,
    dispersion=args.dispersion
)

model.train(
    check_val_every_n_epoch=args.check_val_every_n_epoch,
    max_epochs=args.max_epochs,
    early_stopping=args.early_stopping,
    early_stopping_patience=args.early_stopping_patience,
    early_stopping_monitor=args.early_stopping_monitor,
    plan_kwargs={'lr': args.lr, 'optimizer': args.optimizer}
)

# Save the trained model
model.save(args.output_model_path)

# Retrieve and store latent space representations
latent_representation = model.get_latent_representation()
all_data.obsm["X_scVI"] = latent_representation

# Save the updated AnnData object
all_data.write_h5ad(args.output_adata_path)

