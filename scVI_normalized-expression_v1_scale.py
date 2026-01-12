import anndata
import scvi
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="SCVI Normalized Expression Processing")
    parser.add_argument('--model_path', type=str, required=True)
    parser.add_argument('--adata_path', type=str, required=True)
    parser.add_argument('--output_adata_path', type=str, required=True)
    return parser.parse_args()

args = parse_args()

# Load the AnnData object
all_data = anndata.read_h5ad(args.adata_path)

# Load the model
model = scvi.model.SCVI.load(args.model_path, all_data)

# Retrieve and store normalized expression data
normalized_expression = model.get_normalized_expression(library_size=1e4, n_samples=25)
all_data.layers["scVI_normalized"] = normalized_expression

# Save the AnnData object
all_data.write_h5ad(args.output_adata_path)
