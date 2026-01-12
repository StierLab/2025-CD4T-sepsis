import scanpy as sc
import anndata
import sys

def merge_files(file_paths):
    adatas = [sc.read_h5ad(file_path) for file_path in file_paths]
    return anndata.concat(adatas, join='outer')

if __name__ == "__main__":
    file_paths = sys.argv[1:]  # File paths passed as arguments
    merged_adata = merge_files(file_paths)
    merged_adata.write_h5ad('merged_adata.h5ad')
