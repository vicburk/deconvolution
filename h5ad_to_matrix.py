import scanpy as sc
import scipy.io as io
import subprocess
import glob
from scipy.sparse import issparse, csr_matrix
import os
import gc
import sys


def command_line(command, shell=False):
    result = subprocess.run(command, capture_output=True, text=True, shell=shell)
    if result.returncode == 0:
        print(result.stdout)
    else:
        print(f"Error: {result.stderr}")


# Determine which files to process
if len(sys.argv) > 1:
    # If command-line arguments are provided, use those
    h5ad_files = sys.argv[1:]
else:
    # Otherwise, use glob to find all .h5ad files
    h5ad_files = glob.glob("*.h5ad")

# Process each .h5ad file one by one
for file in h5ad_files:
    print(f"Processing {file}...")
    # Extract filename without extension
    prefix = os.path.splitext(os.path.basename(file))[0]
    # Create a separate directory for each file
    output_dir = f"matrix_files/{prefix}/"
    os.makedirs(output_dir, exist_ok=True)
    # Load file in read mode to reduce memory usage
    adata = sc.read_h5ad(file, backed="r")
    # Write barcodes.tsv
    with open(f"{output_dir}/barcodes.tsv", "w") as f:
        for item in adata.obs_names:
            f.write(item + "\n")
    # Write features.tsv
    with open(f"{output_dir}/features.tsv", "w") as f:
        for item in ["\t".join([x, x, "Gene Expression"]) for x in adata.var_names]:
            f.write(item + "\n")
    # Extract matrix without loading everything into memory
    matrix = adata[:, :].X.T  # Access without full memory load
    if not issparse(matrix):
        matrix = csr_matrix(matrix)
    io.mmwrite(f"{output_dir}/matrix", matrix)
    # Save metadata
    adata.obs.to_csv(f"{output_dir}/metadata.csv")
    # Free memory
    del adata
    gc.collect()

# List created directories
command_line(["ls", "matrix_files"])

# Gzip files only for processed h5ad files
for file in h5ad_files:
    prefix = os.path.splitext(os.path.basename(file))[0]
    output_dir = f"matrix_files/{prefix}/"

    # List of files to gzip
    files_to_gzip = [
        f"{output_dir}/barcodes.tsv",
        f"{output_dir}/features.tsv",
        f"{output_dir}/matrix.mtx",
        f"{output_dir}/metadata.csv",
    ]

    # Gzip each file
    for file_path in files_to_gzip:
        command_line(["gzip", file_path])

# List files after compression
command_line(["ls", "-R", "matrix_files"])
