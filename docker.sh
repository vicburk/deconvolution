#!/bin/bash

# Install Docker (if not already installed)
if ! command -v docker &> /dev/null
then
  echo "Docker is not installed. Installing..."
  curl -fsSL https://get.docker.com -o get-docker.sh
  sudo sh get-docker.sh
  rm get-docker.sh
  echo "Docker installed successfully."
fi


# Pull the CIBERSORTx Docker image
echo "Pulling cibersortx/fractions Docker image..."
docker pull cibersortx/fractions
echo "Docker image pulled successfully."


# Get user credentials securely
read -r -s -p "Enter your CIBERSORTx username: " username
echo ""
read -r -s -p "Enter your CIBERSORTx token: " token
echo ""


# Clone the repository (if needed - adapt this part if you have a different way of getting the input files)
# Example:
# git clone <repository_url> cibersortx_data
# cd cibersortx_data


# Get input and output directory paths
read -r -p "Enter the absolute path to the input directory: " input_dir
read -r -p "Enter the absolute path to the output directory: " output_dir


# Validate input - Check if directories exist
if [[ ! -d "$input_dir" ]]; then
  echo "Error: Input directory '$input_dir' does not exist."
  exit 1
fi

if [[ ! -d "$output_dir" ]]; then
  echo "Error: Output directory '$output_dir' does not exist."
  exit 1
fi



# Run CIBERSORTx Fractions (example command - adapt as needed)
docker run -v "$input_dir":/src/data -v "$output_dir":/src/outdir cibersortx/fractions \
  --username "$username" --token "$token" \
  --single_cell TRUE \
  --refsample /src/data/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt \  # Make sure these paths are correct relative to the input directory
  --mixture /src/data/Fig2b-WholeBlood_RNAseq.txt \                 # Make sure these paths are correct relative to the input directory
  --fraction 0 --rmbatchSmode TRUE

echo "CIBERSORTx Fractions execution completed." 