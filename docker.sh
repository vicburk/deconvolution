#!/bin/bash

# Input file with values
INPUT_FILE="sample_size.txt"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Input file $INPUT_FILE not found"
  exit 1
fi

# Prompt for username and token
read -p "Enter CibersortX username: " USERNAME
read -p "Enter CibersortX token: " TOKEN
echo # new line after token input

# Export variables so they are available to the parallel command
export USERNAME
export TOKEN

if ! command -v parallel &>/dev/null; then
  echo "parallel is not installed, installing..."
  if [ -x "$(command -v apt-get)" ]; then
    sudo apt-get update && sudo apt-get install -y parallel
  elif [ -x "$(command -v yum)" ]; then
    sudo yum install -y parallel
  elif [ -x "$(command -v brew)" ]; then
    brew install parallel
  else
    echo "Unsupported Linux distribution"
    exit 1
  fi
fi

parallel --jobs 4 \
  "mkdir -p $HOME/deconvolution/pseudo_output_{}; \
  docker run \
  -v $HOME/deconvolution:/src/data \
  -v $HOME/deconvolution/pseudo_output_{}:/src/outdir cibersortx/fractions \
  --username \"$USERNAME\" \
  --token \"$TOKEN\" \
  --single_cell TRUE \
  --refsample scRNA_combined_{}.txt \
  --mixture pseudobulk_scRNA_combined_{}.txt \
  --perm 100" \
  <"$INPUT_FILE"
