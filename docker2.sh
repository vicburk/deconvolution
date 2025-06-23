#!/bin/bash

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

# Write script similar to downsample test but use
# "4d3469a7-339f-40b3-92a3-22f7043545f8"
# "26f6625b-e76c-490a-beb1-aea16933cd6d"
# "5c9ab5a5-04a9-4282-9320-3b4d7b95131c"

parallel --jobs 3 \
  "mkdir -p /home/vibu/deconvolution/output_{}; \
  docker run -v /home/vibu/deconvolution:/src/data \
  -v /home/vibu/deconvolution/output_{}:/src/outdir \
  cibersortx/fractions \
  --username \"$USERNAME\" \
  --token \"$TOKEN\" \
  --single_cell TRUE \
  --refsample scRNA_{}.txt \
  --mixture gene_expression.txt \
  --fraction 0 \
  --rmbatchSmode TRUE" \
  ::: "4d3469a7-339f-40b3-92a3-22f7043545f8" "26f6625b-e76c-490a-beb1-aea16933cd6d" "5c9ab5a5-04a9-4282-9320-3b4d7b95131c"
