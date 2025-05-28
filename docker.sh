#!/bin/bash

#!/bin/bash

if ! command -v parallel &> /dev/null; then
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

parallel --jobs 4 "mkdir -p /home/vibu/deconvolution/output_{}; docker run -v /home/vibu/deconvolution:/src/data -v /home/vibu/deconvolution/output_{}:/src/outdir cibersortx/fractions --username burklovt@mail.uc.edu --token 55cf29498748253a278dbbb834eec37f --single_cell TRUE --refsample scRNA_combined_{}.txt --mixture gene_expression.txt --fraction 0 --rmbatchSmode TRUE" ::: 500 1000 1500 2000