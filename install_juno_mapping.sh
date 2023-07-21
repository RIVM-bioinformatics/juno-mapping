# Make sure to run from the main directory of the pipeline
set -euo pipefail
path_to_main=$(realpath $0)
workdir=$(dirname "$path_to_main")
cd "$workdir"

# Install mamba if needed
if ! command -v mamba &> /dev/null
then
    set +euo pipefail
    echo "mamba could not be found, installing now..."
    conda env update -f envs/mamba.yaml
    source activate mamba
    set -euo pipefail    
fi

# Delete previous installations
rm -rf envs/src/

# Install main env
mamba env update -f envs/juno_mapping.yaml