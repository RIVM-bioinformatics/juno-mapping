#!/bin/bash

set -euo pipefail

set -x

#----------------------------------------------#
# User parameters
if [ ! -z "${1}" ] || [ ! -z "${2}" ] || [ ! -z "${irods_input_projectID}" ]
then
   input_dir="${1}"
   output_dir="${2}"
   PROJECT_NAME="${irods_input_projectID}"
   EXCLUSION_FILE=""
else
    echo "One of the parameters is missing, make sure there is an input directory, output directory and project name(param 1, 2 or irods_input_projectID)."
    exit 1
fi

set +u
#check if there is an exclusion file, if so change the parameter
if [ ! -z "${irods_input_sequencing__run_id}" ] && [ -f "/data/BioGrid/NGSlab/sample_sheets/${irods_input_sequencing__run_id}.exclude" ]
then
  EXCLUSION_FILE="/data/BioGrid/NGSlab/sample_sheets/${irods_input_sequencing__run_id}.exclude"
fi
set -u

if [ ! -d "${input_dir}" ] || [ ! -d "${output_dir}" ]
then
  echo "The input directory $input_dir, output directory $output_dir or fastq dir ${input_dir}/clean_fastq does not exist"
  exit 1
else
  input_fastq="${input_dir}/clean_fastq"
fi

#----------------------------------------------#
## make sure conda works

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/mnt/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/mnt/miniconda/etc/profile.d/conda.sh" ]; then
        . "/mnt/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/mnt/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __conda_hashr


#----------------------------------------------#
# Create the environment

# we can use the base installation of mamba to create the environment. 
# Swapping to a parent env is not necessary anymore.
mamba env create -f envs/juno_mapping.yaml --name pipeline_env
conda activate pipeline_env

#----------------------------------------------#
# Run the pipeline

case $PROJECT_NAME in
  adhoc)
    SPECIES="NotProvided"
    ;;
  myco|myco_kncv)
    SPECIES="Mycobacterium_tuberculosis"
    ;;
  *)
    SPECIES="NotProvided"
    ;;
esac

echo -e "\nRun pipeline..."

if [ ! -z ${irods_runsheet_sys__runsheet__lsf_queue} ]; then
    QUEUE="${irods_runsheet_sys__runsheet__lsf_queue}"
else
    QUEUE="bio"
fi


# Setting up the tmpdir for singularity as the current directory (default is /tmp but it gets full easily)
# Containers will use it for storing tmp files when building a container
export SINGULARITY_TMPDIR="$(pwd)"

#without exclusion file
# set time limit to 600 minutes in case large files are present. GATK is slow.
if [ "${EXCLUSION_FILE}" == "" ]
then
    python juno_mapping.py \
        --queue "${QUEUE}" \
        -i "${input_dir}" \
        -o "${output_dir}" \
        -s "${SPECIES}" \
        --prefix "/mnt/db/juno/sing_containers" \
        -tl 600

        result=$?
else
    python juno_mapping.py \
        --queue "${QUEUE}" \
        -i "${input_dir}" \
        -o "${output_dir}" \
        -s "${SPECIES}" \
        --prefix "/mnt/db/juno/sing_containers" \
        -ex "${EXCLUSION_FILE}" \
        -tl 600

        result=$?
fi

# Propagate metadata

set +euo pipefail

SEQ_KEYS=
SEQ_ENV=`env | grep -E 'irods_input_(sequencing|user)'`
for SEQ_AVU in ${SEQ_ENV}
do
    SEQ_KEYS="${SEQ_KEYS} ${SEQ_AVU%%=*}"
done

for key in $SEQ_KEYS irods_input_illumina__Flowcell irods_input_illumina__Instrument \
    irods_input_illumina__Date irods_input_illumina__Run_number irods_input_illumina__Run_Id
do
    if [ ! -z ${!key} ] ; then
        attrname=${key:12}
        attrname=${attrname/__/::}
        echo "${attrname}: '${!key}'" >> ${output_dir}/metadata.yml
    fi
done

set -euo pipefail

exit ${result}
