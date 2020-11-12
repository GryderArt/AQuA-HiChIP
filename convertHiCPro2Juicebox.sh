#!/usr/bin/env bash

module load hicpro/2.11.1
module load juicer

SCRIPT_PATH=/usr/local/Anaconda/envs_app/hicpro/2.11.1/HiC-Pro_2.11.1/bin/utils/hicpro2juicebox.sh
JUICER_TOOLS=/usr/local/apps/juicer/juicer-1.5.6/scripts/juicer_tools.jar
DATA_HOME=/data/khanlab/projects/HiC/projects
SAMPLE_PATH=${DATA_HOME}/${HIC_PROJECT}/HiCpro_OUTPUT/hic_results/data/${HIC_SAMPLE}/${HIC_SAMPLE}.allValidPairs
SAMPLE_FOLDER=${DATA_HOME}/${HIC_PROJECT}/HiCpro_OUTPUT/hic_results/data/${HIC_SAMPLE}/

echo "Genome: $GENOME"
echo "Script: $SCRIPT_PATH"
echo "Sample: $HIC_PROJECT"
echo "Sample: $HIC_SAMPLE"
echo "Sample folder: $SAMPLE_FOLDER"
echo "Sample path: $SAMPLE_PATH"

cd $SAMPLE_FOLDER
echo "$SCRIPT_PATH -i $SAMPLE_PATH -g $GENOME -j $JUICER_TOOLS"
$SCRIPT_PATH -i $SAMPLE_PATH -g $GENOME -j $JUICER_TOOLS