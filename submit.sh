#!/bin/bash
# JOB HEADERS HERE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=48:00:00

create_bam.sh
python ./script/main.py --input_dir=L5_NA24385_test --correct_bam_output_dir=correct_bam_files_test --output_dir=less_covered_regions_test --merge_supplementary=True --num_threads=8
