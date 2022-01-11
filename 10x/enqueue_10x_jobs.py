#!/usr/bin/python
# usage: ./enqueue_10x_jobs.py -h
# # writes to jobs.json

import json
import glob
import sys
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--log-path", required=True, help="directory for log file output")
parser.add_argument("--output-path", required=True, help="directory for job output")
parser.add_argument(
    "--jobs-json-path",
    required=True,
    help="file path to write json describing all scheduled jobs",
)
parser.add_argument(
    "--dry-run",
    action="store_true",
    default=False,
    help="don't create directories or actually run the jobs",
)
parser.add_argument(
    "fastq_file_paths",
    type=str,
    nargs="+",
    help="path to parent folder of sample-specific folders that contain fastqs",
)
args = parser.parse_args()

print("Processing %s FASTQ file paths." % len(args.fastq_file_paths))

###

jobs = []

for fastq_file_path in args.fastq_file_paths:
    # add trailing slash
    fastq_file_path = os.path.join(fastq_file_path, "")

    # extract name of last folder
    sample_group_name = os.path.basename(os.path.normpath(fastq_file_path))

    # iterate over all subdirectories in this group of samples
    for sample_input_directory in glob.glob(fastq_file_path + "*"):
        # extract sample name, e.g. 08162018-WTA
        sample_name = sample_input_directory.split("/")[-1]

        # specify job output location
        # create output directory specifically for this sample
        sample_output_directory = os.path.join(
            args.output_path, sample_group_name, sample_name
        )
        if not args.dry_run and not os.path.exists(sample_output_directory):
            # create intermediate directories as well
            os.makedirs(sample_output_directory)

        # specify slurm log filepaths
        # create logs directory specifically for this sample
        sample_logs_directory = os.path.join(args.log_path, sample_group_name)
        if not args.dry_run and not os.path.exists(sample_logs_directory):
            # create intermediate directories as well
            os.makedirs(sample_logs_directory)
        stdout_file = os.path.join(sample_logs_directory, sample_name + ".out")
        stderr_file = os.path.join(sample_logs_directory, sample_name + ".err")

        # Configure job
        job_details = {
            "stdout_file": stdout_file,
            "stderr_file": stderr_file,
            "sample_input_directory": sample_input_directory,
            "sample_output_directory": sample_output_directory,
            "sample_name": sample_name,
            "script_name": "job_count.sh",
        }
        jobs.append(job_details)

        # run sbatch
        cmd = """sbatch \
            --account=mmdavis --partition=batch \
            --output={stdout_file} --error={stderr_file} \
            --export=sample_name=\"{sample_name}\",fastq_path=\"{sample_input_directory}\" \
            --chdir={sample_output_directory} \
            {script_name};
        """.format(
            **job_details
        )
        print(cmd)
        if not args.dry_run:
            os.system(cmd)

with open(args.jobs_json_path, "w") as json_w:
    json.dump(jobs, json_w)

print("Complete.")
