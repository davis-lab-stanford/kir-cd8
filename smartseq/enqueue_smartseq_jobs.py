#!/usr/bin/python
# usage: ./enqueue_smartseq_jobs.py -h
# writes to jobs.json

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
    help="path to parent folder of *_R1_001.fastq.gz and *_R2_001.fastq.gz files",
)
args = parser.parse_args()

print("Processing %s FASTQ file paths." % len(args.fastq_file_paths))

# output will be to: output_path/[sample group name (one of fastq_file_paths)]/[samplename].{out|err}
# same for logs rooted at log_path

###

f1_pattern = "_R1_001.fastq.gz"
f2_pattern = "_R2_001.fastq.gz"

jobs = []

for fastq_file_path in args.fastq_file_paths:
    # add trailing slash
    fastq_file_path = os.path.join(fastq_file_path, "")

    # extract name of last folder
    sample_group_name = os.path.basename(os.path.normpath(fastq_file_path))

    # iterate over all [samplename]_R1 files in this group of samples
    for files in glob.glob(fastq_file_path + "*%s" % f1_pattern):
        # extract sample name
        file_name = files.split("/")[-1]
        sample_name = str(file_name)[: -len(f1_pattern)]

        # specify input files
        fastq_gz_file1 = os.path.join(fastq_file_path, sample_name + f1_pattern)
        fastq_gz_file2 = os.path.join(fastq_file_path, sample_name + f2_pattern)

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

        # configure job
        job_details = {
            "stdout_file": stdout_file,
            "stderr_file": stderr_file,
            "sample_output_directory": sample_output_directory,
            "fastq_gz_file1": fastq_gz_file1,
            "fastq_gz_file2": fastq_gz_file2,
            "method": "job_quantification.sh",
        }
        jobs.append(job_details)

        # run sbatch
        cmd = """sbatch \
            --account=mmdavis --partition=batch \
            --output={stdout_file} --error={stderr_file} \
            --export=wDir=\"{sample_output_directory}\",input_file1=\"{fastq_gz_file1}\",input_file2=\"{fastq_gz_file2}\" \
            --chdir={sample_output_directory} \
            {method};
        """.format(
            **job_details
        )
        print(cmd)
        if not args.dry_run:
            os.system(cmd)

with open(args.jobs_json_path, "w") as json_w:
    json.dump(jobs, json_w)

print("Complete.")
