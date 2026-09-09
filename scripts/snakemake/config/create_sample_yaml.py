#!/usr/bin/env python
# coding: utf-8

import pathlib
import os
import sys
import argparse

# --- Parse command-line arguments ---
parser = argparse.ArgumentParser(description="Generate a samples.yaml file from a folder of BAM files.")
parser.add_argument("bamfolder", help="Path to the folder containing BAM files")
parser.add_argument(
    "--extension",
    default=None,
    help="Suffix to strip from filenames to get sample names (e.g. 'GRCh38.bwa_meth.coorsort.filt.bam')"
)
args = parser.parse_args()

bamfolder = args.bamfolder
extension = args.extension

# --- Find and sort BAM files ---
bam_files = []
for filepath in pathlib.Path(bamfolder).glob('**/*.bam'):
    bam_files.append(str(filepath.absolute()))

bam_files.sort(key=lambda x: os.path.basename(x))

# --- Helper function to derive sample name from filename ---
def get_sample_name(bam_path, extension):
    basename = os.path.basename(bam_path)

    if extension is not None:
        if basename.endswith(extension):
            return basename[:-len(extension)]
        else:
            print(
                f"WARNING: '{basename}' does not end with expected extension "
                f"'{extension}'. Falling back to stripping only '.bam'.",
                file=sys.stderr
            )

    # Fallback: strip ".bam"
    if basename.endswith(".bam"):
        return basename[:-4]

    return basename

# --- Write samples.yaml ---
with open("samples.yaml", "w") as f:
    f.write('samples:\n')
    for bam_file in bam_files:
        sample_name = get_sample_name(bam_file, extension)
        f.write(f" {sample_name}: {bam_file}\n")

