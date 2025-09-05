#!/usr/bin/env python3

"""
ichorCNA Results Summary Script

Usage:
  python ichorCNA_results_summary.py \
      --results_dir /path/to/results \
      --output_dir /path/to/output \
      [--bam_name_pattern ".bam"] \
      [--create_zips]

Description:
  This script extracts tumor fraction (TF) data and copy number alteration (CNA)
  data from ichorCNA results, performs basic transformations, and outputs several
  text files. Optionally, it can also create ZIP archives of parameter and CNA
  segmentation files.

Required Arguments:
  --results_dir        Path to the directory containing ichorCNA results.
  --output_dir         Path to the directory for output files.

Optional Arguments:
  --bam_name_pattern   Pattern in file names to remove in output. Defaults to ".bam".
  --create_zips        If provided, the script will also create 'params.zip'
                       and 'cna_seg.zip' inside --output_dir.

Example:
  python ichorCNA_post_analysis.py \
      --results_dir /path/to/results \
      --output_dir /path/to/output \
      --bam_name_pattern "_filt.bam" \
      --create_zips
"""

import os
import argparse
import pandas as pd
import numpy as np
import zipfile

def extract_tf_data(results_dir):
    """
    Extracts relevant data from ichorCNA results' params.txt files.

    Parameters:
    - results_dir (str): Path to the directory containing ichorCNA results folders.

    Returns:
    - pd.DataFrame: A Pandas DataFrame containing extracted data.
    """
    columns = ["library", "tumor_fraction", "ploidy", "gender", "ChrY_coverage_fraction", "ChrX_median_log_ratio"]
    data_list = []
    
    for sample_folder in os.listdir(results_dir):
        sample_path = os.path.join(results_dir, sample_folder)
        
        if os.path.isdir(sample_path):
            params_file_path = None
            for file_name in os.listdir(sample_path):
                if file_name.endswith(".params.txt"):
                    params_file_path = os.path.join(sample_path, file_name)
                    break
            
            if not params_file_path:
                print(f"No params.txt file found for sample: {sample_folder}")
                continue
            
            data = {"library": sample_folder}
            with open(params_file_path, "r") as file:
                for line in file:
                    line = line.strip()
                    if "Gender:" in line:
                        data["gender"] = line.split(":")[1].strip()
                    elif "Tumor Fraction:" in line:
                        value = line.split(":")[1].strip()
                        data["tumor_fraction"] = float(value) if value != 'NA' else np.nan
                    elif "Ploidy:" in line:
                        value = line.split(":")[1].strip()
                        data["ploidy"] = float(value) if value != 'NA' else np.nan
                    elif "ChrY coverage fraction:" in line:
                        value = line.split(":")[1].strip()
                        data["ChrY_coverage_fraction"] = float(value) if value != 'NA' else np.nan
                    elif "ChrX median log ratio:" in line:
                        value = line.split(":")[1].strip()
                        data["ChrX_median_log_ratio"] = float(value) if value != 'NA' else np.nan
            
            data_list.append(data)
                    
    result_df = pd.DataFrame(data_list, columns=columns).sort_values('library').reset_index(drop=True)
    return result_df


def extract_cna_data(parent_directory, logR_column_choice="logR_Copy_Number"):
    """
    Extracts CNA data from .cna.seg files for each sample directory.

    Parameters:
    - parent_directory (str): Path to the directory containing subfolders with .cna.seg files.
    - logR_column_choice (str): Which logR column to extract (e.g., 'logR' or 'logR_Copy_Number').

    Returns:
    - pd.DataFrame: A long-format DataFrame with columns
                    [library, chr, start, end, <logR_column_choice>].
    """
    data = []

    for item in os.listdir(parent_directory):
        item_path = os.path.join(parent_directory, item)
        
        if os.path.isdir(item_path):
            for file_name in os.listdir(item_path):
                if file_name.endswith(".cna.seg"):
                    library = file_name.replace(".cna.seg", "")
                    cna_seg_file = os.path.join(item_path, file_name)
                    
                    if os.path.isfile(cna_seg_file):
                        df = pd.read_csv(cna_seg_file, sep="\t")
                        logR_column = f"{library}.{logR_column_choice}"
                        
                        if logR_column in df.columns:
                            required_columns = ['chr', 'start', 'end', logR_column]
                            extracted_data = df.loc[:, required_columns]
                            extracted_data = extracted_data.rename(columns={logR_column: logR_column_choice})
                            extracted_data['library'] = library
                            data.extend(
                                extracted_data[['library', 'chr', 'start', 'end', logR_column_choice]].values.tolist()
                            )

    combined_df = pd.DataFrame(
        data, 
        columns=['library', 'chr', 'start', 'end', logR_column_choice]
    ).sort_values(by=['library', 'chr', 'start', 'end']).reset_index(drop=True)

    return combined_df


def create_params_zip(results_dir, output_zip):
    """
    Creates a zip file containing all .params.txt files from ichorCNA results.

    Parameters:
    - results_dir (str): Path to the directory containing ichorCNA results folders.
    - output_zip (str): Path to the output zip file.
    """
    with zipfile.ZipFile(output_zip, 'w') as zipf:
        for sample_folder in os.listdir(results_dir):
            sample_path = os.path.join(results_dir, sample_folder)
            
            if os.path.isdir(sample_path):
                for file_name in os.listdir(sample_path):
                    if file_name.endswith(".params.txt"):
                        params_file_path = os.path.join(sample_path, file_name)
                        zipf.write(params_file_path, arcname=file_name)
                        break
                else:
                    print(f"No params.txt file found for sample: {sample_folder}")


def create_cna_seg_zip(results_dir, output_zip):
    """
    Creates a zip file containing all .cna.seg files from ichorCNA results.

    Parameters:
    - results_dir (str): Path to the directory containing ichorCNA results folders.
    - output_zip (str): Path to the output zip file.
    """
    with zipfile.ZipFile(output_zip, 'w') as zipf:
        for sample_folder in os.listdir(results_dir):
            sample_path = os.path.join(results_dir, sample_folder)
            
            if os.path.isdir(sample_path):
                for file_name in os.listdir(sample_path):
                    if file_name.endswith(".cna.seg"):
                        params_file_path = os.path.join(sample_path, file_name)
                        zipf.write(params_file_path, arcname=file_name)
                        break
                else:
                    print(f"No cna.seg file found for sample: {sample_folder}")


def main():
    parser = argparse.ArgumentParser(
        description="Post-analysis script for ichorCNA output.",
        usage="%(prog)s --results_dir RESULTS_DIR --output_dir OUTPUT_DIR [--bam_name_pattern PATTERN] [--create_zips]"
    )
    parser.add_argument("--results_dir", required=True,
                        help="Path to the directory containing ichorCNA results.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to which output files will be written.")
    parser.add_argument("--bam_name_pattern", default=".bam",
                        help="Pattern in file names to remove in output. Default: '.bam'.")
    parser.add_argument("--create_zips", action="store_true", default=False,
                        help="If provided, also create params.zip and cna_seg.zip in the output directory.")

    args = parser.parse_args()

    # Manual check for required args
    if not args.results_dir or not args.output_dir:
        parser.print_help()
        parser.error("\nError: Missing one or more required arguments.\n")

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Extract TF data
    tf_data = extract_tf_data(args.results_dir)
    tf_data['library'] = tf_data['library'].str.replace(args.bam_name_pattern, '', regex=False)

    # 2. Extract CNA data (logR)
    cna_data_logR = extract_cna_data(args.results_dir, logR_column_choice="logR")
    cna_data_logR['library'] = cna_data_logR['library'].str.replace(args.bam_name_pattern, '', regex=False)

    # 3. Extract CNA data (logR_Copy_Number)
    cna_data_logR_Copy_Number = extract_cna_data(args.results_dir, logR_column_choice="logR_Copy_Number")
    cna_data_logR_Copy_Number['library'] = cna_data_logR_Copy_Number['library'].str.replace(
        args.bam_name_pattern, '', regex=False)

    # 4. Pivot to form matrix versions
    cna_matrix_logR = cna_data_logR.pivot(index=["chr", "start", "end"],
                                         columns="library",
                                         values="logR").reset_index()
    cna_matrix_logR.columns.name = None

    cna_matrix_logR_Copy_Number = cna_data_logR_Copy_Number.pivot(index=["chr", "start", "end"],
                                                                  columns="library",
                                                                  values="logR_Copy_Number").reset_index()
    cna_matrix_logR_Copy_Number.columns.name = None

    # 5. Write out results
    tf_data.to_csv(os.path.join(args.output_dir, "tf.txt"), sep="\t", index=False)
    cna_data_logR.to_csv(os.path.join(args.output_dir, "cna_logR_long.txt"), sep="\t", index=False)
    cna_data_logR_Copy_Number.to_csv(os.path.join(args.output_dir, "cna_logR_Copy_Number_long.txt"), sep="\t", index=False)
    cna_matrix_logR.to_csv(os.path.join(args.output_dir, "cna_logR_matrix.txt"), sep="\t", index=False)
    cna_matrix_logR_Copy_Number.to_csv(os.path.join(args.output_dir, "cna_logR_Copy_Number_matrix.txt"), sep="\t", index=False)

    # 6. Create zips if requested
    if args.create_zips:
        create_params_zip(args.results_dir, os.path.join(args.output_dir, "params.zip"))
        create_cna_seg_zip(args.results_dir, os.path.join(args.output_dir, "cna_seg.zip"))


if __name__ == "__main__":
    main()