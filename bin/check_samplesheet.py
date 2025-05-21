#!/usr/bin/env python

import os
import sys
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Check and validate input samplesheet")
    parser.add_argument("samplesheet", help="Input samplesheet CSV file")
    parser.add_argument("output", help="Output CSV file")
    return parser.parse_args()

def check_samplesheet(file_in, file_out):
    """
    Check that the samplesheet follows the expected format:
    sample,fastq_1,fastq_2,vcf
    """
    required_columns = ["sample", "fastq_1", "fastq_2", "vcf"]
    optional_columns = ["strandedness"]

    with open(file_in, "r") as fin:
        reader = csv.DictReader(fin)
        if not reader.fieldnames:
            print("ERROR: Input samplesheet is empty or malformed", file=sys.stderr)
            sys.exit(1)

        # Check header
        if not all(col in reader.fieldnames for col in required_columns):
            print(f"ERROR: Input samplesheet must contain these column headers: {', '.join(required_columns)}", file=sys.stderr)
            sys.exit(1)

        # Check each row
        rows = []
        for i, row in enumerate(reader):
            # Check required fields
            for col in required_columns:
                if not row[col]:
                    print(f"ERROR: Missing value for '{col}' in row {i+1}", file=sys.stderr)
                    sys.exit(1)

            # Check file existence
            for col in ["fastq_1", "fastq_2", "vcf"]:
                if not os.path.exists(row[col]):
                    print(f"WARNING: File does not exist: {row[col]}", file=sys.stderr)

            # Set default strandedness if not provided
            if "strandedness" in reader.fieldnames and not row["strandedness"]:
                row["strandedness"] = "unstranded"

            # Add row to output
            rows.append(row)

    # Write validated samplesheet
    with open(file_out, "w") as fout:
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Create meta file for each sample - write to current directory
    for row in rows:
        sample = row["sample"]
        meta_file = f"{sample}.meta.csv"  # Write to current directory
        with open(meta_file, "w") as fout:
            writer = csv.DictWriter(fout, fieldnames=["id", "single_end"])
            writer.writeheader()
            writer.writerow({"id": sample, "single_end": "false"})

def main():
    args = parse_args()
    check_samplesheet(args.samplesheet, args.output)

if __name__ == "__main__":
    main()
