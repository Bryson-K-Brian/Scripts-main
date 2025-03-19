import os
import subprocess
import argparse

def check_existing_files(output_dir, sra_id):
    """Check if the SRA file has already been downloaded and converted."""
    read1 = os.path.join(output_dir, f"{sra_id}_1.fastq")
    read2 = os.path.join(output_dir, f"{sra_id}_2.fastq")
    single_read = os.path.join(output_dir, f"{sra_id}.fastq")
    
    if os.path.exists(read1) and os.path.exists(read2):
        return True  # Paired-end files exist
    if os.path.exists(single_read):
        return True  # Single-end file exists
    return False  # Needs downloading

def run_fasterq_dump(sra_file, output_dir, threads):
    """Run fasterq-dump in parallel using GNU parallel."""
    
    # Read SRA IDs from the input file
    with open(sra_file, "r") as f:
        sra_list = [line.strip() for line in f if line.strip()]

    if not sra_list:
        print("❌ No valid SRA IDs found in the input file.")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Filter out already downloaded files
    sra_list = [sra for sra in sra_list if not check_existing_files(output_dir, sra)]
    
    if not sra_list:
        print("✅ All files are already downloaded.")
        return

    print(f"🚀 Starting parallel downloads for {len(sra_list)} SRA IDs...")

    # GNU parallel command
    parallel_cmd = f'parallel -j {threads} "fasterq-dump --split-files --progress --outdir {output_dir} {{}}" ::: {" ".join(sra_list)}'

    # Run the command
    try:
        subprocess.run(parallel_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during parallel execution: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SRA files in parallel using fasterq-dump.")
    parser.add_argument("sra_file", help="File containing SRA IDs (one per line).")
    parser.add_argument("output_dir", help="Directory to store downloaded FASTQ files.")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel downloads (default: 4).")

    args = parser.parse_args()
    run_fasterq_dump(args.sra_file, args.output_dir, args.threads)