import os
import subprocess
import pandas as pd
import requests
import tqdm
import time
from concurrent.futures import ThreadPoolExecutor

def is_sra_accession_valid(sra_id):
    """
    Check if an SRA accession ID is valid by querying the NCBI E-utilities API.
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {'db': 'sra', 'term': sra_id, 'retmode': 'json'}
    
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        return int(data['esearchresult']['count']) > 0
    return False

def get_sra_file_size(sra_id):
    """
    Estimate the file size of an SRA accession from ENA database.
    """
    ena_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}&result=read_run&fields=read_count,base_count,size_MB&format=tsv"
    response = requests.get(ena_url)
    
    if response.status_code == 200:
        lines = response.text.split("\n")
        if len(lines) > 1:
            size_mb = lines[1].split("\t")[-1]
            return float(size_mb) if size_mb else None
    return None

def check_existing_files(output_dir, sra_id):
    """
    Check if the SRA file has already been downloaded or extracted.
    """
    sra_cache_path = os.path.expanduser(f"~/ncbi/public/sra/{sra_id}.sra")
    read1 = os.path.join(output_dir, f"{sra_id}_1.fastq.gz")
    read2 = os.path.join(output_dir, f"{sra_id}_2.fastq.gz")

    if os.path.exists(read1) and os.path.exists(read2):
        return "FASTQ"
    if os.path.exists(sra_cache_path):
        return "SRA"
    return None

def download_sra(sra):
    """
    Download the SRA file using prefetch with resume capability.
    """
    cmd = f"prefetch {sra}"
    subprocess.run(cmd, shell=True)

def extract_fastq(sra, output_dir):
    """
    Convert SRA to FASTQ using fasterq-dump with progress bar.
    """
    cmd = f"fasterq-dump --split-files --progress --outdir {output_dir} ~/ncbi/public/sra/{sra}.sra"
    subprocess.run(cmd, shell=True)

def run_parallel_downloads(sra_file, output_dir, threads):
    """
    Handles parallel downloads using prefetch and fasterq-dump.
    """
    df = pd.read_csv(sra_file)
    sra_list = df.iloc[:, 0].dropna().astype(str).tolist()

    if not sra_list:
        print("❌ No valid SRA IDs found in the file.")
        return
    
    os.makedirs(output_dir, exist_ok=True)

    # Separate valid & new SRA IDs
    sra_to_prefetch = [sra for sra in sra_list if check_existing_files(output_dir, sra) is None]
    sra_to_extract = [sra for sra in sra_list if check_existing_files(output_dir, sra) == "SRA"]

    # Progress bar for validation
    print("🔎 Validating SRA IDs...")
    valid_sra = [sra for sra in tqdm.tqdm(sra_to_prefetch) if is_sra_accession_valid(sra)]

    if not valid_sra and not sra_to_extract:
        print("✅ All files are already processed.")
        return

    # Step 1: Prefetch SRA files in parallel
    if valid_sra:
        print(f"🚀 Downloading {len(valid_sra)} SRA files...")
        with ThreadPoolExecutor(max_workers=threads) as executor:
            list(tqdm.tqdm(executor.map(download_sra, valid_sra), total=len(valid_sra)))

    # Step 2: Run fasterq-dump in parallel
    sra_ready_for_fasterq = valid_sra + sra_to_extract
    if sra_ready_for_fasterq:
        print(f"📂 Extracting FASTQ files for {len(sra_ready_for_fasterq)} SRA IDs...")
        with ThreadPoolExecutor(max_workers=threads) as executor:
            list(tqdm.tqdm(executor.map(lambda sra: extract_fastq(sra, output_dir), sra_ready_for_fasterq), total=len(sra_ready_for_fasterq)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Download and extract SRA files in parallel.")
    parser.add_argument("sra_file", help="CSV file containing SRA IDs.")
    parser.add_argument("output_dir", help="Directory to store downloaded FASTQ files.")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel downloads (default: 4).")

    args = parser.parse_args()
    run_parallel_downloads(args.sra_file, args.output_dir, args.threads)