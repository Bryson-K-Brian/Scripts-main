import os
import readline
import requests
import subprocess
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

BASE_URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq"

# 🔹 Enable tab completion (Unix Only)
def complete_path(text, state):
    """Tab-complete file and directory paths."""
    matches = [x for x in os.listdir('.') if x.startswith(text)]
    return matches[state] if state < len(matches) else None

readline.set_completer(complete_path)
readline.parse_and_bind("tab: complete")

def prompt_for_file(prompt_text):
    """Prompt user for a file with tab completion."""
    while True:
        file_path = input(prompt_text).strip()
        if os.path.isfile(file_path):
            return file_path
        print("❌ File not found. Try again.")

def prompt_for_directory(prompt_text):
    """Prompt user for a directory with tab completion."""
    while True:
        dir_path = input(prompt_text).strip()
        if os.path.isdir(dir_path):
            return dir_path
        print("❌ Directory not found. Try again.")

def prompt_for_threads():
    """Prompt user for number of parallel downloads."""
    while True:
        try:
            threads = int(input("🔹 Enter number of parallel downloads (default 5): ") or "5")
            if threads > 0:
                return threads
        except ValueError:
            pass
        print("❌ Invalid input. Enter a positive number.")

def get_sra_accessions(file_path):
    """Extract SRA accession numbers from TXT or CSV files."""
    if file_path.endswith(".txt"):
        with open(file_path, "r") as file:
            return [line.strip() for line in file if line.strip()]
    elif file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
        print("\n📄 CSV File Columns Found:")
        for i, col in enumerate(df.columns):
            print(f"  {i+1}. {col}")
        while True:
            try:
                col_index = int(input("\n🔹 Select the column number containing SRA accessions: ")) - 1
                if 0 <= col_index < len(df.columns):
                    return df.iloc[:, col_index].dropna().astype(str).tolist()
            except ValueError:
                pass
            print("❌ Invalid selection. Try again.")
    else:
        print("❌ Unsupported file format. Use a .txt or .csv file.")
        exit(1)

def file_exists(url, output_path):
    """Check if the file exists and is complete."""
    if os.path.exists(output_path):
        remote_size = int(requests.head(url).headers.get("Content-Length", 0))
        local_size = os.path.getsize(output_path)
        if remote_size == local_size:
            return True  # File is fully downloaded
    return False

def download_with_resume(url, output_path):
    """Download file with resume support and progress bar."""
    if file_exists(url, output_path):
        return f"✅ {os.path.basename(output_path)} (Already Downloaded)"

    temp_path = output_path + ".part"
    existing_size = os.path.getsize(temp_path) if os.path.exists(temp_path) else 0
    headers = {"Range": f"bytes={existing_size}-"} if existing_size else {}

    try:
        with requests.get(url, headers=headers, stream=True) as response:
            response.raise_for_status()
            mode = "ab" if existing_size else "wb"
            total_size = int(response.headers.get("Content-Length", 0)) + existing_size

            with open(temp_path, mode) as file, tqdm(
                total=total_size, unit="B", unit_scale=True,
                desc=f"⬇️ {os.path.basename(output_path)}", initial=existing_size
            ) as progress:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
                    progress.update(len(chunk))

        os.rename(temp_path, output_path)
        return f"✅ {os.path.basename(output_path)} (Downloaded from ENA)"
    except requests.RequestException as e:
        return f"❌ Failed: {os.path.basename(output_path)} → {e}"

def download_from_ncbi(sra_id, output_dir):
    """Use NCBI's fasterq-dump if ENA fails."""
    output_path_r1 = os.path.join(output_dir, f"{sra_id}_1.fastq.gz")
    output_path_r2 = os.path.join(output_dir, f"{sra_id}_2.fastq.gz")

    if os.path.exists(output_path_r1) and os.path.exists(output_path_r2):
        return f"✅ {sra_id} (Already Downloaded from NCBI)"

    try:
        subprocess.run(
            ["fasterq-dump", "--split-files", "--gzip", "-O", output_dir, sra_id],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return f"✅ {sra_id} (Downloaded from NCBI)"
    except subprocess.CalledProcessError:
        return f"❌ {sra_id} (Failed on both ENA & NCBI)"

def process_sra(sra_id, output_dir):
    """Try ENA first, then fallback to NCBI if not found."""
    subdir = f"{sra_id[:6]}/{sra_id}"
    if len(sra_id) > 9:
        subdir = f"{sra_id[:6]}/00{sra_id[-1]}/{sra_id}"

    urls = [
        (f"{BASE_URL}/{subdir}/{sra_id}_1.fastq.gz", os.path.join(output_dir, f"{sra_id}_1.fastq.gz")),
        (f"{BASE_URL}/{subdir}/{sra_id}_2.fastq.gz", os.path.join(output_dir, f"{sra_id}_2.fastq.gz")),
    ]

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(lambda x: download_with_resume(*x), urls))

    if all("❌" in res for res in results):
        return download_from_ncbi(sra_id, output_dir)
    return f"🔄 {sra_id}: {results[0]}, {results[1]}"

def download_from_ena_or_ncbi(accession_list, output_dir, max_threads=5):
    """Download multiple SRA files in parallel."""
    os.makedirs(output_dir, exist_ok=True)

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        results = list(tqdm(executor.map(process_sra, accession_list, [output_dir] * len(accession_list)), total=len(accession_list), desc="🚀 Downloading SRA files"))

    print("\n".join(results))

# ** Interactive Prompts **
print("\n🚀 Welcome to the SRA Downloader! 🚀\n")

accession_file = prompt_for_file("🔹 Enter the SRA accession file (TXT/CSV, TAB to autocomplete): ")
output_dir = prompt_for_directory("🔹 Enter the output directory (TAB to autocomplete): ")
max_threads = prompt_for_threads()

accession_list = get_sra_accessions(accession_file)
download_from_ena_or_ncbi(accession_list, output_dir, max_threads)