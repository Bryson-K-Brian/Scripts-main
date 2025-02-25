import os
import glob
import requests
import readline
from tqdm import tqdm

# Enable tab completion for file and directory paths
def complete_path(text, state):
    return (glob.glob(text + '*') + [None])[state]

readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete_path)

def get_file_size(url):
    """Check the total file size on the server."""
    try:
        response = requests.head(url, allow_redirects=True)
        if "Content-Length" in response.headers:
            return int(response.headers["Content-Length"])
    except requests.RequestException:
        return None
    return None

def download_with_resume(url, output_path):
    """Download a file with resume capability and progress bar."""
    total_size = get_file_size(url)
    temp_path = output_path + ".part"

    if os.path.exists(output_path):
        # Skip download if the file is fully downloaded
        if os.path.getsize(output_path) == total_size:
            print(f"✅ File already exists: {output_path} (Skipping)")
            return
        else:
            print(f"⚠️ File is incomplete, resuming: {output_path}")

    existing_size = os.path.getsize(temp_path) if os.path.exists(temp_path) else 0

    headers = {"Range": f"bytes={existing_size}-"} if existing_size else {}

    try:
        with requests.get(url, headers=headers, stream=True) as response:
            response.raise_for_status()
            mode = "ab" if existing_size else "wb"

            with open(temp_path, mode) as file, tqdm(
                total=total_size, unit="B", unit_scale=True,
                desc=f"⬇️ {os.path.basename(output_path)}",
                initial=existing_size, ascii=True
            ) as progress:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
                    progress.update(len(chunk))

        os.rename(temp_path, output_path)  # Rename after successful download
        print(f"✅ Downloaded: {output_path}")
    except requests.RequestException as e:
        print(f"❌ Failed to download {url}: {e}")

def download_from_ena(accession_file, output_dir):
    base_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq"

    # Read SRA IDs from file
    with open(accession_file, "r") as file:
        accession_list = [line.strip() for line in file]

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    for sra_id in accession_list:
        # Determine ENA directory structure
        subdir = f"{sra_id[:6]}/{sra_id}"
        if len(sra_id) > 9:  # Handles long IDs (SRR123456789)
            subdir = f"{sra_id[:6]}/00{sra_id[-1]}/{sra_id}"

        # File URLs for Read 1 and Read 2
        url_r1 = f"{base_url}/{subdir}/{sra_id}_1.fastq.gz"
        url_r2 = f"{base_url}/{subdir}/{sra_id}_2.fastq.gz"

        # Output file paths
        output_r1 = os.path.join(output_dir, f"{sra_id}_1.fastq.gz")
        output_r2 = os.path.join(output_dir, f"{sra_id}_2.fastq.gz")

        print(f"🔄 Processing {sra_id} ...")

        # Download Read 1
        download_with_resume(url_r1, output_r1)

        # Download Read 2
        download_with_resume(url_r2, output_r2)

# Prompt user for inputs with tab completion
accession_file = input("📂 Enter the path to the text file with SRA accessions (Tab to autocomplete): ").strip()
output_dir = input("📁 Enter the output directory for FASTQ files (Tab to autocomplete): ").strip()

# Run the downloader
download_from_ena(accession_file, output_dir)