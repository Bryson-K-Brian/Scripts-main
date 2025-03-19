import requests
import csv
import os
import readline
import glob

# Enable Tab Completion for File Input
def complete_filename(text, state):
    """Provides tab-completion for file paths."""
    return (glob.glob(text + '*') + [None])[state]

readline.set_completer(complete_filename)
readline.parse_and_bind("tab: complete")

# NCBI E-utilities URLs
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ELINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_sra_from_biosample(biosample_id):
    """Fetch SRA Run (SRR) IDs linked to a BioSample (SAMN)."""
    print(f"🔍 Searching NCBI SRA for BioSample: {biosample_id}")

    # Step 1: Search BioSample in NCBI
    params = {"db": "biosample", "term": biosample_id, "retmode": "json"}
    response = requests.get(NCBI_ESEARCH, params=params)
    response.raise_for_status()
    data = response.json()
    uid_list = data.get("esearchresult", {}).get("idlist", [])

    if not uid_list:
        print(f"❌ No records found for {biosample_id}")
        return []

    # Step 2: Link BioSample to SRA
    link_params = {"dbfrom": "biosample", "db": "sra", "id": ",".join(uid_list), "retmode": "json"}
    link_response = requests.get(NCBI_ELINK, params=link_params)
    link_response.raise_for_status()
    link_data = link_response.json()

    sra_ids = []
    for linkset in link_data.get("linksets", []):
        for link in linkset.get("linksetdbs", []):
            if link.get("dbto") == "sra":
                sra_ids.extend(link.get("links", []))

    if not sra_ids:
        print(f"❌ No linked SRA records found for {biosample_id}")
        return []

    # Step 3: Fetch SRR numbers
    fetch_params = {"db": "sra", "id": ",".join(sra_ids), "rettype": "runinfo", "retmode": "text"}
    fetch_response = requests.get(NCBI_EFETCH, params=fetch_params)
    fetch_response.raise_for_status()

    srr_list = []
    lines = fetch_response.text.split("\n")
    header = lines[0].split(",")  # CSV Header

    if "Run" in header:
        run_index = header.index("Run")
        for line in lines[1:]:
            fields = line.split(",")
            if len(fields) > run_index:
                srr_list.append(fields[run_index].strip())  # Remove spaces and blank lines

    return [srr for srr in srr_list if srr]  # Remove empty strings

def main():
    # Prompt user for input/output file with Tab Completion
    input_file = input("📂 Enter input file with BioSample IDs (Tab Complete Enabled): ").strip()

    if not os.path.exists(input_file):
        print(f"❌ File not found: {input_file}")
        return

    biosample_sra_file = "biosample_sra.csv"
    sra_only_file = "sra_only.csv"

    # Read BioSample IDs from the input file
    with open(input_file, "r") as file:
        biosample_ids = [line.strip() for line in file if line.strip()]

    results = []
    sra_only = []

    for biosample in biosample_ids:
        srrs = fetch_sra_from_biosample(biosample)
        if srrs:
            results.append([biosample, ", ".join(srrs)])
            sra_only.extend([[srr] for srr in srrs])  # Store SRRs separately

    # Write BioSample → SRA Mapping
    with open(biosample_sra_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["BioSample", "SRA Run IDs"])
        writer.writerows(results)
    print(f"✅ BioSample to SRA mapping saved: {biosample_sra_file}")

    # Write Only SRR IDs (Remove blanks)
    sra_only = [row for row in sra_only if row[0].strip()]  # Remove empty rows
    with open(sra_only_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["SRA Run ID"])
        writer.writerows(sra_only)
    print(f"✅ Only SRA Run IDs saved: {sra_only_file}")

if __name__ == "__main__":
    main()