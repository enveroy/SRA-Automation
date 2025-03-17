from Bio import Entrez
import subprocess
import os
import xml.etree.ElementTree as ET
import time
import zipfile
from tqdm import tqdm
import re


Entrez.email = "alanjoshy101@gmail.com"  # NCBI requires an email
Entrez.api_key = "e377eb28f7818aa92fb7eb66ef9960467c08"# Optional for higher access limits
sra_toolkit_path = "/workspaces/162095528/sratoolkit.3.1.1-ubuntu64/bin"

def search_sra(query, max_results=1):
    """
    To search SRA database on a specified number of passes query or specific Accession Number
    """
    handle = Entrez.esearch(db="sra", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    b= record['IdList']
    sra_acc=[]

    if b:
        # To read the data in XML
        handle=Entrez.efetch(db="sra",id=b, rettype="xml")
        xml_data=handle.read()
        handle.close()

        root=ET.fromstring(xml_data)# to parse the XML data

        for run in root.findall(".//RUN"): sra_acc.append(run.attrib['accession'])

        return sra_acc

def download_sra_data(accession):
    cmd = [
        os.path.join(sra_toolkit_path, "fastq-dump"),
        "--split-files",
        "--gzip",
        accession
    ]

    start_time = time.time()
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        while True:
            # Check if the process has finished
            if process.poll() is not None:
                break

            # Calculate the elapsed time in seconds
            elapsed_time = time.time() - start_time
            # Display the timer in the format "Downloading... 00:00:00"
            print(f"\rDownloading... {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}", end="")
            time.sleep(1)

        stdout, stderr = process.communicate()# Ensure process completes

        if process.returncode != 0:
            print(f"\nError downloading data for accession {accession}: {stderr.decode()}")
        else:
            print(f"\nDownload completed for accession: {accession}")

    except KeyboardInterrupt:
        process.terminate()
        print("Download interrupted.")

def run_fastqc(file1, file2, output_dir=None):
    """
    Run FastQC on paired-end FASTQ files.
    """
    cmd = ["fastqc", file1, file2]
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        cmd.extend(["-o", output_dir])

    try:
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print(f"FastQC completed for {file1} and {file2}. Reports saved to {output_dir or 'current directory'}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running FastQC: {e.stderr.decode()}")

def find_paired_files(directory):
    """
    Find paired-end FASTQ files in the specified directory.

    Parameters:
        directory (str): Path to the directory containing FASTQ files.

    """
    # Regex pattern to match filenames in the format: prefix_1.fastq.gz or prefix_2.fastq.gz
    pattern = re.compile(r"(.*)_([12])\.fastq\.gz$")

    # List all files in the directory
    files = os.listdir(directory)

    # Dictionary to store paired files with the prefix as the key
    paired_files = {}

    # Iterate through each file in the directory
    for file in files:
        # Check if the file matches the pattern
        match = pattern.match(file)
        if match:
            # Extract the prefix and read number from the matched filename
            prefix, read_num = match.groups()  # prefix is the sample name, read_num is "1" or "2"

            # Initialize a new entry for this prefix if not already present
            if prefix not in paired_files:
                paired_files[prefix] = [None, None]

            # Place the file in the correct position (0 for read 1, 1 for read 2)
            paired_files[prefix][int(read_num) - 1] = os.path.join(directory, file)

    # Collect only complete pairs into the final list
    paired_file_list = []
    for pair in paired_files.values():
        if pair[0] and pair[1]:  # Ensure both read 1 and read 2 files are present
            paired_file_list.append((pair[0], pair[1]))

    return paired_file_list

def unzip_and_check_quality(zip_file):
    """
    Unzip a FastQC report, check per-base sequence quality, and return True if trimming is needed.
    """
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        extract_dir = zip_file.replace(".zip", "")
        zip_ref.extractall(extract_dir)

    # Locate the `fastqc_data.txt` file
    nested_dir = os.path.join(extract_dir, os.path.basename(extract_dir))
    fastqc_data_path = os.path.join(nested_dir, "fastqc_data.txt")

    if os.path.exists(fastqc_data_path):
        return check_per_base_quality(fastqc_data_path)
    else:
        print(f"fastqc_data.txt not found in {nested_dir}")
        return False

def check_per_base_quality(fastqc_data_path):
    """
    Check the per-base sequence quality from the FastQC data file.
    Return True if any base quality is below 30.
    """
    with open(fastqc_data_path, 'r') as file:
        lines = file.readlines()

    in_quality_section = False
    low_quality_detected = False

    for line in lines:
        if line.startswith(">>Per base sequence quality"):
            in_quality_section = True
        elif in_quality_section:
            if line.startswith(">>END_MODULE"):
                break
            if not line.startswith("#") and line.strip():
                # Parse the columns to get the position and quality score
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        position = parts[0]
                        quality = float(parts[1])
                        if quality < 30:
                            print(f"Low quality detected at position {position} with score {quality}")
                            low_quality_detected = True
                    except ValueError:
                        continue

    return low_quality_detected

def run_trimmomatic(file1, file2, output_dir):
    trimmed_file1 = os.path.join(output_dir, os.path.basename(file1).replace(".fastq.gz", "_trimmed.fastq.gz"))
    trimmed_file2 = os.path.join(output_dir, os.path.basename(file2).replace(".fastq.gz", "_trimmed.fastq.gz"))

    cmd = [
        "java", "-jar", "/usr/share/java/trimmomatic.jar", "PE", "-threads", "4",
        file1, file2,
        trimmed_file1, "/dev/null",
        trimmed_file2, "/dev/null",
        "SLIDINGWINDOW:4:30", "MINLEN:36"
    ]

    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if process.returncode != 0:
        print(f"Error running Trimmomatic: {process.stderr.decode()}")
    else:
        print(f"Trimming completed for {file1} and {file2}.")

def run_fastqc_on_trimmed_files(trimmed_files, output_dir):
    """
    Run FastQC on a specified list of trimmed FASTQ files and save the reports in the output directory.

    Parameters:
        trimmed_files (list): List of tuples containing paired trimmed FASTQ file paths.
        output_dir (str): Directory to save the FastQC reports.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    if not trimmed_files:
        print("No trimmed FASTQ files found.")
        return

    for file1, file2 in tqdm(trimmed_files, desc="Running FastQC on trimmed files"):
        run_fastqc(file1, file2, output_dir)

    print("FastQC analysis completed for all trimmed files.")


def main():
    query = input("Query: ")
    accessions = search_sra(query)
    downloaded_files = []
    for accession in accessions:
        print(f"Downloading {accession}")
        download_sra_data(accession)
        downloaded_files.extend([
            os.path.join("/workspaces/162095528", f"{accession}_1.fastq.gz"),
            os.path.join("/workspaces/162095528", f"{accession}_2.fastq.gz")
        ])
    print(accessions)

    # Process only the downloaded files
    print(f"\nDownloaded files: {downloaded_files}")
    paired_files = [(downloaded_files[i], downloaded_files[i + 1]) for i in range(0, len(downloaded_files), 2)]

    output_dir = "/workspaces/162095528/fastqR"
    trimmed_dir = "/workspaces/162095528/trimmed"
    trimmed_fastqc_dir = "/workspaces/162095528/trimmed_fastqc"

    # Step 1: Run FastQC on the downloaded files
    if not paired_files:
        print("No paired FASTQ files found.")
    else:
        for file1, file2 in tqdm(paired_files, desc="Running FastQC"):
            run_fastqc(file1, file2, output_dir)

    # Step 2: Unzip FastQC reports, check quality, and trim if needed
    print("\nUnzipping FastQC reports and checking per-base quality...")
    trimmed_files = []
    for file1, file2 in paired_files:
        prefix_name = os.path.basename(file1).split('_')[0]
        zip1 = os.path.join(output_dir, f"{prefix_name}_1_fastqc.zip")
        zip2 = os.path.join(output_dir, f"{prefix_name}_2_fastqc.zip")

        if os.path.exists(zip1) and os.path.exists(zip2):
            needs_trimming1 = unzip_and_check_quality(zip1)
            needs_trimming2 = unzip_and_check_quality(zip2)

            # Step 3: Trim reads if low-quality bases are detected
            if needs_trimming1 or needs_trimming2:
                print(f"Trimming reads for {prefix_name}...")
                run_trimmomatic(file1, file2, trimmed_dir)
                trimmed_files.append((
                    os.path.join(trimmed_dir, os.path.basename(file1).replace(".fastq.gz", "_trimmed.fastq.gz")),
                    os.path.join(trimmed_dir, os.path.basename(file2).replace(".fastq.gz", "_trimmed.fastq.gz"))
                ))
                print("Done")
            else:
                print(f"No trimming needed for {prefix_name}.")

    # Step 4: Run FastQC on trimmed files
    print("\nRunning FastQC on trimmed files...")
    run_fastqc_on_trimmed_files(trimmed_files, trimmed_fastqc_dir)

if __name__ == "__main__":
    main()
