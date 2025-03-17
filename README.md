    # Cleaning Bases Using Bioinformatics Tools
    #### Video Demo:  https://youtu.be/ydWB4gvSWf8
    #### Description:
    This Python script automates the workflow of retrieving, quality checking, and processing sequencing data from the SRA (Sequence Read Archive). The pipeline includes searching for SRA datasets, downloading FASTQ files, running quality checks with FastQC, and trimming reads with Trimmomatic based on quality thresholds.

    Software Prerequistes:
    SRA Toolkit
    Trimmomatic
    FastQC

    Python Packages:
    biopython
    tqdm

    Usage:
    Modify the following paths in the script to match your system:
    Entrez.email: Set your email for NCBI Entrez queries.
    sra_toolkit_path: Path to the SRA Toolkit installation.
    Paths for output directories: output_dir, trimmed_dir, and trimmed_fastqc_dir.

    Run the program and enter a query or accession number when prompted. Do make sure that it is a paired end read.

    The script sets up a local port for viewing the FastQC .html reports in your browser, enabling detailed visualization of quality metrics.

    Configuration:
    FastQC:
    Default settings are used. Reports include per-base quality, sequence duplication levels, and GC content.

    Trimmomatic:
    SLIDINGWINDOW:4:30: Trims reads if the average quality within a 4-base sliding window drops below 30.
    MINLEN:36: Discards reads shorter than 36 bases.

    Functions:

    search_sra(query, max_result):
    The search_sra function allows users to search the SRA (Sequence Read Archive) database using a query term (e.g., a study name, experiment, or organism) or a specific SRA accession number. The function retrieves metadata for the matching sequencing runs by querying the SRA database and parsing the XML data returned. It then extracts and returns a list of SRA accession numbers for the sequencing runs that match the search criteria. This function is useful for automating the retrieval of sequencing data from SRA based on user-defined search terms or known accession numbers.
    Parameters:
    query: A string to search the SRA database (e.g., a term like "human genome" or a specific accession number).
    max_results: The maximum number of results to return from the search (default is 1).

    download_sra_data(accession)
    The download_sra_data function automates the process of downloading sequencing data from the SRA (Sequence Read Archive) using the SRA Toolkit. It uses the fastq-dump command to retrieve paired-end FASTQ files for a specified accession number. The function tracks the download progress, displaying the elapsed time until completion. If the download is interrupted, it handles the exception gracefully and informs the user.
    Parameters:
    accession: The SRA accession number to download data for.


    run_fastqc(file1, file2, output_dir=None)
    The run_fastqc function runs the FastQC tool on paired-end FASTQ files to generate quality reports. The function executes the FastQC command-line tool and saves the output reports either in the current directory or in a specified output directory. If any errors occur during the process, they are displayed to the user.
    Parameters:
    file1: Path to the first FASTQ file.
    file2: Path to the second FASTQ file.
    output_dir: (Optional) Directory to save the FastQC reports.


    find_paired_files(directory)
    The find_paired_files function identifies paired-end FASTQ files in a given directory by matching filenames that follow a specific pattern (e.g., sample_1.fastq.gz and sample_2.fastq.gz). It returns a list of tuples where each tuple contains the paths to the paired files. If a complete pair is found, it is included in the returned list.
    Parameters:
    directory: Path to the directory containing FASTQ files.
    Returns:
    A list of tuples, where each tuple contains the paths to a pair of FASTQ files.


    unzip_and_check_quality(zip_file)
    The unzip_and_check_quality function unzips a FastQC report, then checks the per-base sequence quality to determine if trimming is necessary. It calls the check_per_base_quality function to identify if any base quality is below a defined threshold (e.g., 30). The function returns True if trimming is needed and False otherwise.
    Parameters:
    zip_file: Path to the FastQC report ZIP file.
    Returns:
    True if trimming is needed based on per-base quality; False otherwise.

    check_per_base_quality(fastqc_data_path)
    The check_per_base_quality function checks the per-base sequence quality from a FastQC data file (typically fastqc_data.txt). It parses the file to evaluate the quality score at each base position. If any base has a quality score below 30, it returns True, indicating that trimming is needed. If no low-quality bases are found, it returns False.
    Parameters:
    fastqc_data_path: Path to the FastQC data file containing sequence quality information.
    Returns:
    True if any base quality is below 30; False otherwise.


    run_trimmomatic(file1, file2, output_dir)
    The run_trimmomatic function uses the Trimmomatic tool to trim paired-end FASTQ files based on quality thresholds. It applies a sliding window approach to trim reads if the average quality within a 4-base window falls below 30, and it discards reads shorter than 36 bases. The trimmed files are saved in the specified output directory.
    Parameters:
    file1: Path to the first FASTQ file.
    file2: Path to the second FASTQ file.
    output_dir: Directory to save the trimmed FASTQ files.


    run_fastqc_on_trimmed_files(trimmed_files, output_dir)
    The run_fastqc_on_trimmed_files function runs FastQC on a list of trimmed FASTQ files to generate quality reports for each pair. It ensures that the output directory exists and saves the reports there. This function is useful for assessing the quality of the trimmed reads.
    Parameters:
    trimmed_files: A list of tuples containing paths to trimmed FASTQ files.
    output_dir: Directory to save the FastQC reports.








