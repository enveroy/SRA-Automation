import pytest
from unittest.mock import patch, MagicMock
from project import search_sra, check_per_base_quality, find_paired_files

def test_search_sra():
    # Mock the response from Entrez to avoid actual network calls
    mock_accession = ["SRR1234567"]
    with patch("project.Entrez.esearch") as mock_esearch, \
         patch("project.Entrez.read") as mock_read, \
         patch("project.Entrez.efetch") as mock_efetch:

        # Setup the mock to return specific values
        mock_esearch.return_value = MagicMock()
        mock_read.return_value = {"IdList": ["12345678"]}
        mock_efetch.return_value.read.return_value = (
            "<EXPERIMENT_PACKAGE_SET><RUN_SET><RUN accession='SRR1234567'></RUN></RUN_SET></EXPERIMENT_PACKAGE_SET>"
        )

        # Run the function
        result = search_sra("test_query")

        # Verify the result matches the mock
        assert result == mock_accession, f"Expected {mock_accession}, but got {result}"

def test_check_per_base_quality(tmp_path):
    # Create a mock FastQC data file with quality scores
    fastqc_data_path = tmp_path / "fastqc_data.txt"
    fastqc_data_content = (
        ">>Per base sequence quality\n"
        "#Base    Mean\n"
        "1        32.0\n"
        "2        28.0\n"
        "3        30.0\n"
        ">>END_MODULE\n"
    )

    # Write content to file
    fastqc_data_path.write_text(fastqc_data_content)

    # Call function and check if trimming is recommended (due to quality < 30 at position 2)
    result = check_per_base_quality(fastqc_data_path)

    # Expected: trimming is needed due to low quality at base 2
    assert result == True, "Expected low quality detected for base 2"

def test_find_paired_files(tmp_path):
    # Create mock files in the tmp_path directory
    file1 = tmp_path / "sample_1.fastq.gz"
    file2 = tmp_path / "sample_2.fastq.gz"
    file1.touch()  # Create an empty file
    file2.touch()  # Create an empty file

    # Run the function
    result = find_paired_files(tmp_path)

    # Expected: a list with a tuple of the paired files
    expected_result = [(str(file1), str(file2))]
    assert result == expected_result, f"Expected {expected_result}, but got {result}"
