import os
import re

# define processor and year to construct paths
year = "2022postEE"
processor = "WWtoMuEle"

# Define paths
sample_dir = f"/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/condor/{processor}/{year}"  # Condor logs directory
output_dir = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}"  # Output directory with ROOT files
log_dir = f"/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/condor/logs/{processor}/{year}"  # Condor log directory
expected_samples_file = "samples_list.txt"  # File with expected samples (if manually created)
found_samples_file = "samples.txt"  # File to store detected sample names
missing_samples_file = "missing_samples.txt"  # File to store missing samples

# Switch between directory-based and file-based sample retrieval
use_directory = True  # Set to False if using manually created text file

def get_samples_from_directory(sample_dir):
    """Retrieve sample names from Condor logs directory, extract base names and IDs, and save them to a file."""
    samples = []
    
    for dirname in sorted(os.listdir(sample_dir)):  # Get subdirectories
        parts = dirname.split("_")  # Split name by underscores
        if parts[-1].isdigit():  # Ensure last part is a number
            sample_name = "_".join(parts[:-1])  # Reconstruct base sample name
            sample_id = parts[-1]  # Extract ID
            samples.append((sample_name, sample_id))

    # Save found sample names to a file
    with open(found_samples_file, "w") as f:
        for sample_name, sample_id in samples:
            f.write(f"{sample_name}_{sample_id}\n")

    print(f"Saved {len(samples)} sample names to {found_samples_file}")
    return samples  # Return as list of tuples


def read_sample_list(file_path):
    """Extracts sample names from the text file."""
    samples = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip().split()[-1]  # Extract last word (sample name)
            parts = line.split("_")  # Split by underscores
            if parts[-1].isdigit():  # Check if last part is a number
                sample_name = parts[-2]  # Base sample name
                sample_id = parts[-1]  # Sample number
                samples.append((sample_name, sample_id))
    return samples

def get_latest_error(log_dir, sample_name, sample_id):
    """Find the latest .err file for a missing sample and extract relevant error information."""
    sample_log_dir = os.path.join(log_dir, f"{sample_name}_{sample_id}")
    if not os.path.exists(sample_log_dir):
        return "No .err file", ""

    err_files = [f for f in os.listdir(sample_log_dir) if f.endswith(".err")]
    if not err_files:
        return "No .err file", ""

    latest_err_file = sorted(err_files, key=lambda x: int(x.split(".")[-3]))[-1]
    err_path = os.path.join(sample_log_dir, latest_err_file)
    
    with open(err_path, "r") as f:
        err_content = f.read()
    
    # Extract the first relevant error message
    error_matches = re.findall(r"([a-zA-Z]+Error): (.+)", err_content)
    error_msg = error_matches[-1][0] + ": " + error_matches[-1][1] if error_matches else "Unknown error"

    # Extract XRootD storage site endpoint
    xrootd_matches = re.findall(r"root://[a-zA-Z0-9\-.]+(?:[:]\d+)?", err_content)
    print(xrootd_matches)
    xrootd_site = xrootd_matches[0] if xrootd_matches else ""

    return error_msg, xrootd_site

def check_missing_samples(base_dir, log_dir, samples, missing_samples_file):
    """Checks which samples exist in the output directory and extracts errors for missing ones."""
    missing_samples = []

    for sample_name, sample_id in samples:
        sample_dir = os.path.join(base_dir, sample_name)
        sample_file = f"{sample_name}_{sample_id}.root"
        file_path = os.path.join(sample_dir, sample_file)

        if os.path.exists(file_path):
            print(f"\033[1m{sample_file} is present\033[0m")  # Bold output
        else:
            print(f"{sample_file} is MISSING")
            error_msg, xrootd_site = get_latest_error(log_dir, sample_name, sample_id)
            missing_samples.append((sample_file.strip().split(".")[0], xrootd_site, error_msg))  # strip of root extension of sample file

    # Save missing samples to file
    with open(missing_samples_file, "w") as f:
        f.write(f"{'Sample':<80}{'XRootD Site':<40}{'Error Message'}\n")
        f.write("="*140 + "\n")
        for sample, site, error in missing_samples:
            f.write(f"{sample:<80}{site:<40}{error}\n")

    print(f"Saved missing samples list with errors to {missing_samples_file}")

# Get samples based on the selected method
if use_directory:
    detected_samples = get_samples_from_directory(sample_dir)
else:
    detected_samples = read_sample_list(expected_samples_file)

# Check which samples are missing
check_missing_samples(output_dir, log_dir, detected_samples, missing_samples_file)    
