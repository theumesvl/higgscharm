import os

# define processor and year to construct paths
year = "2022postEE"
processor = "WWtoMuEle"

# Define paths
sample_dir = f"/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/condor/{processor}/{year}"  # Condor logs directory
output_dir = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}"  # Output directory with ROOT files
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


def check_missing_samples(base_dir, samples):
    """Checks which samples exist in the output directory."""
    missing_samples = []

    for sample_name, sample_id in samples:
        sample_dir = os.path.join(base_dir, sample_name)
        sample_file = f"{sample_name}_{sample_id}.root"
        file_path = os.path.join(sample_dir, sample_file)

        if os.path.exists(file_path):
            print(f"\033[1m{sample_file} is present\033[0m")  # Bold output
        else:
            print(f"{sample_file} is MISSING")
            missing_samples.append(sample_file.strip().split(".")[0]) # strip of root extension

    # Save missing samples to file
    with open(missing_samples_file, "w") as f:
        for sample in missing_samples:
            f.write(sample + "\n")

    print(f"Missing samples saved to {missing_samples_file}")


# Get samples based on the selected method
if use_directory:
    detected_samples = get_samples_from_directory(sample_dir)
else:
    detected_samples = read_sample_list(expected_samples_file)

# Check which samples are missing
check_missing_samples(output_dir, detected_samples)
    
