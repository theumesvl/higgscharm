import os
import re

# define year and processor
year = "2022postEE"
processor = "WWtoMuEle"

# Define the base directory
base_dir = f"/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/condor/{processor}/{year}"
output_file = "xrootd_path_useds.txt"

# Regular expression to find XRootD paths
xrootd_pattern = re.compile(r"root://[a-zA-Z0-9\-.]+:\d+/store/[^\s]+")

# List to store found paths
xrootd_paths = []

# Iterate over each sample directory
for sample_dir in sorted(os.listdir(base_dir)):
    sample_path = os.path.join(base_dir, sample_dir)
    
    # Ensure it's a directory
    if not os.path.isdir(sample_path):
        continue

    # Find the .sh file
    for filename in os.listdir(sample_path):
        if filename.endswith(".sh"):
            sh_file_path = os.path.join(sample_path, filename)
            
            # Read the .sh file and extract XRootD paths
            with open(sh_file_path, "r") as f:
                content = f.read()
                matches = xrootd_pattern.findall(content)
                xrootd_paths.extend(matches)

# Save extracted paths to a file
with open(output_file, "w") as f:
    for path in sorted(set(xrootd_paths)):  # Remove duplicates and sort
        f.write(path + "\n")

print(f"Extracted {len(set(xrootd_paths))} unique XRootD paths and saved to {output_file}")

