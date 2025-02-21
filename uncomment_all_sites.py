import re

# Path to make_filesets.py
make_filesets_path = "/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/analysis/filesets/make_filesets.py" 

def uncomment_storage_sites(file_path):
    """Uncomments all storage sites in make_filesets.py."""
    with open(file_path, "r") as f:
        lines = f.readlines()

    # Regex to detect commented-out storage sites
    uncommented_lines = []
    for line in lines:
        # Look for lines starting with # followed by "root://" (or similar storage site formats)
        if re.match(r'^\s*#\s*(\"T0|\"T1|\"T2|\"T3)', line):
            uncommented_lines.append(" " * 8 + line.lstrip("# ").lstrip())  # Ensure 8-space indentation
        else:
            uncommented_lines.append(line)

    # Save back the updated file
    with open(file_path, "w") as f:
        f.writelines(uncommented_lines)

    print(f"Uncommented all storage sites in {file_path}")

# Run the function
uncomment_storage_sites(make_filesets_path)

