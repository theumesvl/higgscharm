import os
import time

def remove_old_root_files(base_dir, days=1):
    """Removes .root files in subdirectories of base_dir if they are older than 'days'."""
    cutoff_time = time.time() - days * 86400  # Convert days to seconds

    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)

        if not os.path.isdir(subdir_path):  # Skip if it's not a directory
            continue

        for file in os.listdir(subdir_path):
            file_path = os.path.join(subdir_path, file)

            if file.endswith(".root") and os.path.isfile(file_path):
                file_creation_time = os.path.getctime(file_path)  # Get creation time

                if file_creation_time < cutoff_time:
                    print(f"Deleting: {file_path}")
                    os.remove(file_path)

# Example usage
if __name__ == "__main__":
    processor = "WWtoMuEle"
    year = "2022postEE"
    directory = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}" # Change this to your target directory
    remove_old_root_files(directory, days=1)

