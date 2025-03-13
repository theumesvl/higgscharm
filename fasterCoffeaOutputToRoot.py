import os
import coffea.util as cu
import multiprocessing

def process_coffea_file(coffea_file):
    """Loads a Coffea file and converts it to ROOT."""
    try:
        print(f"Processing {coffea_file}...")
        out = cu.load(coffea_file)
        # Call your ROOT conversion function here
        # write_root(out, save_path, args)
        print(f"Finished processing {coffea_file}")
    except Exception as e:
        print(f"Error processing {coffea_file}: {e}")

def convert_coffea_to_root_parallel(output_dir, num_workers=4):
    """Converts all Coffea files in the given directory using multiprocessing."""
    coffea_files = []

    for sample in os.listdir(output_dir):
        sample_path = os.path.join(output_dir, sample)

        if not os.path.isdir(sample_path):
            continue

        for file in os.listdir(sample_path):
            coffea_file = os.path.join(sample_path, file)

            if os.path.isfile(coffea_file) and file.endswith(".coffea"):
                coffea_files.append(coffea_file)

    # Use multiprocessing to process files in parallel
    with multiprocessing.Pool(processes=num_workers) as pool:
        pool.map(process_coffea_file, coffea_files)

# Example usage
if __name__ == "__main__":
    processor = "WWtoMuEle"
    year = "2022postEE"
    output_dir = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}"
    
    convert_coffea_to_root_parallel(output_dir, num_workers=8)

