import os
import coffea.util as cu
import multiprocessing
import uproot

# Set this flag to True to overwrite existing ROOT files, False to skip them
OVERWRITE_EXISTING = False

def process_coffea_file(coffea_file):
    """Loads a Coffea file and converts it to ROOT, if not already present or if overwriting."""
    try:
        root_file = coffea_file.replace(".coffea", ".root")

        if os.path.exists(root_file) and not OVERWRITE_EXISTING:
            print(f"Skipping {coffea_file}, ROOT file already exists.")
            return

        print(f"Processing {coffea_file}...")

        out = cu.load(coffea_file)

        # Convert and save as ROOT
        with uproot.recreate(root_file) as f:
            for hist_name, histogram in out["histograms"].items():
                for category in histogram.axes["category"]:
                    category_histogram = histogram[{"category": category}]
                    variables = [v for v in category_histogram.axes.name if v != "variation"]

                    for variable in variables:
                        for syst_var in category_histogram.axes["variation"]:
                            variation_histogram = category_histogram[{"variation": syst_var}]
                            f[f"{category}_{variable}_{syst_var}"] = variation_histogram.project(variable)

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

if __name__ == "__main__":
    processor = "WWtoMuEle"
    year = "2022postEE"
    output_dir = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}"

    convert_coffea_to_root_parallel(output_dir, num_workers=3)
    # Change OVERWRITE_EXISTING to True at the top if you want to overwrite existing ROOT files.

