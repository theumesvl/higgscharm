import os
import pickle
import uproot
import coffea.util as cu

# Set this flag to True to overwrite existing ROOT files, False to skip them
OVERWRITE_EXISTING = False

def convert_coffea_to_root(output_dir):
    for sample in os.listdir(output_dir):
        sample_path = os.path.join(output_dir, sample)
        if not os.path.isdir(sample_path):
            continue
        
        for file in os.listdir(sample_path):
            coffea_file = os.path.join(sample_path, file)
            if not os.path.isfile(coffea_file) or not file.endswith(".coffea"):
                continue
            
            # Define output ROOT file path
            root_file_path = coffea_file.replace(".coffea", ".root")

            # Check if the ROOT file already exists
            if not OVERWRITE_EXISTING and os.path.exists(root_file_path):
                print(f"Skipping {coffea_file} (ROOT file already exists)")
                continue  # Skip processing if file exists

            # Load coffea file
            print(f"Processing {coffea_file}...")
            out = cu.load(coffea_file)
            
            # Save histograms to ROOT
            with uproot.recreate(root_file_path) as f:
                for hist_name, histogram in out["histograms"].items():
                    for category in histogram.axes["category"]:
                        category_histogram = histogram[{"category": category}]
                        variables = [v for v in category_histogram.axes.name if v != "variation"]
                        
                        for variable in variables:
                            for syst_var in category_histogram.axes["variation"]:
                                variation_histogram = category_histogram[{"variation": syst_var}]
                                f[f"{category}_{variable}_{syst_var}"] = variation_histogram.project(variable)
            
            print(f"Converted {coffea_file} -> {root_file_path}")

if __name__ == "__main__":
    processor = "WWtoMuEle"
    year = "2022postEE"
    output_dir = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}"
    convert_coffea_to_root(output_dir)
