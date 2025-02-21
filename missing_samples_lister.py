import os
import re

# define processor and year to construct paths
year = "2022postEE"
processor = "WWtoMuEle"

# Switch between directory-based and file-based sample retrieval
use_directory = True  # Set to False if using manually created text file

# Define paths
base_dir = "/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm"
sample_dir = f"/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/condor/{processor}/{year}"  # Condor logs directory
output_dir = f"/eos/user/t/tvanlaer/higgscharm/outputs/{processor}/{year}"  # Output directory with ROOT files
log_dir = f"/afs/cern.ch/user/t/tvanlaer/Hc/higgscharm/condor/logs/{processor}/{year}"  # Condor log directory
expected_samples_file = "samples_list.txt"  # File with expected samples (if manually created)
found_samples_file = "samples.txt"  # File to store detected sample names
missing_samples_file = "missing_samples.txt"  # File to store missing samples
storage_site_report_file = "storage_sites_report.txt"

# dictionary matching storage site names to xrootd endpoints: xrootd name should start with "root://" and go till the next "//" of the file paths (not including the "//")
xrootd_to_site = {
    "root://cmseos.fnal.gov": "T3_US_FNALLPC",
    "root://cmsdcadisk.fnal.gov": "T1_US_FNAL_Disk",
    "root://xrootd-redir1-vanderbilt.sites.opensciencegrid.org:1094": "T2_US_Vanderbilt",
    "root://eos.cms.rcac.purdue.edu": "T2_US_Purdue",
    "root://xrootd-local.unl.edu:1094": "T2_US_Nebraska",
    "root://dcache-cms-xrootd.desy.de:1094": "T2_DE_DESY",
    "root://maite.iihe.ac.be:1095": "T2_BE_IIHE",
    "root://eoscms.cern.ch": "T2_CH_CERN",
    "root://cmsdcache-kit-disk.gridka.de:1094": "T1_DE_KIT_Disk",
    "root://grid-cms-xrootd.physik.rwth-aachen.de:1094": "T2_DE_RWTH",
    "root://ingrid-se02.cism.ucl.ac.be:1094": "T2_BE_UCL",
    "root://rdr.echo.stfc.ac.uk": "T1_UK_RAL_Disk",
    "root://lcgsedr01.jinr.ru:1094": "T2_RU_JINR",
    "root://gfe02.grid.hep.ph.ic.ac.uk:1094": "T2_UK_London_IC",
    "root://ccxrdcms.in2p3.fr:1094": "T1_FR_CCIN2P3_Disk",
    "root://t3se01.psi.ch:1094": "T3_CH_PSI",
    "root://gaexrdoor.ciemat.es:1094": "T2_ES_CIEMAT",
    "root://xrootd-cms.infn.it:1194": "T2_IT_Pisa", # not sure
    "root://k8s-redir.ultralight.org:1094": "T2_US_Caltech", #not sure
    "root://se01.indiacms.res.in": "T2_IN_TIFR",
    "root://cmsio2.rc.ufl.edu:1094": "T2_US_Florida",
    "root://t2dsk0011.cmsaf.mit.edu:1094": "T2_US_MIT",
    "": "T2_FR_GRIF",
    "root://ce-03.recas.ba.infn.it": "T2_IT_Bari", # not sure
    "root://ruhex-osgce.rutgers.edu": "T3_US_Rutgers",
    "root://se.cis.gov.pl:1094": "T2_PL_Swierk",
    "root://cmsrm-cream01.roma1.infn.it": "T2_IT_Rome", # not sure
    "root://lyogrid07.in2p3.fr": "T3_FR_IPNL", # not sure
    "root://deepthought.crc.nd.edu": "T3_US_NotreDame", # not sure
    "root://ce04-lcg.cr.cnaf.infn.it": "T1_IT_CNAF_Disk", # not sure
    "root://osg-se.sprace.org.br:1094": "T2_BR_SPRACE",
    "root://storage01.lcg.cscs.ch:1096": "T2_CH_CSCS",
    "root://cmsxrootd.hep.wisc.edu:1094": "T2_US_Wisconsin",
    "root://cceos.ihep.ac.cn:1094": "T2_CN_Beijing",
    "root://redirector.t2.ucsd.edu:1095": "T2_US_UCSD",
}


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

    return missing_samples

def generate_storage_site_report(missing_samples, xrootd_to_site, storage_site_report_file):
    """Generate a report of all unique storage sites encountered that cause trouble."""
    storage_sites = set()

    for _, xrootd_site, _ in missing_samples:
        if xrootd_site:  # Ensure the site is not empty
            site_name = xrootd_to_site.get(xrootd_site, "Unknown site")  # Map XRootD to site name
            storage_sites.add(site_name)

    # Save storage site report
    with open(storage_site_report_file, "w") as f:
        f.write("Storage Site Name\n")
        f.write("=================\n")
        for site in sorted(storage_sites):
            f.write(f"{site}\n")

    print(f"Saved storage site report to {storage_site_report_file}")

# Get samples based on the selected method
if use_directory:
    detected_samples = get_samples_from_directory(sample_dir)
else:
    detected_samples = read_sample_list(expected_samples_file)

# Check which samples are missing
missing_samples = check_missing_samples(output_dir, log_dir, detected_samples, missing_samples_file)    

# Generate storage site report of troublesome storage sites
generate_storage_site_report(missing_samples, xrootd_to_site, storage_site_report_file)
