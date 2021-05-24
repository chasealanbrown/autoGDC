import re
import logging
import toml
from os import path
from numpy import float64

# Logger for package
logging.getLogger().setLevel(logging.DEBUG)
LOG = logging.getLogger("autoGDC")

###############################################################################
# Load the config and settings files (static/no-logic part of config)
###############################################################################

this_path = path.dirname(path.abspath(__file__))

# Load the local system settings for data paths
with open(path.join(this_path, "..", "config.toml")) as f:
  conf = toml.loads(f.read())
  # Check if the local data directory needs to be set to the default
  if conf["default"]["data_dir"] == "":
    conf["default"]["data_dir"] = path.expanduser(path.join("~", ".autogdc", "data"))

# Load the GDC API key for access to controlled data
with open(path.join(this_path, "..", "GDC_API.key")) as f:
  if f.read() == "":
    gdc_api_key = None
  else:
    l = list(f.readlines())
    if len(l)>1:
      LOG.error("Ensure that the 'GDC_API.key' file\
                has only one line of text")
    gdc_api_key = l[0].strip()

# Load other GDC settings
# - This is typically configurable by devs
with open(path.join(this_path, "settings.toml")) as f:
  settings = toml.loads(f.read())

###############################################################################
# Config: Apply logic
###############################################################################

# Create directory structure for local data directory on different systems
for k in conf:
  data_dir = path.expanduser(conf[k]["data_dir"])
  # Directories and files in the "archive" directory
  arxiv_dir = path.join(data_dir, "archive")
  conf[k]["archive_dir"] = arxiv_dir
  conf[k]["newdata_dir"] = path.join(arxiv_dir, "new")
  conf[k]["rawdata_dir"] = path.join(arxiv_dir, "raw")
  conf[k]["mygene_dir"] = path.join(arxiv_dir, "mygene")
  conf[k]["metadata_dir"] = path.join(arxiv_dir, "metadata")
  conf[k]["metadb_path"] = path.join(arxiv_dir, "metadb.tsv.gz")
  conf[k]["metadb_index_path"] = path.join(arxiv_dir,
                                           "metadb_index.tsv.gz")
  # Create directories and metadata filepaths for all assay types
  conf[k]["assay_dir"] = {}
  conf[k]["metadata_filepath"] = {}
  conf[k]["db_file"] = {}
  for assay in settings["filetype_regex"]:
    conf[k]["assay_dir"][assay] = path.join(arxiv_dir, assay)
    conf[k]["metadata_filepath"][assay] = path.join(arxiv_dir,
                                                    assay,
                                                    "feature_metadata.tsv")
    # A dictionary containing the file paths of the databases for each assay.
    # (combined gdc files with the same assay type, indexed by assay type keys)
    conf[k]["db_file"][assay] = path.join(data_dir, f"{assay}.h5")

###############################################################################
# Settings: Apply logic
###############################################################################

# Several ways of reporting a none/null type in GDC exist in string form
settings["null_type_strings"] = ["not reported",
                                 "unknown",
                                 "--",
                                 "na",
                                 "n/a",
                                 "null",
                                 "none"]

settings["meth_dtypes"] = {"Composite Element REF": str,
               "Beta_value": float64,
               "Chromosome":str,
               "Start":int,
               "End":int,
               "Gene_Symbol":str,
               "Gene_Type":str,
               "Transcript_ID":str,
               "Position_to_TSS":str,
               "CGI_Coordinate":str,
               "Feature_Type":str}

settings["mirna_dtypes"] = {"miRNA_ID": str,
                "isoform_coords": str,
                "read_count": int,
                "reads_per_million_miRNA_mapped": float64,
                "cross-mapped": str,
                "miRNA_region": str}

# If you have a GDC API key for controlled access,
#   place it in the GDC_API.key file
settings["api_key"] = gdc_api_key

settings["read_params"] = {}
for assay, regex in settings["filetype_regex"].items():
  # Compile the regular expressions for the filetypes
  settings["filetype_regex"][assay] = re.compile(regex)
  # Parameters for reading different assay types
  params = {"value_name": "beta_value" if "DNAm" in assay
                                       else "value",
            "index_name": "loci" if "DNAm" in assay
                                 else "index",
            "dtypes": settings["meth_dtypes"] if "DNAm" in assay
                                              else (settings["mirna_dtypes"]
                                                    if assay in ["RNA_miRNA",
                                                                 "RNA_isoforms"]
                                                    else None),
            "header": None if assay in ["RNA_counts",
                                        "RNA_FPKM",
                                        "RNA_FPKM-UQ"]
                           else 0,
            "subset_cols": [0,1] if "DNAm" in assay
                                 else ([0,3] if assay == "RNA_isoforms"
                                             else([0,2] if assay == "RNA_miRNA"
                                                        else None))
            }
  settings["read_params"][assay] = params

# Store all of these settings in a dictionary
SETTINGS = {config_key: dict(settings, **vals)
            for config_key, vals in conf.items()}
