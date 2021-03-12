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

# Load the local system settings for data paths
with open(path.join("..", "config.toml")) as f:
  conf = toml.loads(f.read())
  # Check if the local data directory needs to be set to the default
  if conf["default"]["data_dir"] == "":
    conf["default"]["data_dir"] = path.join("~", ".autogdc", "data")

# Load the GDC API key for access to controlled data
with open(path.join("..", "GDC_API.key")) as f:
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
with open("settings.toml") as f:
  settings = toml.loads(f.read())

###############################################################################
# Config: Apply logic
###############################################################################

# Create directory structure for local data directory on different systems
for k in conf:
  data_dir = conf[k]["data_dir"]
  # Directories and files in the "archive" directory
  arxiv_dir = path.join(data_dir, "archive")
  conf[k]["newdata_dir"] = path.join(arxiv_dir, "new")
  conf[k]["rawdata_dir"] = path.join(arxiv_dir, "raw")
  conf[k]["mygene_dir"] = path.join(arxiv_dir, "mygene")
  conf[k]["metadb_path"] = path.join(arxiv_dir, "metadb.tsv.gz")
  conf[k]["metadb_index_path"] = path.join(arxiv_dir,
                                                 "metadb_index.tsv.gz")
  # Create directories and metadata filepaths for all assay types
  conf[k]["assay_dirs"] = {}
  for assay in settings["filetype_regex"]:
    conf[k]["assay_dirs"][assay] = path.join(arxiv_dir, assay)
    conf[k][f"{assay}_metadata_filepath"] = path.join(arxiv_dir,
                                                      assay,
                                                      "feature_metadata.tsv")

###############################################################################
# Settings: Apply logic
###############################################################################

# If you have a GDC API key for controlled access,
#   place it in the GDC_API.key file
settings["api_key"] = gdc_api_key
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
