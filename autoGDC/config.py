import re
import os
import toml
from os import path
from numpy import float64
import diskcache
from wrenches.general import calling_module_name, pkg_logger, pkg_default_cache_directory

PKGNAME = calling_module_name()
LOG = pkg_logger(name=PKGNAME)#(name=None, level=logging.INFO, console_level=logging.DEBUG):


###############################################################################
# Load the config and settings files (static/no-logic part of config)
###############################################################################

this_path = path.dirname(path.abspath(__file__))

# Load the local system settings for data paths
with open(path.join(this_path, "..", "config.toml")) as f:
  conf = toml.loads(f.read())
  # Default and test directories have a default if empty string is given
  for conf_type in ["default", "test"]:
    # Check if the local data directory needs to be set to the default
    if conf[conf_type]["data_dir"] == "":
      data_dir = path.expanduser(path.join("~", f".{PKGNAME}", "data"))
      conf[conf_type]["data_dir"] = data_dir

    if not os.path.exists(conf[conf_type]["data_dir"]):
      os.makedirs(conf[conf_type]["data_dir"])

    # Check if the local cache directory needs to be set to the default
    if conf[conf_type]["cache_dir"] == "":
  #    cache_path = ["~", ".cache", PKGNAME]
  #    cache_path = path.expanduser(path.join(*cache_path))
      conf[conf_type]["cache_dir"] = pkg_default_cache_directory(name=PKGNAME)
    if not os.path.exists(conf[conf_type]["cache_dir"]):
        os.makedirs(conf[conf_type]["cache_dir"])

CACHE = diskcache.Cache(dir=conf["default"]["cache_dir"])

# Load the GDC API key for access to controlled data
with open(path.join(this_path, "..", "GDC_API.key")) as f:
  if f.read() == "":
    gdc_api_key = None
  else:
    l = list(f.readlines())
    if len(l)>1:
      LOG.error("ERROR! The GDC_API.key file is improperly formatted!\
                Ensure that the 'GDC_API.key' file\
                has only one line of text!")
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
  conf[k]["metadb_path"] = path.join(arxiv_dir, "metadb.h5")#.tsv.gz")
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

settings["metadata_dtypes"] = {"id":str,
                               "filename":str,
                               "md5":str,
                               "size":int,
                               "state":str,
                               "acl":str,
                               "project_id":str,
                               "category":str,
                               "data_type":str,
                               "primary_site":str,
                               "disease_type":str,
                               "case_id":str,
                               "submitter_id":str,
                               "tumor_descriptor":str,
                               "sample_id":str,
                               "sample_type":str,
                               "tissue_type":str,
                               "aliquot_id":str,
                               "gender":str,
                               "race":str,
                               "age_at_diagnosis":float,
                               "tumor_stage":str,
                               "age":float,
                               "cd4_count":float,
                               "risk_factor":str,
                               "bmi":float,
                               "diabetes_treatment_type":str,
                               "downloaded":bool,
                               "organized":bool}

nested_index_list = ["cases",
                     "portions",
                     "samples",
                     "analytes",
                     "aliquots",
                     "follow_ups",
                     "diagnoses"]
nested_index_dtypes = {f"nested_{k}_index":float for k in nested_index_list}
settings["metadata_dtypes"].update(nested_index_dtypes)

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
            "dtype": settings["meth_dtypes"] if "DNAm" in assay
                                              else (settings["mirna_dtypes"]
                                                    if assay in ["RNA_miRNA",
                                                                 "RNA_isoforms"]
                                                    else None),
            "header": None if assay in ["RNA_counts",
                                        "RNA_FPKM",
                                        "RNA_FPKM-UQ"]
                           else 0,
            "subset_col": [0,1] if "DNAm" in assay
                                 else ([0,3] if assay == "RNA_isoforms"
                                             else([0,2] if assay == "RNA_miRNA"
                                                        else None))
            }
  settings["read_params"][assay] = params

# Store all of these settings in a dictionary
SETTINGS = {config_key: dict(settings, **vals)
            for config_key, vals in conf.items()}
