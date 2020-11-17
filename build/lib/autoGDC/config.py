from os import path, listdir

# Infinium feature metadata
#DNAm_27k_url = "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanmethylation27/productsupportfiles/illumina_humanmethylation27_content.xlsx"
# Wrenlab feature metadata
DNAm_27k_url = "https://www.dropbox.com/s/o36utk95fl1euy1/HumanMethylation27_feature_metadata.tsv?dl=0"
#DNAm_450k_url = "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv"
DNAm_450k_url = "https://www.dropbox.com/s/wyck0lsa6941utw/HumanMethylation450_feature_metadata.tsv?dl=0"

# Absolute path for this directory
this_dir = path.dirname(path.realpath(__file__))

# Binary program released by GDC for downloading files
gdc_client_path = path.join(this_dir, "bin", "gdc-client")

# Different system configurations for data directory locations
default_data_dir = path.join(this_dir, "data", "archive")
wrenlab_data_dir = path.join("data", "autogdc", "archive")

# Different system configurations for Mygene information database
default_mg_dir = path.join(this_dir, "data", "mygene")
wrenlab_mg_dir = path.join("data", "databases", "mygene")

# Load the GDC API key for access to controlled data
with open("GDC_API.key") as f:
  s = f.read()
  if s == "":
    gdc_api_key = None
  else:
    gdc_api_key = s

###############################################################################
# These values can, and will be changed often (between each new study)
#   Also:
#     What is written here will only be used if no parameters are provided
#       the `DataSet` object.
#     Otherwise, the parameters provided there overwrite these values.
###############################################################################
#dataset_settings = {
#    # This is a parameter that defines the data that will be gathered
#    #   via downloading and/or selecting from the data archive stored locally.
#    #   It has this strange json-like form because it is used directly
#    #   in the REST API for GDC via the gdc-client binary
#    "filt":
#            # Default to Glioma RNA and DNA
#            {"op":'and',
#              "content":
#                [{"op":"IN",
#                  "content":
#                      {"field": 'cases.disease_type',
#                       "value": ['Gliomas',
#                                 'Neuroblastoma']}},
#                 {"op":"IN",
#                  "content":
#                      {"field": 'files.analysis.workflow_type',
#                       "value": ["HTSeq - FPKM",
#                                 "Liftover"]}}]},
#
#    # Upper limit for the number of file to load/download
#    #   We default to a non-huge number here to save people time when testing
#    "size": 10**2,
#
#    # The different contrasts to be considered in differentially expressed
#    # gene analysis
#    "contrasts": ["sample_type"],
#
#    # Whether or not to pair data on RNA and DNAm axis
#    #   TODO: This needs to be improved
#    "paired_assay": False
#
#}



###############################################################################
# The following values will not be changed often (only specific cases)
###############################################################################
gdc_settings = {
    # If you have a GDC API key for controlled access,
    #   place it in the GDC_API.key file
    "api_key": gdc_api_key,

    # These are the sample metadata features that we will pull
    #   from the GDC REST API for determining what data we have on disk,
    #   and for organizing the data
    "fields":
            ["id",
             "file_id",
             "file_name",
             "md5sum",
             "file_size",
             "data_type",
             "data_category",
             "access",
             "submitter_id",
             "primary_site",
             "archive.state",
             "files.data_category",
             "analysis.workflow_type",
             "cases.submitter_id",
             "cases.case_id",
             "cases.primary_site",
             "cases.disease_type",
             "cases.demographic.race",
             "cases.demographic.gender",
             "cases.diagnoses.tumor_stage",
             "cases.diagnoses.age_at_diagnosis",
             "cases.samples.tumor_descriptor",
             "cases.samples.tissue_type",
             "cases.samples.sample_type",
             "cases.samples.sample_id",
             "cases.samples.portions.analytes.aliquots.aliquot_id",
             "cases.project.project_id"
             "cases.follow_ups.comorbidity",
             "cases.follow_ups.diabetes_treatment_type",
             "cases.follow_ups.bmi",
             "cases.follow_ups.risk_factor",
             "cases.follow_ups.risk_factor_treatment",
             "cases.follow_ups.cd4_count"
             "diagnoses.vital_status"],

	# DNA methylation feature dtypes
    #   (needed for faster parsing of files by Cython)
    "meth_dtypes":
            {"Composite Element REF": str,
             "Beta_value": np.float64,
             "Chromosome": str,
             "Start": int,
             "End": int,
             "Gene_Symbol": str,
             "Gene_Type": str,
             "Transcript_ID": str,
             "Position_to_TSS": str,
             "CGI_Coordinate": str,
             "Feature_Type": str},

    "mirna_dtypes":
            {"miRNA_ID": str,
             "isoform_coords": str,
             "read_count": int,
             "reads_per_million_miRNA_mapped": np.float64,
             "cross-mapped": str,
             "miRNA_region": str},

    "filetype_regexs":
            {"DNAm_450":"HumanMethylation450",
             "DNAm_27":"HumanMethylation27",
             "RNA_FPKM":"FPKM.txt",
             "RNA_FPKM-UQ":"FPKM-UQ",
             # Files can have either '-' or '.' in name
             "RNA_counts":"htseq[\.\-]counts.gz",
             "RNA_miRNA":"mirnas.quantification",
             "RNA_isoforms":"isoforms.quantification"},

    "null_type_strings":["not reported", "--", "na", "n/a", "null", "none"],

    "gdc_client_path": gdc_client_path,
    "DNAm_450k_url": DNAm_450k_url,
    "DNAm_27k_url": DNAm_27k_url,
}

###############################################################################
# These values keep track of computer configurations
#   i.e. if you want to store the downloaded data at another location on disk,
#   you can add your own configuration key here
#   (shown here as 'default' or 'wrenlab')
#   This will likely only be changed once, during initial setup
#   In order to use your configuration key, set the `config_key` parameter
#   in the DataSet object
###############################################################################
local_settings = {
                  # Default configuration uses the `data` directory within
                  #   this git repo, as any other assumptions of file storage
                  #   may not work for the user
                  "default":
                    {"data_dir": default_data_dir,
                     "mygene_dir": default_mg_dir,
                     "assay_dirs":
                     {assay_dir: path.join(default_data_dir, assay_dir)
                      for assay_dir in gdc_settings["filetype_regexs"].keys()},
                     "DNAm_450k_metadata_filepath":
                     path.join(default_data_dir, "DNAm_450", "feature_metadata.tsv"),
                     "DNAm_27k_metadata_filepath":
                     path.join(default_data_dir, "DNAm_27", "feature_metadata.tsv"),
                     # Whether or not to keep raw downloaded files in archive
                     "keep_raw": False},

                  # Wren lab specific configuration
                  #   (All our computers have a `/data/` directory
                  #    for storage of larger files)
                  "wrenlab":
                    {"data_dir": wrenlab_data_dir,
                     "mygene_dir": wrenlab_mg_dir,
                     "assay_dirs":
                     {assay_dir: path.join(wrenlab_data_dir, assay_dir)
                      for assay_dir in gdc_settings["filetype_regexs"].keys()},
                     "DNAm_450k_metadata_filepath":
                     path.join(wrenlab_data_dir, "DNAm_450", "feature_metadata.tsv"),
                     "DNAm_27k_metadata_filepath":
                     path.join(wrenlab_data_dir, "DNAm_27", "feature_metadata.tsv"),
                     "keep_raw": True}
                 }


# Store all of these settings in a dictionary, keyed by the local_settings key
#   (Such that we can obtain the correct settings
#       given a config_key, such as 'wrenlab', across all classes)
SETTINGS = {config_key: gdc_settings\
                            .update(lclvals)#\
#                            .update(dataset_settings)\
            for config_key, lclvals in local_settings.items()}

