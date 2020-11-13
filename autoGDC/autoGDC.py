# General dependencies
import os
import io
import re
import gc
import json
import shutil
import logging
import hashlib
import requests
import subprocess
import numpy as np
import pandas as pd

# Specific function imports
from os import path
from tqdm import tqdm
from joblib import Memory
from geode import chdir as characteristic_direction

# Local module imports
from . import df_utils
from .R.wrapper import pydeseq

# Gene information
import mygene
mg = mygene.MyGeneInfo()
mg.set_caching(cache_db="/data/databases/mygene/")

# Logger for autoGDC
logging.getLogger().setLevel(logging.INFO)
LOG = logging.getLogger(__name__)


meth_dtypes = {"Composite Element REF": str,
               "Beta_value": np.float64,
               "Chromosome":str,
               "Start":int,
               "End":int,
               "Gene_Symbol":str,
               "Gene_Type":str,
               "Transcript_ID":str,
               "Position_to_TSS":str,
               "CGI_Coordinate":str,
               "Feature_Type":str}

mirna_dtypes = {"miRNA_ID": str,
                "isoform_coords": str,
                "read_count": int,
                "reads_per_million_miRNA_mapped": np.float64,
                "cross-mapped": str,
                "miRNA_region": str}

# Absolute path for this file
this_path = path.dirname(os.path.realpath(__file__))

# Gene name information
#gene_info_path = path.join(this_path, "data", "mart_export.txt")
#gene_info = pd.read_csv(gene_info_path, sep = "\t")
#gene_id_name = gene_info[["Gene stable ID", "Gene name"]]\
#                        .drop_duplicates()\
#                        .set_index("Gene stable ID")\
#                        .iloc[:,0]
#gene_ensg_set = set(gene_id_name.index.tolist())

# Binary program released by GDC for downloading files
gdc_client_path = path.join(this_path, "bin", "gdc-client")


def init_cache():
  """
    Creates a cache to store files after large computations
  """
  usr_dir = path.expanduser("~")
  CACHE_DIRECTORY = path.join(usr_dir, ".cache", "autoGDC", "cache")
  os.makedirs(CACHE_DIRECTORY, exist_ok=True)
  memory = Memory(cachedir = CACHE_DIRECTORY, compress = True)
  return memory

memoize = init_cache().cache


def read_and_filter(filepath: str,
                    subset_features: list = None):
  """
  Summary:
    Reads a tsv file, formats any gene names, and renames it to GDC's `file_id`

  Arguments:
    filepath:
      Path to tsv downloaded from GDC

    subet_features:
      Subset of the features as a list
  """

  # Read the series
  #   (which was saved and compress during download phase)
  series = pd.read_csv(filepath,
                       sep = "\t",
                       index_col = 0,
                       header = None,
                       engine = "c",
                       squeeze = True)

  # Rename it to be the filename / file_id
  file_id = path.splitext(path.basename(filepath.rstrip(".gz")))[0]
  series.name = file_id

  series = series[~series.index.duplicated(keep="first")]

  # Select subset of features
  if subset_features:
    series = series.reindex(subset_features)
  return series


@memoize
def multifile_df(file_paths: list,
                 subset_features: list = None):
  """
  Summary:
    Reads a list of tsv files, formats any gene names,
      and renames it to GDC's `file_id`

  Arguments:
    filepaths:
      List of paths to tsv downloaded from GDC

    subet_features:
      Subset of the features as a list
  """

#  kwargs = {"subset_features": subset_features}

  # IO bound
  # Parallel imap with progress bar
  # series_iter = p_imap(partial(read_and_filter, **kwargs), file_paths)
  series_list = [read_and_filter(fp, subset_features)
                      for fp in tqdm(file_paths)]

  # Concatenate all series
  df = pd.DataFrame({s.name:s for s in series_list})

  df.columns.name = "file_id"
  return df.dropna(how = "all")


def contains_all_substrings(main_list, substr_list):
  """
  Summary:
    Determine if all substrings within a list of substrings are represented
      within a large list of strings

  Arguments:
    main_list:
      List of strings that will be tested to ensure that all substrings are
        represented within this list.

    substr_list:
      List of substrings that need to be represented to return `True`

  Example:
    substr_list = ["RNA", "DNA_methylation"]
    main_list = ["file_1_RNA", "file_1_DNA_methylation", "file_1_RNA_2"]
    constains_all_substrings(main_list, substr_list)
    # >>> True

    main_list = ["file_1_RNA", "file_1_RNA_2"]
    contains_all_substrings(main_list, substr_list)
    # >>> False
  """
  return all(True if any((substr in fullstr for fullstr in main_list))
             else False
             for substr in substr_list)


def subset_paired_assay(mdf):
  """
  Summary:
    Returns a subset of the metadata dataframe wherein all of the cases
      contain all of the biological assays considered for the study

  Arguments:
    mdf:
      Metadata dataframe from the autoGDC object
  """
  # Get all of the assay types:
  assay_types = mdf["workflow_type"].unique()

  # Only select cases for which have paired assays
  mdf_file_lists = mdf.groupby("case_id")["workflow_type"].apply(list)
  mdf_paired_case_bool = mdf_file_lists.apply(
                            lambda x: contains_all_substrings(x, assay_types))
  paired_cases = mdf_paired_case_bool[mdf_paired_case_bool].index

  mdf = mdf[mdf["case_id"].isin(paired_cases)]
  return mdf


class Dataset(object):
  """
  Summary:
    Class to automatically download data from GDC
      and store it in a more managable set of dataframes

  Arguments:
    Many of these arguments and the acceptable values can be seen at:
      https://docs.gdc.cancer.gov/Data_Dictionary/viewer/

    filt:
      A filter dictionary of data to search on GDC.

    fields:
      Fields to include in metadata

    size:
      Maximum number of samples requested

    contrasts:
      The label for which a model will be trained and
        differential expression / methylation will be reported
        i.e. ["age", "tumor_stage"]

    paired_assay:
      Force the data to contain only cases
        wherein all assay types are represented

    data_dir:
      full path to the directory containing the GDC downloaded data

    methyl_loci_subset:
      List of features that will be extracted from the overall GDC dataset

  Attributes:
    dataframe:
      A mulit-indexed dataframe containing
        both RNA and methylation data for each gene

    metadata:
      Defines what is the stage with the highest expected reward.
      Values: Initially random

    ddx:
      Differential expression / methylation results

  Example Usage:
    from autoGDC import autoGDC

    # This is the type of format that the gdc-client accepts
    filt = {"op":'and',
            "content":
              [{"op":"IN",
                "content":
                    {"field": 'cases.disease_type',
                    "value": ['Gliomas',
                              'Neuroblastoma']}},
                {"op":"IN",
                "content":
                    {"field": 'files.analysis.workflow_type',
                    "value": ["HTSeq - Counts",
                              "HTSeq - FPKM",
                              "Liftover",
                              "BCGSC miRNA Profiling"]
                    }}]}

    # Main Study Oject
    study = autoGDC(filt = filt,
                    size = 100000)

    df = study.dataframe
    ages = study.metadata["age"]

    differential_methylation = study.ddx()

  """
  # The GDC api url is used to obtain data via the REST API
  #   Useful urls are:
  #     {gdc_api_url}/cases
  #     {gdc_api_url}/data
  #     {gdc_api_url}/files
  api_url = "https://api.gdc.cancer.gov"


  def __init__(self,
               filt: dict = None,
               fields: list = None,
               size: int = 10,
               contrasts: list = None,
               paired_assay: bool = False,
               data_dir: str = None,
               methyl_loci_subset: list = None):

    self.downloader = collate.Downloader()

    self.paired_assay = paired_assay
    self.contrasts = contrasts

    if fields is None:
      fields = ["file_name",
           "md5sum",
           "file_size",
           "archive.state",
           "access",
           "disease_type",
           "submitter_id",
           "primary_site",
           "analysis.workflow_type",
           "cases.case_id",
           "cases.diagnoses.tumor_stage",
           "cases.diagnoses.age_at_diagnosis",
           "cases.samples.tissue_type",
           "cases.samples.sample_type",
           "cases.project.project_id"
           "cases.disease_type",
           "diagnoses.vital_status"]

    if filt is None:
      # Default to Glioma RNA and DNA
      filt = {"op":'and',
              "content":
                [{"op":"IN",
                  "content":
                      {"field": 'cases.disease_type',
                       "value": ['Gliomas',
                                 'Neuroblastoma']}},
                 {"op":"IN",
                  "content":
                      {"field": 'files.analysis.workflow_type',
                       "value": ["HTSeq - FPKM",
                                 "Liftover"]}}]}

    self.params = {"filters" : json.dumps(filt),
             "format" : "tsv",
             "fields" : ",".join(fields),
             "size" : str(size)}

    self.file_types = {"HumanMethylation450":[],
               "HumanMethylation27":[],
               "FPKM.txt":[],
               "FPKM-UQ":[],
               "htseq.counts.gz":[],
               "mirnas.quantification":[],
               "isoforms.quantification":[]}

    # Location of all downloaded data storage
    self.data_dirs = self._create_data_dirs(data_dir)

    self.methyl_loci_subset = methyl_loci_subset

    # Calculates data and metadata
    LOG.debug("Initializing metadata, data, and frame...")
    self._metadata = None
    self._data = None
    self._frame = None
    LOG.debug("Initialization complete.")


  @property
  def metadata(self):
    if self._metadata is None:
      self._metadata = self._get_metadata()
    return self._metadata


  @metadata.setter
  def metadata(self, value):
    self._metadata = value


  @property
  def data(self):
    if self._data is None:
      d = self._get_data()
      self._data = d #self._change_id_columns(data_df = d,
                     #                        meta_df = self.metadata,
                     #                        id_type = "case_id")
    return self._data


  @property
  def frame(self):
    if self._frame is None:
      self._frame = self._get_frame()
    return self._frame


  def _create_data_dirs(self, data_dir):
    # Create / assign download directory
    #   Structured to have a directory for each assay type
    data_dirs = {}

    # Main directory
    if data_dir is None:
      data_dir = path.join("/", "data", "autoGDC")
    data_dir = path.join(path.expanduser(data_dir))
    os.makedirs(data_dir, exist_ok = True)
    data_dirs["main"] = data_dir

    # Raw data / temporary download storage directory
    raw_dir = path.join(data_dir, "raw")
    os.makedirs(raw_dir, exist_ok = True)
    data_dirs["raw"] = raw_dir

    # Directory for each biological assay
    for f in self.file_types.keys():
      f = f.replace(".", "_")
      filetype_dir = path.join(data_dir, f)
      os.makedirs(filetype_dir, exist_ok = True)
      data_dirs[f] = filetype_dir

    return data_dirs


  def _feature_metadata_check(self):

    LOG.debug("Checking feature metadata...")

    # Check if feature metadata exists for each assay
    for assay, assay_dir in self.data_dirs.items():
      if any(s in assay for s in ["Methylation27", "Methylation450"]):
        metadata_path = path.join(assay_dir, "feature_metadata.tsv")

        # If it isn't there, put it there
        #   (using a local git file now - need to use git LFS)
        if not path.exists(metadata_path):
          LOG.info("Moving feature metadata to data directories...")
          std_metadata_fname = f"{assay}_feature_metadata.tsv"
          std_metadata_path = path.join(this_path, "data", std_metadata_fname)
          shutil.copy(std_metadata_path, metadata_path)

    return True


  def _get_metadata(self):
    LOG.info("Searching for data that meets criteria...")

    # Download metadata from GDC
    response = requests.get(f"{self.api_url}/files", params = self.params)
    response_file = io.StringIO(response.content.decode("utf-8"))
    try:
      metadata = pd.read_csv(response_file, sep = "\t")
    except EmptyDataError as e:
      LOG.warn(f"It appears that there were no samples that fit the criteria.")
      LOG.exception(e)


    # Simplify the columns
    metadata.columns = [col.split(".")[-1] for col in metadata.columns]

    # Drop columns containing null string data
    metadata = metadata.applymap(lambda x:
                                 np.nan if any(s == str(x).lower()
                                           for s in ["not reported", "--"])
                                 else x)
    metadata = metadata.dropna(axis = 1, how = "all")

    # Drop columns containing > 90% actual NaN data:
    threshold = 0.9
    num_nan_thresh = threshold * len(metadata)
    metadata = metadata.dropna(axis = 1,
                               thresh = num_nan_thresh)

    # Age is in days - change it to years
    if "age_at_diagnosis" in metadata.columns:
      metadata["age"] = metadata["age_at_diagnosis"]/365

    # The file UUID as index
    #   for comparison when determining download requirements
    metadata = metadata.set_index("id")

    if self.paired_assay:
      metadata = subset_paired_assay(metadata)

    LOG.info(f"Found {metadata.shape[0]} files which match the criteria.")
    return metadata



  def _download(self, file_ids):
    LOG.info("Downloading data... (This may hang and take a while!)")

    # Create manifest file for downloading the files
    manifest = self.metadata.loc[file_ids][["file_name", "md5sum", "file_size"]]
    manifest["state"] = "validated"
    manifest.columns = ["filename", "md5", "size", "state"]

    # Process the manifest in chunks for download
    LOG.info(f"Starting download of {manifest.shape[0]} files")
    chunk_size = 50
    for i in range(0, len(manifest), chunk_size):
      pct_dn = round( i/len(manifest) * 100.0, 2)
      msg = f"Downloading chunk of {chunk_size} files... ({pct_dn}% finished)"
      LOG.debug(msg)

      chunk = manifest.iloc[i: i + chunk_size]

      manifest_file = path.join(self.data_dirs["raw"], "manifest.txt")
      manifest.to_csv(manifest_file, sep = "\t")
      subprocess.call([gdc_client_path,
                       "download",
                       "-m",
                       manifest_file,
                       "--dir",
                       self.data_dirs["raw"]])

    LOG.info("All files have been downloaded.")

    # Now move and rename each file - then extract metadata and gzip
    LOG.info("Reorganizing...")

    # Each directory in the raw directory
    rawpath = self.data_dirs["raw"]
    for d in tqdm(os.listdir(rawpath)):
      dpath = path.join(rawpath, d)

      if not path.isfile(dpath):

        # For each of the files
        #   (The directory should only have 1 or 2 files)
        for f in os.listdir(dpath):
          fpath = path.join(dpath, f)

          # Just the bioassay files (not log directory files)
          if path.isfile(fpath):

            # Find the assay type (in the filename)
            #   make the new filename the file_id
            #   and store the file in the assay directory
            for typ in self.file_types.keys():
              if typ in f:

                # Read files and remove the metadata
                if "Methylation" in typ:
                  # Methylation data is a lot larger due to metadata
                  # Use dtypes and the C engine for the reading of these files
                  dtypes = meth_dtypes
                  header = 0
                  subset_cols = [0, 1]
                  values_name = "beta_value"
                  index_name = "loci"

                elif "quantification" in typ:
                  dtypes = mirna_dtypes
                  header = 0
                  if "isoform" in typ:
                    subset_cols = [0,3]
                  else:
                    subset_cols = [0,2]
                  values_name = "value"
                  index_name = "index"

                else:
                  dtypes = None
                  header = None
                  subset_cols = None
                  values_name = "value"
                  index_name = "index"

                LOG.debug(f"Reading {fpath}")
                try:
                  # Load data (without extra metadata) into memory
                  series = pd.read_csv(fpath,
                                       sep = "\t",
                                       index_col = 0,
                                       header = header,
                                       usecols = subset_cols,
                                       dtype = dtypes,
                                       engine = "c").iloc[:,0]

                  # Change to gene name instead of ENSG
                  if series.index[0].startswith("ENSG"):
                    series.index = [ix.split(".")[0] for ix in series.index]
                    # series.index = series.index.map(gene_id_name)
                    series.index = [gene_id_name.loc[transcript]
                                    if transcript in gene_ensg_set
                                    else transcript
                                    for transcript in series.index]
                    series = series.loc[series.index.dropna()]

                  series.name = values_name
                  series.index.name = index_name

                  # Remove the non compressed file
                  os.remove(fpath)

                  # Save the file as gzipped without extra metadata
                  gz_path = path.basename(path.splitext(
                                            fpath.rstrip(".gz"))[0]) +".tsv.gz"
                  LOG.debug(f"Saving {gz_path}")
                  series.to_csv(gz_path,
                                sep='\t',
                                compression='gzip')

                  # Define new place to move the file
                  #   `d` should be the UUID of the file
                  new_filename = d + ".tsv.gz"
                  new_dir = self.data_dirs[typ.replace(".", "_")]
                  new_path = path.join(new_dir, new_filename)

                  # Move the file
                  LOG.debug(f"Moving {gz_path} to {new_path}")
                  shutil.move(gz_path, new_path)

                  # Remove the whole directory that was downloaded
                  shutil.rmtree(dpath)

                except:
                    LOG.info(f"Error processing file: {fpath}")

    self._feature_metadata_check()

    return True


  def _get_data(self):
    data_dict = {"Transcriptome Profiling": {"miRNA": None,
                                             "isoforms": None,
                                             "counts": None,
                                             "FPKM": None,
                                             "FPKM-UQ": None},
                 "DNA Methylation": {"450": None,
                                     "27": None,
                                     "meta450": None,
                                     "meta27": None}}

    # All file ids in dataframe
    study_ids = set(self.metadata.index.tolist())
    LOG.debug(f"UUIDs for entire study are: {study_ids}")

    # File ids already available locally
    owned_ids = set([path.splitext(file_id.rstrip(".gz"))[0]
                             for assay, assaydir in self.data_dirs.items()
                             if assay not in ["raw", "main"]
                             for file_id in os.listdir(assaydir)
                             if 'metadata' not in file_id])

    LOG.debug(f"Total UUIDs already downloaded are: {owned_ids}")

    # File ids required to download
    download_ids = study_ids - owned_ids
    LOG.info("Downloading these file_ids:" + str(download_ids))
    if download_ids:
      self._download(list(download_ids))

    LOG.info("Creating dataframes...")

    # Get a dataframe for each assay type
    for filetype in self.file_types.keys():
      assaydir_name = filetype.replace(".", "_")
      assaydir = self.data_dirs[assaydir_name]

      # Get all of the files pertaining to this study
      for filename in os.listdir(assaydir):
        filepath = path.join(assaydir, filename)
        uuid = path.basename(path.splitext(filename.rstrip(".gz"))[0])
        if uuid in study_ids:
          self.file_types[filetype] += [filepath]

      # Create dataframe
      LOG.info(f"{assaydir_name} DataFrame...")

      # Caching requires copies of each item
      #   such that inputs are not dependent upon `self`
      fpths = (self.file_types[filetype] + [""])[:-1]
      if "Methylation" in assaydir_name:
        if "450" in assaydir_name:
          assay_subtype = "450"
        else:
          assay_subtype = "27"

        data_dict["DNA Methylation"][assay_subtype] = \
                df_utils.quantile(
                        multifile_df(
                            sorted(fpths),
                            subset_features = self.methyl_loci_subset
                            )
                        )

        # Metadata
        metadata_fpath = path.join(assaydir, "feature_metadata.tsv")
        data_dict["DNA Methylation"]["meta"+assay_subtype] = \
                pd.read_csv(metadata_fpath,
                            sep = "\t",
                            index_col = 0)
      else:
        main_assay_type = "Transcriptome Profiling"
        if "counts" in assaydir_name:
          assay_subtype = "counts"
        elif "FPKM-UQ" == assaydir_name:
          assay_subtype = "FPKM-UQ"
        elif "FPKM_txt" == assaydir_name:
          assay_subtype = "FPKM"
        elif "isoforms" in assaydir_name:
          assay_subtype = "isoforms"
        elif "mirna" in assaydir_name:
          assay_subtype = "miRNA"

        data_dict[main_assay_type][assay_subtype] = \
                                       df_utils.quantile(multifile_df(sorted(fpths)))

    gc.collect()
    LOG.info("Finished constructing data structure...")

    return data_dict


  def _change_id_columns(self, data_df, meta_df, id_type = "case_id"):
    """id_type: the type of id (either case_id or file_id) that you would like your dataframe to use as columns"""
    # Need case IDs to pair data
    # Then we can get the gene expression for each sample as a target
    # label by case id rather than file id
    file_to_case = {file_id:case_id for file_id, case_id in meta_df["case_id"].items()}
    case_to_file = {case_id:file_id for file_id, case_id in meta_df["case_id"].items()}
    for d_type in data_df.keys():
      for platform_type in data_df[d_type].keys():
        if data_df[d_type][platform_type] is not None:
          id_type_now = data_df[d_type][platform_type].columns.name

          if id_type_now == id_type:
            pass
          else:
            if id_type == "case_id":
              data_df[d_type][platform_type] = data_df[d_type][platform_type].rename(columns = file_to_case)
            elif id_type == "file_id":
              data_df[d_type][platform_type] = data_df[d_type][platform_type].rename(columns = case_to_file)
            data_df[d_type][platform_type].columns.name = id_type
    return data_df


  def _filter_loci(self,
                   min_seq = 5,
                   max_seq = 20,
                   collapse_level = "Gene_Symbol", # could be "Transcript_ID"
                   position_col = "Position_to_TSS",
                   pos_bounds = (-1500, 0)):

    # The 450k loci chip should have the coverage we need
    meta = self.data["DNA Methylation"]["meta450"]

    # Assume each loci is associated with the first gene in list
    # Loci with non-null genes
    meta = meta.loc[meta[collapse_level] != ".", :]
    # Loci with only one gene associated
    meta = meta[meta[collapse_level]\
                    .apply(lambda x: len(x.split(";")) == 1)]

    # Change position to integer (not string)
    # If we want all alternative transcripts
    #   with different TSS, we can use the set of genes == 1:
    #     `len(set(x.split(";")))==1`
    meta.loc[:, position_col] = meta[position_col].astype(int)

    # We should also throw out any loci with positions outside of our bounds:
    meta = meta[meta[position_col]\
                     .apply(lambda x:
                             (pos_bounds[0] < x) and (x < pos_bounds[1]))]

    # Use the number of loci within each gene
    #   in order to filter to genes which have a certain coverage of DNAm
    loci_per_gene = meta[collapse_level].value_counts()
    loci_per_gene.name = "number of loci within gene"

    # Filter the genes according to the min/max
    #   desired number of DNAm loci per gene
    filtered_genes = loci_per_gene[(loci_per_gene >= min_seq)
                                   & (loci_per_gene <= max_seq)].index
    filtered_loci = meta[meta[collapse_level].isin(filtered_genes)]

    filt_methdf = self.data["DNA Methylation"]["450"]\
                                .reindex(filtered_loci.index).dropna()

    # change position to be bounded between 0 and 1
#    pos = filtered_loci[position_col]
#    filtered_loci[position_col] = (p - p.min()) / p.max()

    return filt_methdf, filtered_loci


  def _get_frame(self,
                 collapse_level = "Gene_Symbol",
                 agg_func = list,
                 position_col = "Position_to_TSS",
                 min_seq_len = 15,
                 max_seq_len = 25,
                 pos_bounds = (-1500, 1000)):
    """
      This multi-indexed dataframe combines the RNA and DNAm data
    """

    dfs_dict = self.data["Transcriptome Profiling"]
    filt_methdf, filtered_loci = self._filter_loci(min_seq = min_seq_len,
                                               max_seq = max_seq_len,
                                               collapse_level = collapse_level,
                                               position_col = position_col,
                                               pos_bounds = pos_bounds)

    df = df_utils.combined_region_collapsed_frame(dfs_dict = dfs_dict,
                                          main_seq_data = filt_methdf,
                                          seq_feature_metadata = filtered_loci,
                                          region = collapse_level,
                                          agg_func = agg_func,
                                          pos_col = position_col,
                                          max_seq_len = max_seq_len,
                                          value_name = "beta value")
    return df


  def ddx(self,
          contrasts = None,
          formula = None,
          DESeq2 = True):
    if contrasts is None:
      contrasts = self.contrasts
    if formula is None:
      formula = "~"+"+".join(self.contrasts)

    df = self.data["Transcriptome Profiling"]['counts'].astype(int)
    design = self.metadata[contrasts].reindex(df.columns).reset_index()

    if DESeq2:
        DEG = pydeseq(counts = df, design = design, formula = formula)
    else:
      # Characteristic Direction (Multivariate statistical method)
      # 0 excluded, 1 is control, 2 is perturbation
      classes = self.metadata[contrasts]

      # Calculate differential expression / methylation
      DEG = characteristic_direction(data = self.dataframe.values,
                   sampleclass = classes,
                   genes = self.dataframe.index,
                   gamma = 1., # smooths covariance and reduces noise
                   sort = True,
                   calculate_sig = True,
                   nnull = 100,
                   sig_only = True,
                   norm_vector = False)

      DEG = pd.DataFrame(DEG)

    return DEG
