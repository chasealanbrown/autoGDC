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
from store import Archive
from config import SETTING

# Logger for autoGDC
logging.getLogger().setLevel(logging.INFO)
LOG = logging.getLogger(__name__)


class Dataset:
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

  def __init__(self,
               config_key: str = "default",
               filt: dict = None,
               fields: list = None,
               size: int = 10,
               contrasts: list = None,
               paired_assay: bool = False):

    self.conf = SETTING[config_key]
    if fields is None:
      fields = self.conf["fields"]
    else:
      fields = fields
    self.params = {"filters" : json.dumps(filt),
                   "format" : "tsv",
                   "fields" : ",".join(fields),
                   "size" : str(size)}

    self.paired_assay = paired_assay
    self.contrasts = contrasts

    self.archive = Archive(config_key = config_key,
                           params = self.params)

    self._metadata = None
    self._data = None
    self._frame = None


  @property
  def metadata(self):
    if self._metadata is None:
      self._metadata = self._get_metadata()
    return self._metadata


  @metadata.setter
  def metadata(self, value):
    """
    Summary:
      Sample metadata setter is provided in order to allow additional filtering
      steps to occur outside of this class and set here.

      This can be done prior to gathering the data into the main dataframe,
      as the first step is to define the sample metadata, then to use
      this sample metadata to construct the dataframes from the archive
    """
    self._metadata = value


  @property
  def data(self):
    if self._data is None:
      self._data = self._get_data()
    return self._data


  @property
  def frame(self):
    if self._frame is None:
      self._frame = self._get_frame()
    return self._frame


  def _get_data(self):
    data_dict = {"RNA": {"miRNA": None,
                         "isoforms": None,
                         "counts": None,
                         "FPKM": None,
                         "FPKM-UQ": None},
                 "DNAm": {"450": None,
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

