# General dependencies
import os
#import io
#import gc
#import json
import shutil
import tarfile
#import logging
#import hashlib
#import requests
#import subprocess
#import numpy as np
import pandas as pd

## Specific function imports
from os import path
from tqdm import tqdm
from io import StringIO
#from joblib import Memory
#from geode import chdir as characteristic_direction

## Local module imports
from .config import SETTINGS, LOG
from .download import Downloader

# Logger for collator
#logging.getLogger().setLevel(logging.INFO)
#LOG = logging.getLogger(__name__)


class Collator(Downloader):
  def __init__(self,
               config_key: str = "default",
               params: dict = None,
               paired_assay: bool = False):

    super().__init__(config_key = config_key,
                     params = params,
                     paired_assay = paired_assay)
    self.conf = SETTINGS[config_key]


  @property
  def database_files(self):
    """
    Returns:
      A dictionary containing the file paths of the databases for each assay.
        (combined gdc files with the same assay type)
        (indexed by assay type keys)
    """
    return self.conf["db_files"]


  def _update(self):
    self._download_data()
    newdata = self._organize()
    self._update_databases(newdata)
    return True


  def _update_databases(self, newdata: dict):

    for assay_type, df in newdata.items():
      db_path = self.database_files[assay_type]

      # First initialize the databases (if needed), storing the
      #   associated metadata in the database first
      metadata_fp = self.conf[f"{assay}_metadata_filepath"]
      if not path.isfile(db_path):
        feature_metadata = pd.read_csv(metadata_fp, sep='\t')
        # Metadata for features should be static
        #   so this will only be stored once
        # The mode is "w" for 'write' here, opposed to "a" (append)
        #   becuase the "write" is applied file-wide. (i.e. rewrites the file)
        feature_metadata.to_hdf(db_path,
                                key = "feature_metadata",
                                mode = "w",
                                data_columns = True,
                                format = "table")

        # This mode *can* be dangerous to use elsewhere
        #   because "a" will rewrite certain dataframes (overwrites keys)
        append = False
        mode = "a"

      else:
        # The "r+" mode lets rows to be appended to keyed dataframes in the H5
        append = True
        mode = "r+"

      # Samples are read one at a time and added as columns
      # Here, we need to have samples as the rows for easier indexing in HDF5
      # TODO: Make sure this is appropriate, or use better logic
      df = df.T

      # The data is indexed by file ids, similar to metadb
      sample_metadata = self.metadb.loc[df.index]

      df.to_hdf(db_path,
                key = "data",
                append = append,
                mode = mode,
                data_columns = True,
                format = "table")
      sample_metadata.to_hdf(db_path,
                             key = "sample_metadata",
                             append = append,
                             mode = mode,
                             data_columns = True,
                             format = "table")
    return True


  def _organize(self) -> dict:
    """
    Summary:
      This method iterates through the newly downloaded files,
        which are provided as a convoluted directory structure from GDC.
        The method does two steps simulateously (to save disk IO time):
          1) Reads the data files, organizing them in memory into a
             dictionary of dataframes (by assay type)
          2) Moves (or deletes*) the read files from the `newdata` location
             into the `rawdata` location.

      * If the configuration file has `keep_raw` as False, the raw files
        wil be deleted, and only the HDF5 files saved from the resulting
        resulting dataframes will be kept

    Returns:
      Dictionary of newly downloaded data.
      Provided as a dictionary of dataframes for each assay type.
    """

    newdata = {assaytype: None for assaytype in self.conf["filetype_regex"]}

    # Now move/delete each file - and get the
    LOG.info("Organizing downloaded files...")

    # Each directory in the raw directory
    newdata_path = self.conf["newdata_dir"]
    for d in tqdm(os.listdir(newdata_path)):
      dpath = path.join(newdata_path, d)

      # Tarred and g-zipped files hold several series to be extracted
      if path.isfile(dpath) and dpath.endswith("tar.gz"):
        with tarfile.open(dpath, "r:gz") as tar:
          for member in tar.getmembers():
            f = tar.extractfile(member)
            if f is not None:
              series = read_series(fpath = member.name,
                                   file_content = f.read())
              # Add series to dataframe dictionary of new data
              if newdata[assay_type] is None:
                newdata[assay_type] = series.to_frame()
              else:
                newdata[assay_type].join(series.to_frame())

      # Directories (non-compressed) are treated differently
      if not path.isfile(dpath):
        # For each of the files
        #   (The directory should only have 1 or 2 files)
        for fname in os.listdir(dpath):
          fpath = path.join(dpath, fname)
          # Just the bioassay files (not log directory files)
          if path.isfile(fpath):
            series = read_series(fpath)
            # Add series to dataframe dictionary of new data
            if newdata[assay_type] is None:
              newdata[assay_type] = series.to_frame()
            else:
              newdata[assay_type].join(series.to_frame())

      # Remove or move the whole directory that was downloaded
      if not self.config["keep_raw"]:
        LOG.debug(f"Removing {dpath}")
        shutil.rmtree(dpath)
      else:
#         assay_dir = self.conf["assay_dirs"][assay_type]
#         new_path = path.join(assay_dir, dpath)
        new_path = path.join(self.conf["rawdata_dir"], dpath)
        LOG.debug(f"Moving {dpath} to {new_path}")
        shutil.move(dpath, new_path)

    return newdata


def read_series(fpath: str, file_content: str = None) -> pd.Series:
  """
  Read a series from a file or a tarfile compressed file (member)

  Arguments:
    fpath:
      This filepath argument is actually a bit nefariously named
      - it can either be a file path directly, or a path *within a tarfile.

    file_content:
      If the fpath is derived from a tar, then this will hold the content of
      the file
  """
  if file_content is not None:
    fpath = StringIO(file_content)
    fname = fpath
  else:
    fname = path.basename(fpath)

  # Find the assay type (in the filename)
  #   make the new filename the file_id
  #   and store the file in the assay directory
  try:
    assay_type = [assay_type for assay_type, regex
                             in self.conf["filetype_regexs"].items()
                             if regex.match(fname)][0]
  except ValueError as e:
    LOG.warn("""The regular expression for identifiying assay type
             from file names appears to have failed""")
    LOG.error(f"Error reading file: {fpath}\n{e}", exc_info=True)

  # Get parameters used for reading data from config
  readparams = self.conf["assaytype_readparams"][assay_type]
  dtypes = readparams["dtypes"]
  header = readparams["header"]
  subset_cols = readparams["subset_cols"]
  values_name = readparams["values_name"]
  index_name = readparams["index_name"]

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

    file_id = self.metadb[self.metadb.file_name == fname].index[0]
    series.name = file_id
    series.index.name = index_name

  except Exception as e:
    LOG.error(f"Error reading file: {fpath}\n{e}", exc_info=True)
  return series

