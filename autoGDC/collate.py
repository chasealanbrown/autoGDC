# General dependencies
#import os
#import io
#import gc
#import json
#import shutil
import logging
#import hashlib
#import requests
#import subprocess
#import numpy as np
#import pandas as pd
#
## Specific function imports
#from os import path
#from tqdm import tqdm
#from joblib import Memory
#from geode import chdir as characteristic_direction
#
## Local module imports
from config import SETTING
#from . import df_utils
from .store import Archive

# Gene information
import mygene

# Logger for collator
logging.getLogger().setLevel(logging.INFO)
LOG = logging.getLogger(__name__)


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
    series = series.reindex(subset_featuees)
  return series


#@memoize
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


class Downloader:
  # The GDC api url is used to obtain data via the REST API
  #   Useful urls are:
  #     {gdc_api_url}/cases
  #     {gdc_api_url}/data
  #     {gdc_api_url}/files
  api_url = "https://api.gdc.cancer.gov"

  def __init__(self,
               data_dir: str = None,
               keep_raw: bool = False):
    self.archive = Archive(self.data_dir)

class Collator:
    def __init__(self, config_key: str = "default"):
      self.config = CONFIG[config_key]
      self.mg = mygene.MyGeneInfo()
      self.mg.set_caching(cache_db = self.config["mygene_dir"])
