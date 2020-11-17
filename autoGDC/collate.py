# General dependencies
#import os
#import io
#import gc
#import json
#import shutil
#import logging
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
from .config import SETTING
from .autoGDC import LOG
from .download import Downloader


class Collator(Downloader):
  def __init__(self,
               config_key: str = "default",
               params: dict = None):

    super().__init__(config_key = config_key,
                     params = params)
    self.conf = SETTING[config_key]


  def _organize_files(self):
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

    return

