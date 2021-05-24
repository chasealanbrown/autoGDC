import os
#import logging
import pandas as pd
from os import path, listdir

from .config import SETTINGS, LOG
from .collate import Collator

# Logger for collator
#logging.getLogger().setLevel(logging.INFO)
#LOG = logging.getLogger(__name__)


class Archive(Collator):
  """
  Summary:
    Class that keeps track of the inventory of downloaded files

  Arguments:
    data_dir:
      The main directory where all of the downloaded GDC files will be stored
  """

  def __init__(self,
               config_key: str = "default",
               params: dict = None,
               paired_assay: bool = False):

    self.conf = SETTINGS[config_key]

    # Initialize directory structure for data first
    self._init_dirs()

    super().__init__(config_key = config_key,
                     params = params)

#    self._databases = None

  def _read_dataframe(self, assay):
    db_path = self.database_files[assay]
    query = f"file_id in {self.file_ids}"
    data = pd.read_h5(db_path,
                      key = "data",
                      where = query)
    sample_metadata = pd.read_h5(db_path,
                                 key = "sample_metadata",
                                 where = query)
    feature_metadata = pd.read_h5(db_path,
                                  key = "feature_metadata")
    return pd.DataFrame(data.values,
                        index = feature_metadata,
                        columns = sample_metadata)

#  @property
#  def files(self):
#    """
#    Returns:
#      A dictionary containing the lists of files
#        indexed by filetype (directory) keys that have been downloaded
#    """
#    file_map = {}
#
#    # First, check the directories for actual raw files
#    assay_types = self.filetype_regexs.keys()
#    for assay_type in assay_types:
#      assay_dir = self.data_dirs[assaytype]
#      file_map[assay_dir] = [path.join(assay_dir, f)
#                             for f in listdir(assay_dir)]
#
#    # Then add the files that are on the unprocessed dataframes
#    #   for each data type
#    return
#
#
#  @property
#  def uuids(self):
#    """
#    Returns:
#      A dictionary containing the lists of uuids
#        indexed by filetype (directory) keys that have been downloaded
#    """
#    uuid_map = {}
#
#    # This should probably go into Collator, which can then move the files
#    #   and delete them (according to keep_raw parameter)
#    #   Then all that the Archive needs to do is just search the h5 file
#    #   for the associated metadata, once the Downloader and Collator
#    #   are finished
#
##    # First, check the directories for actual raw files
##    assay_types = self.filetype_regexs.keys()
##    for assay_type in assay_types:
##      assay_dir = self.data_dirs[assaytype]
##      uuid_map[assay_type] = [basename(splitext(]
##
##    for assaydir in file_map
##    for filename in listdir(
##    # Then add the files that are on the unprocessed dataframes
##    #   for each data type
#    return

  def _init_dirs(self):
    LOG.info("This appears to be the first run.\nInitializing directories...")
    # Initialize main directories
    dirs = ["data_dir", "archive_dir", "newdata_dir",
            "rawdata_dir", "mygene_dir", "metadata_dir"]
    for d in dirs:
      os.makedirs(self.conf[d], exist_ok = True)

    # Directory for each biological assay
    for assay, fp in self.conf["assay_dir"].items():
      os.makedirs(fp, exist_ok = True)

    return
