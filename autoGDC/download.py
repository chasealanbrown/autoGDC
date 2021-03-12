"""
Module for downloading GDC data
"""
import re
from os import path
from io import StringIO
import json
import requests
import requests_cache
import pandas as pd
import numpy as np
from tqdm import tqdm
from urllib.parse import urlparse

# Gene information
import mygene

# Local modules
from .config import SETTINGS, LOG
from .df_utils import column_expand, subset_paired_assay, metadata_json_to_df

import logging
LOG.setLevel(logging.DEBUG)

# Feature metadata urls for methylation
drpbx = "https://www.dropbox.com/s/"
_end_url = "_feature_metadata.tsv?dl=1"
DNAm_27k_url = f"{drpbx}o36utk95fl1euy1/HumanMethylation27{_end_url}"
DNAm_450k_url = f"{drpbx}wyck0lsa6941utw/HumanMethylation450{_end_url}"

def download_url(url, filepath, params = None):
  """
  Summary:
    Download a file from a specific url with a progress bar.
  Source:
    https://stackoverflow.com/questions/56795227/how-do-i-make-progress-bar-while-downloading-file-in-python#56796119
  """
  parsed = urlparse(url)
  filename = path.basename(parsed.path)
  LOG.debug(f"Downloading {filepath} from {url} ...")
  # TODO:
  #   If the url points to a directory, this may need to be handled
  #     by a different function to get all of the files.
  #     However - here - we *assume* that if `filepath` is a directory,
  #     then this is simply the location to store the data,
  #     and we use the filename from the url to save it.
  if path.isdir(filepath):
    file_dir = filepath
    filepath = path.join(file_dir, filename)

  with requests.get(url, stream=True, params = params) as r:
    r.raise_for_status()
    with open(filepath, 'wb') as f:
      datasize = r.headers.get("Content-Length")
      datasize = r.headers.get("content-length") if datasize is None else datasize
      pbar = tqdm(total = None if datasize is None else int(datasize),
                  desc = f"Downloading {filename}...")
      for chunk in r.iter_content(chunk_size = 8192):
        if chunk:  # filter out keep-alive new chunks
          f.write(chunk)
          pbar.update(len(chunk))
  return filepath


def download_post(url: str, data: dict, headers: dict, download_dir: str = "data/"):
  """
  Summary:
    Download data using POST on a REST API with a progress bar.
  """
  LOG.info(f"Downloading {data} from {url} ...")
  try:
    with requests.post(url, data = json.dumps(data), headers = headers) as r:
      r_head_cd = r.headers.get("Content-Disposition")
      r_head_cd = r.headers.get("content-disposition") if r_head_cd is None else r_head_cd
      file_name = re.findall(r"filename=(.+)", r_head_cd)[0]

      fp = path.join(download_dir, file_name)
      with open(fp, "wb") as output_file:
        output_file.write(r.content)
  except Exception as e:
    LOG.error(f"Critical Error when downloading {data} from {url}")
    LOG.error(e)
  return True



class Downloader(object):
  """
  Summary:
    Class to download data from GDC

  Arguments:
    Many examples of how the data is accessed can be seen here:
      https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/

    config_key:
      A string which denotes the type of local system configuration desired.
      Typical values are 'default' or 'wrenlab'.
      Other configurations can be created in the file `config.py`

    params:
      Parameters to query GDC via REST API

  Attributes:
    file_ids:
      file_ids to be downloaded from GDC

    metadata:
      The table of GDC metadata for this specific download

    manifest:
      The manifest table for this download.  This includes all file_ids that
        are desired for this specific download object.
        (This table is split up by the _download() function into chunks)

    metadatabase:
      The full table of all GDC metadata.
  """


  # The GDC api url is used to obtain data via the REST API
  #   Useful urls are:
  #     {api_url}/cases
  #     {api_url}/data
  #     {api_url}/files
  #     {api_url}/status
  #     {api_url}/gql/_mapping
  api_url = "https://api.gdc.cancer.gov"

  def __init__(self,
               config_key: str = "default",
               params: dict = None,
               paired_assay: bool = False):

    self.conf = SETTINGS[config_key]
    self.status_filepath = path.join(self.conf["data_dir"], "gdc_status")
    requests_cache.install_cache(path.join(self.conf["newdata_dir"], "autoGDC_dl"),
                             backend='sqlite',
                             expire_after=60*60*24*7) # A week in seconds

    self.params = params
    self.paired_assay = paired_assay
    self.updated_index = False
    self._gdc_filelist_url = None
    self._gdc_status = None
    self._gdc_fields = None
    self._file_ids = None
    self._owned_file_ids = None
    self._new_file_ids = None
    self._metadata = None
    self._manifest = None
    self._metadatabase = None
    self._metadatabase_index = None

    self.mg = mygene.MyGeneInfo()
    self.mg.set_caching(cache_db = path.join(self.conf["mygene_dir"], "mygene"),
                        verbose = False)
    self.mg.delay = 0

  @property
  def gdc_status(self):
    """
    Summary:
      The status of the remote GDC data.
        Keeping track of this will allow us to download less data from GDC,
        as we can check the local metadatabase if our version is up to date
    """
    if self._gdc_status is None:
      # TODO: Add an 'expiration date' to only check the status every so often,
      #       rather than every single time a Download object is instantiated.
      with requests.get(f"{self.api_url}/status") as r:
        status = r.json()
      self._gdc_status = status
    return self._gdc_status


  @property
  def local_status(self):
    """The status of the local GDC data."""
    try:
      with open(self.status_filepath) as f:
        status = json.loads(f.read())
    except FileNotFoundError:
      # The reason that this file does not exist, is likely because
      #   the project is just being initiated - therefore
      #   the solution is likely to download the master index table from GDC
      # However, we should not begin that process here, as it would cause an
      #   infinite loop - when retrieving the master index file,
      #   this property (`local_status`) is used
      # So instead we just give None
      return None
    return status


  @property
  def gdc_filelist_url(self):
    """
    Summary:
      The url for GDC's master list of all available files.
        The table contains mostly just index info, such as:
        Filenames, UUIDs, md5sum, etc.
    """
    if self._gdc_filelist_url is None:
      _base_url = "https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/"
      # The first table (index 0 in this list) from the url above
      #   gives the gdc version releases info
      release_tbl = pd.read_html(_base_url)[0]
      # The lastest info is on top
      newest = release_tbl.iloc[0]
      # remove 'v' (given as e.g. 'v27.0')
      _version = newest.Version.strip("v")
      _date = pd.to_datetime(newest_info.Date).strftime("%Y%m%d")
      _file_str = f"gdc_manifest_{_date}_data_release_{_version}_active.tsv.gz"
      self._gdc_filelist_url = f"{_base_url}{_file_str}"
    return self._gdc_filelist_url


  # This should be a part of config? or GDC_config?
  # Perhaps a GDCArchive class or something?
  @property
  def gdc_fields(self):
    if self._gdc_fields is None:
      url = f"{self.api_url}/gql/_mapping"
      with requests.get(url) as r:
        df = pd.DataFrame(json.loads(r.content.decode("utf-8")))
      self._gdc_fields = df
    return self._gdc_fields


  @property
  def file_ids(self):
    """
    Summary:
      The file_ids associated with the dataset.
        This includes both the file_ids that will be downloaded (new_file_ids),
        as well as file_ids that are already locally stored (owned_file_ids)
    """
    if self._file_ids is None:
      self._file_ids = self._get_file_ids()
    return self._file_ids


  @property
  def owned_file_ids(self):
    """
    Summary:
      The file_ids already downloaded from GDC and available locally
    """
    if self._owned_file_ids is None:
      self._owned_file_ids = self.metadatabase[self.metadatabase.downloaded].index.tolist()
    return self._owned_file_ids


  @property
  def new_file_ids(self):
    """
    Summary:
      The file_ids to be downloaded from GDC
    """
    if self._new_file_ids is None:
      # TODO: Add ability to use GDC_API.key

      # Remove file_ids that are not open access
      #   We may also use this if `acl` column: == "[b'open']"]
      open_id_bool = self.metadatabase.acl.apply(lambda x: "open" in x)
      locked_ids = self.metadatabase[~open_id_bool].index.tolist()

      self._new_file_ids = list(set(self.file_ids) # File_ids in this study
                                - set(self.owned_file_ids) # Already downloaded
                                - set(locked_ids)) # Controlled access files

      # TODO: Not sure what is happening here, but it seems that somtimes
      #    ther is a NaN or None in this list.  So this is another band-aid
      self._new_file_ids = [fid for fid in self._new_file_ids
                            if type(fid) is str]
    return self._new_file_ids


  @property
  def metadata(self):
    """
    Summary:
      The table of GDC metadata for this specific download
    """
    if self._metadata is None:
      self._metadata = self.metadatabase.query("id == @self.file_ids")#self.metadatabase.loc[self.file_ids]
      if self.paired_assay:
        # TODO: This needs to be fixed by placing it in params
        #       - if size is ~50, you won't get anything - becuase there aren't
        #         enough file to get a match.  This is wrong.
        #         We have to increase the size before getting to this point
        self._metadata = subset_paired_assay(self._metadata)

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
    LOG.info("""Be aware - setting the metadata prior to calling the
    `study.dataframes` property will alter the filtering and
    subsetting of `data.dataframes`.  You can use this functionality to reduce
    the data downloaded/subset further than the original filter.""")
    self._metadata = value


  @property
  def manifest(self):
    """
    Summary:
      The manifest table for this download.  This includes all file_ids that
        are desired for this specific download object.
        (This table is split up by the _download() function into chunks)
    """
    if self._manifest is None:
      self._manifest = self._get_manifest()
    return self._manifest


  @property
  def metadatabase(self):
    """
    Summary:
      The full table of all GDC metadata.
    """
    if self._metadatabase is None:
      self._metadatabase = self._get_full_metadatabase()
    return self._metadatabase


  @property
  def metadatabase_index(self):
    """
    Summary:
      The index table for the GDC metadata.
      New versions are released every so often.
    """
    if self._metadatabase_index is None:
      self._metadatabase_index = self._get_metadatabase_index()
    return self._metadatabase_index


  def _get_metadatabase_index(self) -> pd.DataFrame:
    """
    Summary:
      A helper function to download the main index of files in GDC as a table.
    """

    # TODO:
    #   Check version of GDC database and compare it ours
    outdated = self.local_status != self.gdc_status
    if outdated:
      LOG.info("Updating the index for all GDC files...")

      # Table provided by GDC of filenames, UUIDs, md5sum, etc.
      download_url(url = self.gdc_filelist_url,
                   filepath = self.conf["metadatabase_index_path"])

      # Update local status to the value of GDC status
      # TODO:
      #   This should only be run if we get confirmation of proper download
      #   e.g. md5sum check?
      with open(self.status_filepath, 'w') as f:
        json.dump(self.gdc_status, f)

      self.updated_index = True

    metadatabase_ix = pd.read_csv(self.conf["metadatabase_index_path"],
                                  sep = '\t',
                                  compression = "gzip").set_index("id")
    return metadatabase_ix


  def metadata_and_index_synced(self):
    """Check to see if the indices are synced for the full metadatabase
    and the metadata master index table"""
    return self.metadatabase_index.index == self.metadtabase.index


  def _get_full_metadatabase(self) -> pd.DataFrame:
    """
    Summary:
      The full table of all GDC metadata.
    """
    # TODO: Figure out how to use the json filtration parameters to query
    #       locally stored metadatabase, rather than querying gdc every time.


    # If we just updated the master index file, or if the file doesn't exist,
    #   then we need to create this full metadata database
    index_missing = not path.exists(self.conf["metadatabase_index_path"])
    metadb_missing = not path.exists(self.conf["metadatabase_path"])
    if index_missing | metadb_missing | self.updated_index:
      # Now download clinical and demographic metadata for each file
      mdbix = self.metadatabase_index

      sample_metadata_frames = []
      # This must be done in chunks (server errors arise if larger sizes used)
      chunk_size = 10000
      total = mdbix[mdbix.category == "data_file"].shape[0]
      for pos in tqdm(range(0, total, chunk_size), desc = "Updating metadata"):

        # Payload for query
        #   json is used instead of tsv, because there are some files that
        #   end up with lists of values that cause huge sparse tables
        #   (e.g. cases.249,submitter_id, cases.250.submitter_id, etc.)
        #   This is fixed after getting the json to put into a table.
        full_metadata_params = {"format": "json",
                                "fields": ",".join(self.conf["fields"]),
                                "size": str(chunk_size),
                                "from": str(pos)}
        to = pos+chunk_size
        fp = path.join(self.conf["data_dir"],
                       "metadata",
                       f"metadata{pos}-{to}.json")
        download_url(f"{self.api_url}/files", fp, full_metadata_params)

        metadatachunk = metadata_json_to_df(fp,
                           null_type_strings = self.conf["null_type_strings"])
        sample_metadata_frames += [metadatachunk]

      LOG.info("Joining the metadata chunks together...")
      sample_metadata_frames = pd.concat(sample_metadata_frames)
      # This needs to be merge - it allows overlapping columns to 'merge' into
      #   the same column without the need for suffixes and prefixes for
      #   'right' or 'left' columns
      #   The use of outer just allows the places where a perfect match in
      #   mostly matching columns will just add NaNs to the additional columns
      # TODO:
      #   There is a problem even after using `merge` with having columns
      #     such as `acl` (which should be exactly the same between the
      #     metadata master index table and the metadata json information
      #     retruned by GDC's REST API) having NaNs in (presumably) the json
      #     information.
      #   The quick fix for this is just to remove the overlapping columns from
      #   the json information dataframes and assume the metadata master index
      #   table is correct and not the json information from the REST API, but
      #   this certainly needs a closer look.
      ix_cols = set(mdbix.columns.tolist())
      mdf_cols = set(sample_metadata_frames.columns.tolist())
      overlap_cols = ix_cols.intersection(mdf_cols) - {"id"}
      indices = ["id"] + [col for col in sample_metadata_frames
                          if "nested" in col]
      # The index for mdbix is already set to `id`
      #    The `.reset_index()` is needed to use `.set_index()` later
      #    and include `id`
      metadatabase = mdbix.join(sample_metadata_frames.reset_index()\
                                    .drop(overlap_cols, axis=1)\
                                    .set_index(indices),
                                 how="outer")

      # Include tracker of what is downloaded and not
      # TODO:
      #   This is clearly wrong - if we just update,, the data will still be
      #   there ... but it could be re-indexed in a completely different way?
      #   Probably should be handled by store.py somehow?
      #   Redownloading data every time an update to status occurs is not ideal
      metadatabase["downloaded"] = False

      # Correct `acl` Access Control List
      metadatabase["acl"] = metadatabase.acl.apply(lambda x:
                                                   set([i.decode("utf-8")
                                                        for i in eval(x)]))

      LOG.info("Saving the new metadata database...")
      metadatabase.to_csv(self.conf["metadatabase_path"],
                          sep = '\t',
                          compression = "gzip")

    # Otherwise, the file does exist, so we can just load it
    try:
      # TODO:
      #   The real solution to this is to have a whole module that can track
      #     the versioning of GDC and match this important index file with
      #     the main GDC database.
      LOG.info("Reading the metadata database...")
      metadatabase = pd.read_csv(self.conf["metadatabase_path"],
                                          sep = '\t',
                                          compression = "gzip")

      # Set indices to file UUID and enumeration of nested column vals
      indices = ["id"] + [col for col in metadatabase.columns
                          if "nested" in col]
      metadatabase = metadatabase.set_index(indices)

      LOG.debug("The metadata database is loaded:\n {metadatabase}")

      # If the database for metadata has been properly constructed, it will
      #   have the proper columns
      assert "case_id" in metadatabase.columns

    except:
      LOG.warn("Major Error: Metadatabase is not found or is in wrong format.")
      raise

    return metadatabase


  def _get_file_ids(self):
    """
    Summary:
      Get file_ids (from associated study) from GDC
    """
    LOG.info("Searching GDC for file_ids that meet criteria...")

    # Change the fields to only include the file_id
    #   - We will get the rest of the metadata from the full metadata database
    #     (metadatabase property)
    params = self.params
    params["fields"] = "file_id"

    # TODO:
    #     This is a band-aid fix to problem small sizes on paired data
    if self.paired_assay:
      params["size"] = f"{10**9}"
    response = requests.get(f"{self.api_url}/files", params = params)
    response_file = StringIO(response.content.decode("utf-8"))
    try:
      response_df = pd.read_csv(response_file, sep = "\t")
      LOG.debug(f"""Query for GDC file_ids returned the following dataframe:
                \n{response_df}""")
      file_ids = response_df["id"].tolist()
    except ValueError as e:
      LOG.warn("It appears that there were no samples that fit the criteria.")
      LOG.exception(e)
      file_ids = []
    return file_ids


  def _get_manifest(self):
    """
    Summary:
      Retrieves the manifest table to store as a property for this download.
    """
    # Create *complete* manifest file for downloading the files
#    boolmask = self.metadatabase.ids.isin(self.new_file_ids)
    ix = self.metadatabase.query("id == @self.new_file_ids").index
    manifest = self.metadatabase.loc[ix][["file_name",
                                          "md5sum",
                                          "file_size"]]
    #manifest = self.metadatabase.loc[self.new_file_ids][["file_name",
    #                                                     "md5sum",
    #                                                     "file_size"]]
    manifest["state"] = "validated"
    manifest.columns = ["filename", "md5", "size", "state"]
    return manifest


  def _download_feature_metadata(self):
    """
    Summary:
      Downloads feature metadata (such as gene features for CpG loci)
      if it is not already present
    """
    LOG.debug("Check for DNA methylation feature metadata...")

    # Check if feature metadata exists for each assay
    for assay, assay_dir in self.conf["assay_dirs"].items():
      if "DNAm" in assay:
        if assay == "DNAm_450":
          url = DNAm_450k_url
          metadata_path = self.conf["DNAm_450k_metadata_filepath"]
        elif assay == "DNAm_27":
          url = DNAm_27k_url
          metadata_path = self.conf["DNAm_27k_metadata_filepath"]
        else:
          LOG.error(f"Warning: Unforseen error - assay of {assay}")

        # If it isn't there, put it there
        #   (using a local git file now - need to use git LFS or ipfs)
        if not path.exists(metadata_path):
          LOG.info("Downloading feature metadata to data directories...")
          download_url(url, metadata_path)
      elif "RNA" in assay:
        # TODO:
        #   Download information from mygene and store it
        #   Use MyGene to gather data and store it in a pd.DataFrame
        self.mg
    return True


  def _download_data(self, chunk_size: int = 50):
    """
    Summary:
      Download list of files using the gdc-client using manifest files.
        This is done by iteratively making a series of manifest files
        as chunks of a larger manifest.

    Arguments:
        chunk_size:
          An integer to denote how many files to include in each iteration of
            the download process.
    """
    LOG.info("Downloading data... (This may hang and take a while!)")

    # Process the manifest in chunks for download
    LOG.info(f"Starting download of {self.manifest.shape[0]} files")

    for position in tqdm(range(0, len(self.new_file_ids), chunk_size)):
      chunk = self.new_file_ids[position: position + chunk_size]
      file_path = download_post(url = f"{self.api_url}/data",
                    data = {"ids": chunk},
                    headers = {"Content-Type": "application/json"},
                                download_dir=self.conf["newdata_dir"])

#    for position in tqdm(range(0, len(self.manifest), chunk_size)):
#      # Debugging
#      pct_dn = round( position/len(self.manifest) * 100.0, 2)
#      msg = f"Downloading chunk of {chunk_size} files... ({pct_dn}% finished)"
#      LOG.debug(msg)
#
#      chunk = self.manifest.iloc[position: position + chunk_size]
#
#      temp_manifest_file = path.join(self.config["data_dir"],
#                                     "temporary_manifest.txt")
#      chunk.to_csv(temp_manifest_file, sep = "\t")
#      subprocess.call([gdc_client_path,
#                       "download",
#                       "-m",
#                       temp_manifest_file,
#                       "--dir",
#                       self.config["newdata_dir"]])
#
    LOG.info("All files have been downloaded.")
    self._download_feature_metadata()

    _ix = self.metadatabase.query("id == @self.new_file_ids").index
    self.metadatabase["downloaded"].loc[_ix] = True
    self.metadatabase.to_csv(self.conf["metadatabase_path"],
                             sep = '\t',
                             compression = "gzip")
    return True

