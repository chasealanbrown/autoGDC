#import re
#import json
#import requests
#import pandas as pd
#import numpy as np
#from tqdm import tqdm
#from os import path
#from io import StringIO
#
## Gene information
#import mygene
#
## Local modules
#from .config import SETTINGS, LOG
#from .df_utils import column_expand, subset_paired_assay
#
#
#def download_url(url, filepath, params = None):
#  """
#  Summary:
#    Download a file from a specific url with a progress bar.
#  Source:
#    https://stackoverflow.com/questions/56795227/how-do-i-make-progress-bar-while-downloading-file-in-python#56796119
#  """
#  LOG.info(f"Downloading {filepath} from {url} ...")
#  with requests.get(url, stream=True, params = params) as r:
#    r.raise_for_status()
#    with open(filepath, 'wb') as f:
#      datasize = r.headers.get("Content-Length")
#      datasize = r.headers.get("content-length") if datasize is None else datasize
#      pbar = tqdm(total = None if datasize is None else int(datasize))
#      for chunk in r.iter_content(chunk_size = 8192):
#        if chunk:  # filter out keep-alive new chunks
#          f.write(chunk)
#          pbar.update(len(chunk))
#  return
#
#
#def download_post(url: str, data: dict, headers: dict):
#  """
#  Summary:
#    Download data using POST on a REST API with a progress bar.
#  """
#  LOG.info(f"Downloading {data} from {url} ...")
#  with requests.post(url, data = json.dumps(data), headers = headers) as r:
#    r_head_cd = r.headers.get("Content-Disposition")
#    file_name = re.findall("filename=(.+)", r_head_cd)[0]
#
#    with open(file_name, "wb") as output_file:
#      output_file.write(r.content)
#  return True



class DownloadMixin(object):
  """
  Summary:
    Class to add a trait that will allow other objects to download GDC data

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
  #     {gdc_api_url}/cases
  #     {gdc_api_url}/data
  #     {gdc_api_url}/files
  api_url = "https://api.gdc.cancer.gov"

  def __init__(self,
               config_key: str = "default",
               params: dict = None,
               paired_assay: bool = False):

    self.conf = SETTINGS[config_key]
    self.params = params
    self.paired_assay = paired_assay
#    self._version = None
    self._file_ids = None
    self._owned_file_ids = None
    self._new_file_ids = None
    self._metadata = None
    self._manifest = None
    self._metadatabase = None

    self.mg = mygene.MyGeneInfo()
    self.mg.set_caching(cache_db = path.join(self.conf["mygene_dir"], "mygene"))


# Need to have gdc_version and archive_version (local) to compare
#  @property
#  def version(self):
#    """
#    Summary:
#      The version of the GDC data.  Keep track of this will allow us to
#        download less data from GDC, as we can use the locally downloaded
#        metadatabase if our version is up to date
#    """
#    status_filepath = path.join(self.conf["data_dir"], "gdc_version")
#    with open(status_filepath) as f:
#      statusnow = json.loads(f.read())
#
#    if self._version is None:
#      # TODO: Add an 'expiration date' to only check the status every so often,
#      #       rather than every single time a Download object is instantiated.
#      with requests.get(f"{self.gdc_api_url}/status") as r:
#        statustxt = r.content.decode("utf-8")
#        status = json.loads(statustxt)
#        with open(status_filepath) as f:
#          f.write(statustxt)
#      self._version = version
#
#      # TODO: Figure out how to use the json filtration parameters to query
#      #       locally stored metadatabase, rather than querying gdc every time.
#    return self._version


  # This should be a part of config? or GDC_config?
  # Perhaps a GDCArchive class or something?
  @property
  def gdc_fields(self):
    import requests

    url = f"{self.api_url}/gql/_mapping"

    with requests.get(url) as r:
      df = pd.DataFrame(json.loads(r.content.decode("utf-8")))
    return df

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
      locked_ids = self.metadatabase[self.metadatabase.acl == "[b'open']"].index.tolist()
      self._new_file_ids = list(set(self.file_ids) - set(self.owned_file_ids)
                                - set(locked_ids))
    return self._new_file_ids


  @property
  def metadata(self):
    """
    Summary:
      The table of GDC metadata for this specific download
    """
    if self._metadata is None:
      self._metadata = self.metadatabase.loc[self.file_ids]
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


  def _get_full_metadatabase(self,
                             nan_threshold: float = 1.0) -> pd.DataFrame:
    """
    Summary:
      The full table of all GDC metadata.
    """
    # If we have already pulled the database, just open it
    #   otherwise download and store it
    if path.exists(self.conf["metadatabase_path"]):
      possible_metadatabase = pd.read_csv(self.conf["metadatabase_path"],
                         sep = '\t',
                         compression = "gzip").set_index("id")

      # If the database for metadata has been properly constructed, it will
      #   have the proper columns
      if "cases" in possible_metadatabase.columns:
        return possible_metadatabase


    # Otherwise, create the missing database

    # First we get the table provided by GDC
    #   of filenames, UUIDs, md5sum, etc.
    LOG.info("Downloading list of all files in the GDC...")
    download_url(url = self.conf["gdc_filelist_url"],
                 filepath = self.conf["metadatabase_path"])
    metadatabase = pd.read_csv(self.conf["metadatabase_path"],
                               sep = '\t',
                               compression = "gzip").set_index("id")

    # Now download clinical and demographic metadata for each file
    LOG.info("Downloading metadata for all files in the GDC...")

    sample_metadata_frames = []
    # This must be done in chunks (server errors arise if larger sizes used)
    chunk_size = 10000
    total = metadatabase[metadatabase.category == "data_file"].shape[0]
    for pos in tqdm(range(0, total, chunk_size)):

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
      fp = path.join(self.conf["data_dir"], "metadata", f"metadata{pos}-{to}.json")
#        download_url(f"{self.api_url}/files", fp, full_metadata_params)

      metadatachunk = self._metadata_json_to_df(fp, nan_threshold = nan_threshold)
      sample_metadata_frames += [metadatachunk]

    sample_metadata_frames = pd.concat(sample_metadata_frames)
    metadatabase = metadatabase.join(sample_metadata_frames, how="outer")

    # Include tracker of what is downloaded and not
    metadatabase["downloaded"] = False
    metadatabase.to_csv(self.conf["metadatabase_path"],
                        sep = '\t',
                        compression = "gzip")
    return metadatabase

  # TODO
  #   This can be moved outside of the class, as it does not depend on the
  #   class
  def _metadata_json_to_df(self,
                           json_filepath: str,
                           nan_threshold: float = 1.0) -> pd.DataFrame:
    """
    Summary:
      Converting the full GDC metadata from json format to a Dataframe
    """
    LOG.info("Converting json to tsv for dataframe storage...")
    # First, convert the json to a dataframe
    with open(json_filepath) as f:
      # The field "data.hits" is of most interest
      df = pd.DataFrame(json.loads(f.read())["data"]["hits"])

    # The file UUID as index is necessary to join with metadatabase
    df = df.set_index("id")

	# Expand all of the dictionaries and lists within each column
    df = column_expand(df, "cases")
    df = column_expand(df, "demographic")
    df = column_expand(df, "project")
    df = column_expand(df, "diagnoses")
    df = column_expand(df, "follow_ups")
    df = column_expand(df, "analysis")
    df = column_expand(df, "samples")
    df = column_expand(df, "portions")
    df = column_expand(df, "analytes")
    df = column_expand(df, "aliquots")
    df = column_expand(df, "archive")

    # Convert strings denoting nulls into proper null values
	#   TODO: This can likely be much faster with things like `.fillna()`
    #         or at least not using `applymap()`
    df = df.applymap(lambda x:
                     # Perhaps this could be pd.NA?
                     np.nan if any(s == str(x).lower()
                                   for s in self.conf["null_type_strings"])
                            else x)

    # Drop columns containing NaN data
    #   that is more than the provided threshold (defaults to 100%)
    num_nan_thresh = nan_threshold * len(df)
    df = df.dropna(axis = 1, thresh = num_nan_thresh)

    # Age is in days - change it to years
    if "age_at_diagnosis" in df.columns:
      df["age"] = np.round(df["age_at_diagnosis"]/365, 2)

    LOG.info(f"Converted json database to tsv for of {df.shape[0]} files.")
    return df


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
      params["size"] = "1000000"
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
    manifest = self.metadatabase.loc[self.new_file_ids][["file_name",
                                                         "md5sum",
                                                         "file_size"]]
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
        if assay == "DNAm_450k":
          url = self.conf["DNAm_450k_url"]
          metadata_path = self.conf["DNAm_450k_metadata_filepath"]
        elif assay == "DNAm_27k":
          url = self.conf["DNAm_27k_url"]
          metadata_path = self.conf["DNAm_27k_metadata_filepath"]

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
      download_post(url = f"{self.api_url}/data",
                    data = {"ids": chunk},
                    headers = {"Content-Type": "application/json"})

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

    self.metadatabase["downloaded"].loc[self.new_file_ids] = True
    self.metadatabase.to_csv(self.conf["metadatabase_path"],
                             sep = '\t',
                             compression = "gzip")
    return True

