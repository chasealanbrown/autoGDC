import json
import requests
import subprocess
from tqdm import tqdm

# Gene information
import mygene


from .config import SETTING
from .autoGDC import LOG
from .df_utils import column_expand


def download_url(url, filepath, params):
  """
  Summary:
    Download a file from a specific url with a progress bar.
  Source:
    https://stackoverflow.com/questions/56795227/how-do-i-make-progress-bar-while-downloading-file-in-python#56796119
  """
  with requests.get(url, stream=True, params = params) as r:
    r.raise_for_status()
    with open(filepath, 'wb') as f:
      pbar = tqdm(total = int(r.headers['Content-Length']))
      for chunk in r.iter_content(chunk_size = 8192):
        if chunk:  # filter out keep-alive new chunks
          f.write(chunk)
          pbar.update(len(chunk))
  return


class Downloader:
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
  #     {gdc_api_url}/cases
  #     {gdc_api_url}/data
  #     {gdc_api_url}/files
  api_url = "https://api.gdc.cancer.gov"

  def __init__(self,
               config_key: str = "default",
               params: dict = None,
               file_ids: list = []):

    self.conf = SETTING[config_key]
    self.params = params
    self._file_ids = None
    self._metadata = None
    self._manifest = None
    self._metadatabase = None

    self.mg = mygene.MyGeneInfo()
    self.mg.set_caching(cache_db = self.conf["mygene_dir"])


  @property
  def file_ids(self):
    """
    Summary:
      The file_ids to be downloaded from GDC
    """
    if self._file_ids is None:
      self._file_ids = self._get_file_ids()
    return self._file_ids


  @property
  def metadata(self):
    """
    Summary:
      The table of GDC metadata for this specific download
    """
    if self._metadata is None:
      self._metadata = self.metadatabase.loc[self.file_ids]
    return self._metadata


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
      return pd.read_csv(self.conf["metadatabase_path"])
    else:
      # First we get the table provided by GDC
      #   of filenames, UUIDs, md5sum, etc.
      LOG.info("Downloading list of all files in the GDC...")
      download_url(self.conf["gdc_filelist_url"],
                   self.conf["metadatabase_path"])
      metadatabase = pd.read_csv(self.conf["metadatabase_path"], gzip = True)

      # Now download clinical and demographic metadata for each file
      LOG.info("Downloading metadata for all files in the GDC...")

      # This must be done in chunks (server errors arise if larger sizes used)
      chunk_size = 100000
      for pos in tqdm(range(0, len(metadatabase), chunk_size)):

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
        fp = path.join(self.conf["data_dir"], "metadata{pos}-{to}.json")
        download_url(f"{self.api_url}/files", fp, full_metadata_params)

        metadatachunk = self._metadata_json_to_tsv(fp, nan_threshold = nan_threshold)
        metadatabase.join(metadatachunk)

      metadatabase.to_csv(self.conf["metadatabase_path"])
      return metadatabase

  def _metadata_json_to_tsv(self,
                            json_filepath: str,
                            nan_threshold: float = 1.0) -> pd.DataFrame:
    """
    Summary:
      Converting the full GDC metadata from json format to a Dataframe
    """
    # First, convert the json to a dataframe
    with open(json_filepath) as f:
      # The field "data.hits" is of most interest
      df = pd.DataFrame(json.loads(f.read())["data"]["hits"])

	# Expand all of the dictionaries and lists within each column
    df = df_utils.column_expand(df, "cases")
    df = df_utils.column_expand(df, "demographic")
    df = df_utils.column_expand(df, "project")
    df = df_utils.column_expand(df, "diagnoses")
    df = df_utils.column_expand(df, "follow_ups")
    df = df_utils.column_expand(df, "samples")
    df = df_utils.column_expand(df, "portions")
    df = df_utils.column_expand(df, "analytes")
    df = df_utils.column_expand(df, "aliquots")

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

    # The file UUID as index
    #   for comparison when determining download requirements
    df = df.set_index("id")

    LOG.info(f"Downloaded df database of {df.shape[0]} files.")
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

    response = requests.get(f"{self.api_url}/files", params = params)
    response_file = io.StringIO(response.content.decode("utf-8"))
    try:
      response_df = pd.read_csv(response_file, sep = "\t")
      LOG.debug(f"""Query for GDC file_ids returned the following dataframe:
                \n{response_df}""")
      file_ids = response_df["ids"].tolist()
    except EmptyDataError as e:
      LOG.warn(f"It appears that there were no samples that fit the criteria.")
      LOG.exception(e)
      file_ids = []
    return file_ids


  def _get_manifest(self):
    """
    Summary:
      Retrieves the manifest table to store as a property for this download.
    """
    # Create *complete* manifest file for downloading the files
    boolmask = self.metadatabase.ids.isin(self.file_ids)
    manifest = self.metadatabase.loc[self.file_ids][["file_name",
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
    LOG.info(f"Starting download of {manifest.shape[0]} files")
    for position in tqdm(range(0, len(self.manifest), chunk_size)):
      # Debugging
      pct_dn = round( position/len(self.manifest) * 100.0, 2)
      msg = f"Downloading chunk of {chunk_size} files... ({pct_dn}% finished)"
      LOG.debug(msg)

      chunk = self.manifest.iloc[position: position + chunk_size]

      temp_manifest_file = path.join(self.config["data_dir"],
                                     "temporary_manifest.txt")
      chunk.to_csv(temp_manifest_file, sep = "\t")
      subprocess.call([gdc_client_path,
                       "download",
                       "-m",
                       temp_manifest_file,
                       "--dir",
                       self.config["rawdata_dir"]])

    LOG.info("All files have been downloaded.")
    self._download_feature_metadata()

    return True

