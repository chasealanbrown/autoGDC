import requests
from tqdm import tqdm

from .config import SETTING
from .autoGDC import LOG

def download_url(url, filepath):
  with requests.get(url, stream=True) as r:
    r.raise_for_status()
    with open(filepath, 'wb') as f:
      pbar = tqdm(total = int(r.headers['Content-Length']))
      for chunk in r.iter_content(chunk_size = 8192):
        if chunk:  # filter out keep-alive new chunks
          f.write(chunk)
          pbar.update(len(chunk))
  return


class Downloader:
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
    self.file_ids = file_ids
    self._metadata = None


  @property
  def metadata(self):
    if self._metadata is None:
      self._metadata = self._get_metadata()
    return self._metadata


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
    return self._clean_metadata(metadata)


  def _clean_metadata(self,
                      metadata: pd.DataFrame,
                      nan_threshold: float = 0.9):

    # Simplify the column names
    #   Example: cases.0.follow_ups.0.bmi -> bmi
    #   TODO: If both cases.0.follow_ups.1.bmi and cases.0.follow_ups.0.bmi
    #         are present, how should this be handled?
    metadata.columns = [col.split(".")[-1] for col in metadata.columns]

    # Convert strings denoting nulls into proper null values
    metadata = metadata.applymap(lambda x:
                                 # Perhaps this could be pd.NA?
                                 np.nan if any(s == str(x).lower()
                                           for s in self.conf["null_type_strings"])
                                 else x)

    # Drop columns containing NaN data
    #   that is more than the provided threshold (defaults to 90%)
    num_nan_thresh = nan_threshold * len(metadata)
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


  def _feature_metadata_check(self):
    LOG.debug("Checking if DNA methylation feature metadata exists...")

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


  @property
  def manifest(self):
    # Create *complete* manifest file for downloading the files
    manifest = self.metadata.loc[file_ids][["file_name", "md5sum", "file_size"]]
    manifest["state"] = "validated"
    manifest.columns = ["filename", "md5", "size", "state"]



  def _download(self, file_ids):
    LOG.info("Downloading data... (This may hang and take a while!)")

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


