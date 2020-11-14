from os import path, listdir


class Archive:
  """
  Summary:
    Class that keeps track of the inventory of downloaded files

  Arguments:
    data_dir:
      The main directory where all of the downloaded GDC files will be stored
  """

  def __init__(self,
               data_dir: str = None,
               keep_raw: bool = False):
    if data_dir is None:
      self.data_dir = default_data_dir
    else:
      self.data_dir = data_dir

    self.filetype_regexs = {"DNAm_450":"HumanMethylation450",
                            "DNAm_27":"HumanMethylation27",
                            "RNA_FPKM":"FPKM.txt",
                            "RNA_FPKM-UQ":"FPKM-UQ",
                            # Files can have either '-' or '.' in name
                            "RNA_counts":"htseq[\.\-]counts.gz",
                            "RNA_miRNA":"mirnas.quantification",
                            "RNA_isoforms":"isoforms.quantification"
                           }


  def _read_dataframe(self, assay_filepath):
    data = pd.read_h5(assay_filepath, "data")
    sample_metadata = pd.read_h5(assay_filepath, "sample_metadata")
    feature_metadata = pd.read_h5(assay_filepath, "feature_metadata")
    return pd.DataFrame(data.values,
                        index = feature_metadata,
                        columns = sample_metadata)

  @property
  def dataframe_files(self):
    """
    Returns:
      A dictionary containing the file paths of the dataframes
        (combined gdc files with the same assay type)
        (indexed by assay type keys)
    """
    file_map = {assay:path.join(self.data_dir, f"{assay}.h5")
                for assay in self.filetype_regexs.keys()}
    return file_map


  @property
  def dataframes(self):
    """
    Summary:
      A simple property that allows easy access to all data that has been
        downloaded
    """
    dataframe_map = {assay:pd.read_h5(fp, "data")
                     for assay,filepath in self.dataframe_files.items()}
    return

  @property
  def files(self):
    """
    Returns:
      A dictionary containing the lists of files
        indexed by filetype (directory) keys that have been downloaded
    """
    file_map = {}

    # First, check the directories for actual raw files
    assay_types = self.filetype_regexs.keys()
    for assay_type in assay_types:
      assay_dir = self.data_dirs[assaytype]
      file_map[assay_dir] = [path.join(assay_dir, f)
                             for f in listdir(assay_dir)]

    # Then add the files that are on the unprocessed dataframes
    #   for each data type
    return


  @property
  def uuids(self):
    """
    Returns:
      A dictionary containing the lists of uuids
        indexed by filetype (directory) keys that have been downloaded
    """
    uuid_map = {}

    # First, check the directories for actual raw files
    assay_types = self.filetype_regexs.keys()
    for assay_type in assay_types:
      assay_dir = self.data_dirs[assaytype]
      uuid_map[assay_type] = [basename(splitext(]

    for assaydir in file_map
    for filename in listdir(
    # Then add the files that are on the unprocessed dataframes
    #   for each data type
    return

  @property
  def data_dirs(self):
    if self._all_data_dirs is None:
      # Create / assign download directory
      #   Structured to have a directory for each assay type
      self._all_data_dirs = {}

      # Main directory
      self.data_dir = path.join(path.expanduser(self.data_dir))
      os.makedirs(self.data_dir, exist_ok = True)
      self._all_data_dirs["main"] = self.data_dir

      # Directory for each biological assay
      for assaytype in self.file_types.keys():
        assaytype_dir = path.join(self.data_dir, assaytype)
        os.makedirs(assaytype_dir, exist_ok = True)
        self._all_data_dirs[assaytype] = assaytype_dir

      return self._all_data_dirs
    else:
      return self._all_data_dirs

