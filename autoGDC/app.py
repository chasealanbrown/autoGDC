"""The main Dataset object to define a meta-analysis is in this file."""

# General dependencies
import json
import pandas as pd
from tqdm import tqdm

# Specific function imports
from geode import chdir as characteristic_direction

# Local module imports
from .config import SETTINGS, LOG
from .store import Archive
#from .R import wrappers as r_wrappers
from .df_utils import quantile_normalize, combined_region_collapsed_frame


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
               contrasts: list = ["sample_type"],
               paired_assay: bool = False,
			   quantile: bool = True):

    self.conf = SETTINGS[config_key]
    self.params = {"filters" : json.dumps(filt),
                   "format" : "tsv",
                   "fields" : ",".join(fields) if fields is not None else None,
                   "size" : str(size)}

    self.archive = Archive(config_key = config_key,
                           params = self.params,
                           paired_assay = paired_assay)

    self.contrasts = contrasts
    self._data = None

  @property
  def metadata(self):
    return self.archive.metadata


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
    self.archive.metadata = value


  @property
  def data(self):
    """
    Summary:
      A simple property that allows easy access to all
        of the data for this dataset
    """
    if self._data is None:
      self._data = self._get_data()
    return self._data


  def _get_data(self):
    data = {}

    # All file ids in dataframe
    LOG.debug(f"The file_ids for the entire dataset are: {self.archive.file_ids}")

    # File ids already available locally
    LOG.debug(f"Total file_ids already downloaded are: {self.archive.owned_file_ids}")

    # File ids required to download
    LOG.info("Downloading these file_ids:\n{self.archive.new_file_ids}")
    if self.archive.new_file_ids:
       # This one line does the follwing:
       #   1) Downloads all relevant data
       #   2) Organizes the new data
       #        - Moves/deletes files
       #        - Creates/Adds to HDF5 databases
       self.archive._update()

    LOG.info("Querying databases for dataset frames...")
    for assay in tqdm(self.archive.databasefiles):
      LOG.info(f"Querying {assay}...")
      data[assay] = self.archive._read_dataframe(assay)
      # TODO: Figure out a way to store the quantile normed data
      #         - could be just in the HDF5, but needs to have additional
      #           parameter to recalc with more data every so often?
      if self.quantile:
        LOG.info(f"Quantile normalization of {assay}...")
        data[assay] = quantile_normalize(data[assay])

    LOG.info("Constructing paired RNA and DNA methylation dataframe...")
    if self.paired_assay:
      data["RNA_DNAm"] = self._paired_rna_dnam_dataframe()

    return data


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
    meta = self.data["DNAm_450"].columns.to_frame()

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

    filt_methdf = self.data["DNAm_450"]\
                                .reindex(filtered_loci.index).dropna()

    # change position to be bounded between 0 and 1
#    pos = filtered_loci[position_col]
#    filtered_loci[position_col] = (p - p.min()) / p.max()

    return filt_methdf, filtered_loci


  def _paired_rna_dnam_dataframe(self,
                                 collapse_level = "Gene_Symbol",
                                 agg_func = list,
                                 position_col = "Position_to_TSS",
                                 min_seq_len = 15,
                                 max_seq_len = 25,
                                 pos_bounds = (-1500, 1000)):
    """
      This multi-indexed dataframe combines the RNA and DNAm data
    """

    filt_methdf, filtered_loci = self._filter_loci(min_seq = min_seq_len,
                                               max_seq = max_seq_len,
                                               collapse_level = collapse_level,
                                               position_col = position_col,
                                               pos_bounds = pos_bounds)

    dfs_dict = self.data["RNA_counts"]
    df = combined_region_collapsed_frame(dfs_dict = dfs_dict,
                                          main_seq_data = filt_methdf,
                                          seq_feature_metadata = filtered_loci,
                                          region = collapse_level,
                                          agg_func = agg_func,
                                          pos_col = position_col,
                                          max_seq_len = max_seq_len,
                                          value_name = "beta value")
    return df


  def ddx(self,
          contrasts: list = None,
          formula: str = None,
          deseq2: bool = True):
    if contrasts is None:
      contrasts = self.contrasts
    if formula is None:
      formula = "~"+"+".join(contrasts)

    df = self.data["RNA_counts"].astype(int)
    design = self.metadata[contrasts].reindex(df.columns).reset_index()

    if deseq2:
        DEG = r_wrappers.pydeseq(counts = df, design = design, formula = formula)
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

