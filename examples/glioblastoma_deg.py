import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from autoGDC.app import Dataset
from autoGDC.store import Archive

plt.rcParams['figure.figsize'] = (20, 5)
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

sns.set_style("whitegrid")

def all_GBM_RNA_data(workflow_type:str = "HTSeq - FPKM-UQ"):
  """
  Helper function to grab all samples labeled as GBM RNA expression data
  """
  assert workflow_type in ["HTSeq - FPKM-UQ",
                           "HTSeq - FPKM",
                           "HTSeq - Counts",
                           "STAR - Counts"]

  # Using metadata database on local disk
  archive = Archive()
  metadb = archive.metadatabase
  expr_bool = metadb.data_type == "Gene Expression Quantification"
  disease_bool = metadb.disease_type == "Gliomas"
  workflow_bool = metadb.workflow_type == workflow_type

  # Should be ~900 samples
  gbm_expr = metadb[expr_bool & disease_bool & workflow_bool]

  # This is the type of format that the gdc-client accepts
  #   To discover new combinations easily, try the advanced search on GDC
  #   https://portal.gdc.cancer.gov/query
  filt = { "op":'and', "content":[
              {"op":"IN",
                  "content":{ "field": 'cases.disease_type',
                      "value": ["Gliomas"]}},
              {"op":"IN",
                  "content":{ "field": 'data_type',
                      "value": ["Gene Expression Quantification"]}} ] }

  # Main Study Oject
  study = Dataset(config_key = "default",
                  filt = filt,
                  size = 10**6,
                  contrasts = ["sample_type"],
                  paired_assay = True)
  return study

if __name__ == "__main__":

  study = all_GBM_RNA_data()
  study.ddx()

#  # Paired `primary tumor` and `solid tissue normal` cases
#  cases = study.metadata.groupby("case_id")["sample_type"]\
#                        .apply(lambda x: len(pd.Series(x).unique())>1)
#  study.metadata = study.metadata[study.metadata["case_id"].isin(cases[cases].index.tolist())]
#
#  # Paired cases `HTSeq FPKM` and `Liftover` *WITHIN* the tissue type
#  cases = study.metadata.groupby(["case_id", "sample_type"])["workflow_type"]\
#                        .apply(lambda x: len(pd.Series(x).unique())>1)\
#                        .unstack().sum(axis=1)>1
#  study.metadata = study.metadata[study.metadata["case_id"].isin(cases[cases].index.tolist())]
  print(study.metadata.head(3))
  print("Number of files:", study.metadata.shape[0])

