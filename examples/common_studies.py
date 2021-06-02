from autoGDC.app import Dataset
from autoGDC.store import Archive


def all_GBM_RNA_data(workflow_type:str = "HTSeq - FPKM-UQ"):
  """
  Helper function to grab all samples labeled as GBM RNA expression data
  """
  assert workflow_type in ["HTSeq - FPKM-UQ",
                           "HTSeq - FPKM",
                           "HTSeq - Counts",
                           "STAR - Counts"]

  # Using metadata database on local disk
  study = Dataset(config_key = "default")
  archive = study.archive
#  archive = Archive()
  metadb = archive.metadb
  expr_bool = metadb.data_type == "Gene Expression Quantification"
  disease_bool = metadb.disease_type == "Gliomas"
  #  workflow_bool = metadb.workflow_type == workflow_type

  # Should be ~900 samples
  gbm_expr = metadb[expr_bool & disease_bool] # & workflow_bool]

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
  study.filt = filt
  study.size = 10**6
  study.contrasts = ["sample_type"]
  study.paired_assay = True
#  study = Dataset(config_key = "default",
#                  filt = filt,
#                  size = 10**6,
#                  contrasts = ["sample_type"],
#                  paired_assay = True)
  return study


def paired_dnam_rna(diseases=None, size = 1000):
  if diseases is not None:
    disease_query = {"op":"IN",
                      "content":{"field": 'cases.disease_type',
                                 "value": diseases}} # e.g ["Gliomas"]
  else:
    disease_query = None
  dnam_rna_query = {"op":"IN",
                    "content":
                      {"field": 'data_type',
                       "value": [
                          "Gene Expression Quantification", # ~43k samples
                          "miRNA Expression Quantification", # ~16k samples
                          "Isoform Expression Quantification", # ~16k samples
                          "Methylation Beta Value"] # ~12k samples
                       }}
  if disease_query:
    comb_queries = [dnam_rna_query, disease_query]
  else:
    comb_queries = [dnam_rna_query]

  filt = {"op":'and', "content": comb_queries}

  # Main Study Oject
  study = Dataset(config_key = "default",
                  filt = filt,
                  size = size,
                  contrasts = ["sample_type"],
                  paired_assay = True)
  return study

