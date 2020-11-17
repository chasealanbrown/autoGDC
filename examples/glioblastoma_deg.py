import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from autoGDC.app import Dataset

plt.rcParams['figure.figsize'] = (20, 5)
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

sns.set_style("whitegrid")


if __name__ == "__main__":
  # This is the type of format that the gdc-client accepts
  #   To discover new combinations easily, try the advanced search on GDC
  #   https://portal.gdc.cancer.gov/query
  filt = {
      "op":'and',
          "content":[
              {"op":"IN",
                  "content":{
                      "field": 'cases.disease_type',
                      "value": ["Adenomas and Adenocarcinomas"]}},
              {"op":"IN",
                  "content":{
                      "field": 'files.analysis.workflow_type',
                      "value": ["HTSeq - Counts",
                                "Liftover"]}},
              {"op":"IN",
                  "content":{
                      "field": "cases.samples.sample_type",
                      "value":  ["primary tumor",
                                 "solid tissue normal"]}}
          ]
  }

  # Main Study Oject
  study = Dataset(filt = filt)#config_key = "default",
                  #filt = filt,
                  #size = 500,
                  #contrasts = ["sample_type"],
                  #paired_assay = True)


  # Paired `primary tumor` and `solid tissue normal` cases
  cases = study.metadata.groupby("case_id")["sample_type"]\
                        .apply(lambda x: len(pd.Series(x).unique())>1)
  study.metadata = study.metadata[study.metadata["case_id"].isin(cases[cases].index.tolist())]

  # Paired cases `HTSeq FPKM` and `Liftover` *WITHIN* the tissue type
  cases = study.metadata.groupby(["case_id", "sample_type"])["workflow_type"]\
                        .apply(lambda x: len(pd.Series(x).unique())>1)\
                        .unstack().sum(axis=1)>1
  study.metadata = study.metadata[study.metadata["case_id"].isin(cases[cases].index.tolist())]
  display(study.metadata.head(3))
  print("Number of files:", study.metadata.shape[0])


