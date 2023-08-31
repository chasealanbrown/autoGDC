import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import common_studies

plt.rcParams['figure.figsize'] = (20, 5)
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

sns.set_style("whitegrid")

if __name__ == "__main__":

  study = common_studies.paired_dnam_rna()
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

