import pandas as pd
import unittest
from autoGDC.app import Dataset
from autoGDC.store import Archive

ARCHIVE = Arhive()

#ARCHIVE._download_ids(file_ids=["1a8b404b-36c9-462a-9b97-c9a7f027b6b0",
#                                "cc5435a1-2292-4915-9ba8-348140b1c90e"])

# Source for checking this output:
# https://portal.gdc.cancer.gov/files/1a8b404b-36c9-462a-9b97-c9a7f027b6b0
EXPECTED_SMALLEST_RNA_DATA = pd.Series(
    {"id":"1a8b404b-36c9-462a-9b97-c9a7f027b6b0",
     "filename":"5a76aa88-3e63-47d8-a982-3e99d439cd9b.mirnaseq.isoforms.quantification.txt",
     "md5":"09691e14d9f18ebaa0afbd4807b36bc6",
     "size":6630,
     "state":"released",
     "acl":"{'open'}",
     "project_id":"TARGET-AML",
     "category":"data_file",
     "data_type":"Isoform Expression Quantification",
     "primary_site":"Hematopoietic and reticuloendothelial systems",
     "disease_type":"Myeloid Leukemias",
     "case_id":"caf1c6a5-48c2-4e4f-8e11-aab7d41fe746",
     "submitter_id":"TARGET-20-PAVBMM",
     "tumor_descriptor":pd.NA,
     "sample_id":"4494cdf1-b909-4a92-949b-09215c0fe288",
     "sample_type":"Primary Blood Derived Cancer - Bone Marrow",
     "tissue_type":"Tumor",
     "aliquot_id":"fff263a6-525c-476c-9900-eba12a7df0f8",
     "gender":"female",
     "race":"white",
     "age_at_diagnosis":473,
     "tumor_stage":pd.NA,
     "age":1.3,
     "cd4_count":pd.NA,
     "risk_factor":pd.NA,
     "bmi":pd.NA,
     "diabetes_treatment_type":pd.NA,
     "downloaded":False,
     "organized":False})


def load_smallest_mirna_study():
  """
  Load the smallest Isoform Expression Quantification sample data
  This should be the same file as `EXPECTED_DATA` above
  """
  # Load all Isoform data with file_size <= 7000 bytes
#  filt = {"op":"and", "content":
#           [{"op":"<=", "content":
#              {"field":"file_size",
#               "value":7000}},
#            {"op":"in", "content":
#              {"field":"data_type",
#               "value":"Isoform Expression Quantification"}}
#           ]}
  filt = {"op":"in":
            "content":
              {"field":"uuid",
               "value":[EXPECTED_SMALLEST_RNA_DATA["id"]]
               }
          }
  study = Dataset(filt = filt)
  return study



def get_small_rna_dnam_file_ids():
  # Load the archive and get the full metadata database
  a = Archive()
  m = a.metadb

  # Filter to DNAm and RNA expression data
  mgm = m[m.data_type.isin(["Gene Expression Quantification",
                            "Methylation Beta Value"])]
  # Which cases have both Gene expression *and* DNAm data?
  s = mgm.groupby("case_id").data_type.nunique()
  cases = s[s==2].index

  # Filter metadatabase to cases which have both Gene expr and DNAm
  mgm = m[m.case_id.isin(cases)]

  # Get the top few cases that have small file sizes for both RNA and DNAm
  sz = mgm.groupby("case_id")["size"].sum()
  sz = sz.sort_values()
  cases = sz.iloc[:2].index

  # Filter to these small file cases that have both DNAm and RNA expr
  mgm = mgm[mgm.case_id.isin(cases)]
  # Get just the counts (not FPKM, FPKM-UQ, etc) and DNAm
  mgm = mgm[mgm.filename.str.contains("(htseq)|(Methylation)")]

  # Now select just the file_ids for these files
  file_ids = mgm.ids.tolist()
  return file_ids


def small_paired_study():
  file_ids = SMALL_PAIRED_FILE_IDS
  filt = {"op":"IN", "content": {"field": "uuid", "value": file_ids}}
  study = Dataset()
  return study


SMALL_PAIRED_FILE_IDS = get_small_rna_dnam_file_ids()
SMALL_PAIRED_STUDY = small_paired_study()
SMALLEST_RNA_STUDY = load_smallest_mirna_study()


class DownloadTest(unittest.TestCase):
  def test_metadata(self):
    # Does the metadata download correctly?
    meta = ARCHIVE._download_metadata(file_ids=["1a8b404b-36c9-462a-9b97-c9a7f027b6b0"])
    self.assertIsInstance(meta, pd.DataFrame)

    # Now, if we download a specific file_id,
    #   does the file_id that we expect to get match?
    #   Also, do the case_id's match what we expect?
    #   Make sure that we don't accidentally have the following situation:
    #   case_id = file_id #!?!
    #   file_id = case_id #!?!
    self.assertTrue() # download 1 file only
    self.assertTrue() # check it matches expected

    m = study.archive.metadb
    m.xs("1a8b404b-36c9-462a-9b97-c9a7f027b6b0", drop_level=False).T
    m[(m.data_type == "Isoform Expression Quantification") & (m["size"]<=7000)].T
    self.assertTrue(study.metadata.shape[0] > 0)

  def test_download(self):
    # Does the RNA/DNAm data download correctly?
    self.assertTrue(ARCHIVE._download(file_ids=SMALL_PAIRED_FILE_IDS))

  def test_gdc_remote_status(self):
    remote_status = ARCHIVE.gdc_status["status"]
    self.assertTrue(remote_status == "OK")

  def test_local_vs_remote_sync(self):
    # Check to make sure the local database is not outdated
    #   (i.e. the local database reflects GDC's remote version state)
    self.assertTrue(ARCHIVE.outdated == False)


class DataDirectoryStructureTest(unittest.TestCase):
  def test_initialize_directories(self):
    self.assertTrue()

class OrganizationTest(unittest.TestCase):


class ArchiveTest(unittest.TestCase):


class StudyTest(unittest.TestCase):
  def test_correct_file_id(self):
    # TODO: This logic is flawed, as a smaller RNA file could get added
    study = SMALLEST_RNA_STUDY
    exptd = EXPECTED_SMALLEST_RNA_DATA
    self.assertTrue(study.metadata.filename[0] == exptd["filename"])


def test_all()
  unittest.main()

if __name__ == '__main__':
  unittest.main()
