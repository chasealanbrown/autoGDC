# todo:
#   - Download
#      - 
#   - Collate
#      - 
#   - Store
#      - 
#   - App
#      - 

import os
import re
import pytest
import hashlib
from autoGDC.config import PKGNAME, LOG, CACHE, SETTINGS
from autoGDC.download import Downloader


#ARCHIVE = Arhive()
#ARCHIVE._download_ids(file_ids=["1a8b404b-36c9-462a-9b97-c9a7f027b6b0",
#                                "cc5435a1-2292-4915-9ba8-348140b1c90e"])

md5_re = re.compile(r"([0-9a-f]{32})")


@pytest.fixture
def test_golden_standard():
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
  return EXPECTED_SMALLEST_RNA_DATA


@pytest.fixture
def test_dl_object():
  return Downloader(config_key = "test")


@pytest.fixture
def test_download(test_url):
  test_golden_standard["id"]
  filepath = download_url(url=test_url, dirpath=SETTINGS["test"]["cache_dir"])
  real_md5sum_filepath = download_url(url=f"{test_url}.md5", dirpath=SETTINGS["test"]["cache_dir"])
  with open(real_md5sum_filepath, "rb") as f:
    real_md5sum = md5_re.search(f.read().decode("utf-8")).groups()[0]

  with open(filepath, "rb") as f:
    md5sum = hashlib.md5(f.read()).hexdigest()
  assert real_md5sum == md5sum
  return filepath


@pytest.fixture
def test_download_object(test_url, test_dl_object):
#  ftp_url = test_url.strip("https://")
  test_fpath = test_url.replace(f"https://{test_dl_object.ftp_url}/", "")
  success = test_dl_object._download_remote_file(remote_fpath=test_fpath)
  if success:
    filepath = os.path.join(test_dl_object.conf["data_dir"], test_fpath)

    real_md5sum_filepath = download_url(url=f"{test_url}.md5", dirpath=SETTINGS["test"]["cache_dir"])
    with open(real_md5sum_filepath, "rb") as f:
      real_md5sum = md5_re.search(f.read().decode("utf-8")).groups()[0]

    with open(filepath, "rb") as f:
      md5sum = hashlib.md5(f.read()).hexdigest()
    assert real_md5sum == md5sum
    return test_fpath
  else:
    return None



#def test_diskcache
#def test_diskcache_cleanup
#def test_cleanup
