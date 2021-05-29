import unittest
from autoGDC.app import Dataset

study = Dataset()

class DownloadTest(unittest.TestCase):
  def test_metadata(self):
    self.assertTrue(study.metadata.shape[0] > 0)

  def test_download(self):
    file_ids = study.metadata.index[0:2].tolist()
    self.assertTrue(study._download(file_ids))

class DataDirectoryStructureTest(unittest.TestCase):
  def test_initialize_directories(self):
    self.assertTrue()


def test_all()
  unittest.main()

if __name__ == '__main__':
  unittest.main()
