import pytest
import os
from autoGDC.config import PKGNAME, LOG, CACHE, SETTINGS


def test_wrenches():
  assert PKGNAME == "autoGDC"
  assert LOG.name == "autoGDC"
  # TODO: How to test cache directory?

  # Ensure read/write capability
  CACHE.set("test", 1)
  assert CACHE.get("test")

def test_filesystem():
  tdp = str(SETTINGS["test"]["data_dir"])
  tcp = str(SETTINGS["test"]["cache_dir"])
  ddp = str(SETTINGS["default"]["data_dir"])
  dcp = str(SETTINGS["default"]["cache_dir"])
  # Test existance of testing and default data and cache directories
  assert os.path.exists(tdp)
  assert os.path.exists(tcp)
  assert os.path.exists(ddp)
  assert os.path.exists(dcp)

  # Test reading of testing and default data and cache directories
  assert os.access(tdp, os.R_OK)
  assert os.access(tcp, os.R_OK)
  assert os.access(tdp, os.R_OK)
  assert os.access(tcp, os.R_OK)

  # Test writing of testing and default data and cache directories
  assert os.access(tdp, os.W_OK)
  assert os.access(tcp, os.W_OK)
  assert os.access(tdp, os.W_OK)
  assert os.access(tcp, os.W_OK)



#def test_config():
#  conf = "test"
#  assert val

