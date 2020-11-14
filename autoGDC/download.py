

class Downloader:
  # The GDC api url is used to obtain data via the REST API
  #   Useful urls are:
  #     {gdc_api_url}/cases
  #     {gdc_api_url}/data
  #     {gdc_api_url}/files
  api_url = "https://api.gdc.cancer.gov"

  def __init__(self,
               data_dir: str = None,
               keep_raw: bool = False):
    self.archive = Archive(self.data_dir)


