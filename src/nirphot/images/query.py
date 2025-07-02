import logging
from astroquery.mast import Observations
from astropy.table import unique, vstack

from nirphot.paths import DATA_DIR

logger = logging.getLogger(__name__)


class MastQuery:
    def __init__(self, auth_id=None, proposal_id=None):
        self.auth_id = auth_id
        self.proposal_id = proposal_id
        self.download_dir = DATA_DIR

        self.obs_collection = "JWST"
        self.instrument_name = "NIRCAM/IMAGE"

    def authenticate(self):
        logger.info(f"Authenticating with MAST using {self.auth_id}")

        # don't auth if None
        if self.auth_id is None:
            logger.info("No authentication ID provided, skipping authentication")
            return

        try:
            Observations.login(self.auth_id)
        except Exception as e:
            logger.error(f"Failed to authenticate with MAST: {str(e)}")
            raise Exception(f"Failed to authenticate with MAST: {str(e)}")

        return

    def query(self, columns=[]):
        if not columns:
            columns = [
                "dataproduct_type",
                "filters",
                "calib_level",
                "t_exptime",
                "proposal_pi",
                "intentType",
                "obsid",
                "instrument_name",
            ]

        logger.info(f"Querying MAST for PID: {self.proposal_id}")

        matched_obs = Observations.query_criteria(
            obs_collection=self.obs_collection,
            proposal_id=self.proposal_id,
            instrument_name=self.instrument_name,
        )

        return matched_obs

    def retrive_prods(self, size_chunk=5):
        matched_obs = self.query()

        chunks = [
            matched_obs[i : i + size_chunk]
            for i in range(0, len(matched_obs), size_chunk)
        ]

        prod_lists = [Observations.get_product_list(chunk) for chunk in chunks]

        files = unique(vstack(prod_lists), keys="productFilename")

        logger.info(f"Retrieved {len(files)} files")
        logger.info(f"Size of files: {sum(files['size'] / 1e9):.1f} GB")

        return files

    def download(self, curl_flag=True, **kwargs):
        # https://mast.stsci.edu/api/v0/_productsfields.html
        self.authenticate()

        files = self.retrive_prods()

        Observations.download_products(
            files, download_dir=self.download_dir, curl_flag=curl_flag, **kwargs
        )
