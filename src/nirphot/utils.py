import os
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


def generate_filename(base_file, extension, output_dir) -> str:
    date_str = datetime.now().strftime("%Y%m%d")
    filename = f"{base_file}_{date_str}.{extension}"

    # create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created the output directory {output_dir}")

    logger.info(f"Generated the filename {filename}")

    return os.path.join(output_dir, filename)
