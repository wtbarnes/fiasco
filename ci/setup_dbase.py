"""
Only for use on Travis. Sets up the database for testing
"""
import fiasco.util
from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION

paths = fiasco.util.setup_paths()
fiasco.util.download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), paths['ascii_dbase_root'])
fiasco.util.build_hdf5_dbase(paths['ascii_dbase_root'], paths['hdf5_dbase_root'])
