"""
Only for use on Travis. Sets up the database for testing
"""
import fiasco.util
paths = fiasco.util.setup_paths()
fiasco.util.download_dbase(paths['ascii_dbase_root'], ask_before=False)
fiasco.util.build_hdf5_dbase(paths['ascii_dbase_root'], paths['hdf5_dbase_root'], ask_before=False)