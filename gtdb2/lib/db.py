import logging
import os
import shutil

from django.conf import settings


class GeneTackDB:
    """A class to perform basic tasks with GeneTackDB-2 and associated files."""

    def __init__(self):
        """Initialize with a path to gtdb_dir on the local machine."""
        self.gtdb_dir = settings.GTDB_DIR

    def delete_file(self, fn):
        """Argument is the path RELATIVE to the GTDB dir."""
        fn = os.path.join(self.gtdb_dir, fn)
        if os.path.exists(fn):
            os.remove(fn)
        else:
            logging.warning("File '%s' doesn't exist!" % fn)

    def delete_dir(self, path):
        """Argument is the path RELATIVE to the GTDB dir."""
        path = os.path.join(self.gtdb_dir, path)
        if os.path.exists(path):
            shutil.rmtree(path)
        else:
            logging.warning("Folder '%s' doesn't exist!" % path)

