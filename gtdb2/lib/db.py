import logging
import os
import shutil

from django.conf import settings

from gtdb2.models.user import User


class GeneTackDB:
    """A class to perform basic tasks with GeneTackDB-2 and associated files."""

    def __init__(self):
        """Initialize with a path to GTDB root dir on the local machine."""
        self.root_dir = settings.GTDB_DIR

    def get_default_user(self):
        """Returns 'default' user for manage.py commands."""
        user, created = User.objects.get_or_create(
            name='gtdb_manage.py',
            descr='Default GeneTackDB user')
        return user

    def delete_file(self, fn):
        """Argument is the path RELATIVE to the GTDB dir."""
        fn = os.path.join(self.root_dir, fn)
        if os.path.exists(fn):
            os.remove(fn)
        else:
            logging.warning("File '%s' doesn't exist!" % fn)

    def delete_dir(self, path):
        """Argument is the path RELATIVE to the GTDB dir."""
        path = os.path.join(self.root_dir, path)
        if os.path.exists(path):
            shutil.rmtree(path)
        else:
            logging.warning("Folder '%s' doesn't exist!" % path)

