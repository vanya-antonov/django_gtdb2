import logging
import os
import shutil
import tempfile

from django.conf import settings
from django.contrib.auth.models import User


class GeneTackDB:
    "A class to perform basic tasks with GeneTackDB-2 and associated files."

    root_dir = property(
        fget=lambda self: settings.GTDB_DIR,
        doc="Path to GTDB root dir on the local machine.")

    def __init__(self):
        "Creates an empty list for created tmp dirs (if any)."
        self.all_tmp_dirs = []

    def __del__(self):
        "Make sure to delete all created tmp dirs and files."
        for path in self.all_tmp_dirs:
            shutil.rmtree(path)
            logging.debug("Tmp folder '%s' has been removed." % path)

    def _get_or_create_user(self, username, password='123', last_name=None,
                            email='vanya.antonov@gmail.com'):
        """Returns existing user (by username) or creates a new one."""
        user = User.objects.filter(username=username).first()
        if user is None:
            # https://docs.djangoproject.com/en/2.1/topics/auth/default/#topics-auth-creating-users
            user = User.objects.create_user(
                username=username, password=password, email=email)
            user.last_name = last_name
            user.save()
        return user

    def get_or_create_default_user(self):
        """Returns 'default' user for manage.py commands."""
        return self._get_or_create_user(
            username='gtdb2',
            last_name='Default GeneTackDB user')

    def get_or_create_annotation_user(self):
        """Returns a user for automatic annotation."""
        return self._get_or_create_user(
            username='gtdb2_annotation',
            last_name='Objects created by the automatic annotation')

    def get_full_path_to(self, *args):
        "Returns a full path for a GTDB dir relative path."
        return os.path.join(self.root_dir, *args)

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

    def get_or_create_tmp_dir(self):
        """Returns a full path to the tmp folder inside GTDB dir (and creates
        it, if needed).
        """
        if len(self.all_tmp_dirs) == 0:
            subdir = os.path.join(self.root_dir, 'tmp')
            os.makedirs(subdir, exist_ok=True)
            tmp_path = tempfile.mkdtemp(prefix='__gtdb.', dir=subdir)
            self.all_tmp_dirs.append(tmp_path)
        return self.all_tmp_dirs[0]

