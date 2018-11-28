from django.conf import settings
from django.core.management.base import BaseCommand


class AbstractCommand(BaseCommand):
    """This class adds common options and methods."""

    def add_arguments(self, parser):
        """Add --user option to all commands."""
        name = settings.GTDB_USER
        parser.add_argument('--user', metavar='name', default=name,
                            help='[%s] user to create new database units' % name)

    def print_success(txt):
        """A convenience method to pring messages in nice green color."""
        self.stdout.write(self.style.SUCCESS(txt))

