import logging

from django.core.management.base import BaseCommand


class AbstractCommand(BaseCommand):
    """This class adds common options and methods."""

    def handle(self, *args, **options):
        # Make sure the progress messages will be printed during normal run
        if options['verbosity'] > 0:
            logging.getLogger().setLevel(logging.INFO)

    def print_success(txt):
        """A convenience method to pring messages in nice green color."""
        self.stdout.write(self.style.SUCCESS(txt))

